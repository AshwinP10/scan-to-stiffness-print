"""Scan + stiffness-labelled points -> infill-graded print file (proof of concept).

  python scan_to_infill.py --dummy-cube                     # fast sanity check
  python scan_to_infill.py --scan part.glb --n-points 12    # Polycam scan, random stiffness points
  python scan_to_infill.py --scan part.glb --points gantry.csv   # real probe data (x,y,z,stiffness)

Stage A: load + repair + register scan, stiffness point records
Stage B: voxel labelling (column or volumetric Voronoi), region meshes, stiffness -> infill, visualisation
Stage C: 3MF with one infill override per region (modifiers, or split parts above 12 regions)
         + PrusaSlicer 2.5.0 CLI (pinned TPU profile), with a load/re-export check
Stage D: parse sliced G-code, realised infill per region vs target
"""
import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pymeshfix
import shapely
import trimesh
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree
from shapely.geometry import MultiLineString, MultiPoint, Point, Polygon, box
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parent
CMAP = plt.get_cmap("viridis")

K_RANGE = (0.5, 3.0)          # N/mm, placeholder probe range
INFILL_RANGE = (10.0, 90.0)   # %, kept below 100 so slicers keep it as sparse infill

# Fast sanity check: a 40 x 40 x 20 mm cube probed from the top, plus one side contact
# so the partition is genuinely 3D rather than a stack of identical columns.
DUMMY_CUBE_SIZE = (40.0, 40.0, 20.0)
DUMMY_POINTS = [
    (8.0, 8.0, 20.0, 0.6),
    (32.0, 8.0, 20.0, 2.8),
    (20.0, 20.0, 20.0, 1.6),
    (8.0, 32.0, 20.0, 2.2),
    (32.0, 32.0, 20.0, 1.0),
    (40.0, 20.0, 6.0, 2.5),
]


# ============================================================================ STAGE A
def load_scan(path, units="auto", scale=1.0):
    """Load a .glb/.gltf (or any trimesh-readable mesh) into slicer conventions: Z up, millimetres."""
    textured = trimesh.load(path, force="mesh")
    # Geometry only: textured scans duplicate vertices along UV seams, which leaves the surface in
    # thousands of disconnected patches until those duplicates are merged.
    mesh = trimesh.Trimesh(textured.vertices, textured.faces, process=True)
    if Path(path).suffix.lower() in (".glb", ".gltf"):
        # glTF is Y-up; slicers are Z-up.
        mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0]))
    if units == "auto":
        # Polycam exports metres; anything under 5 units across is assumed to be metres.
        scale_units = 1000.0 if float(mesh.extents.max()) < 5.0 else 1.0
    else:
        scale_units = {"m": 1000.0, "mm": 1.0}[units]
    mesh.apply_scale(scale_units * scale)
    print(f"  loaded {path}: {len(mesh.faces)} faces, units x{scale_units:g}, print scale x{scale:g} "
          f"-> extents {np.round(mesh.extents, 1)} mm")
    return mesh


def repair_mesh(mesh):
    """Return a closed, watertight, outward-facing copy of mesh or raise."""
    m = mesh.copy()
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.remove_unreferenced_vertices()
    trimesh.repair.fill_holes(m)
    trimesh.repair.fix_normals(m)
    if not m.is_watertight:
        open_edges = len(trimesh.grouping.group_rows(m.edges_sorted, require_count=1))
        print(f"  trimesh cleanup not watertight ({open_edges} open edges); running MeshFix")
        v, f = pymeshfix.clean_from_arrays(np.asarray(m.vertices), np.asarray(m.faces))
        m = trimesh.Trimesh(v, f, process=True)
        trimesh.repair.fix_normals(m)
    if m.volume < 0:
        m.invert()
    if not (m.is_watertight and m.is_winding_consistent and m.volume > 0):
        raise RuntimeError(f"Mesh repair failed: watertight={m.is_watertight}, "
                           f"winding_consistent={m.is_winding_consistent}, volume={m.volume:.3f}")
    print(f"  repaired: watertight, {len(m.faces)} faces, volume {m.volume:.1f} mm^3")
    return m


def registration_transform():
    """Scan frame -> gantry/turntable frame. Identity placeholder.

    Replace with the transform solved from the turntable origin + calibration point.
    """
    return np.eye(4)


def to_part_frame(mesh, records):
    """Shift mesh and points together so the part's bounding box starts at the origin (z=0 on the bed)."""
    shift = -mesh.bounds[0]
    mesh = mesh.copy()
    mesh.apply_translation(shift)
    moved = [{**r, "x": r["x"] + shift[0], "y": r["y"] + shift[1], "z": r["z"] + shift[2]} for r in records]
    return mesh, moved, shift


def make_point_records(xyz, stiffness):
    return [{"id": i, "x": float(p[0]), "y": float(p[1]), "z": float(p[2]), "stiffness": float(k)}
            for i, (p, k) in enumerate(zip(xyz, stiffness))]


def sample_stiffness_points(mesh, n, rng, k_range=K_RANGE, min_spacing=8.0, top_only=True, inset=0.0,
                            reach_radius=1.0, reach_tol=0.5):
    """Random surface contacts with random stiffness, in the same record format real gantry data will use.

    top_only keeps only surface a top-down probe can reach: upward-facing, and with no higher surface
    within reach_radius (about the probe tip radius) in x/y, which excludes ledges under overhangs.
    """
    pts, face_idx = trimesh.sample.sample_surface(mesh, max(200 * n, 5000), seed=int(rng.integers(1 << 31)))
    normals = mesh.face_normals[face_idx]
    if top_only:
        keep = normals[:, 2] > 0.7
        pts, normals = pts[keep], normals[keep]
        dense, _ = trimesh.sample.sample_surface(mesh, 200000, seed=int(rng.integers(1 << 31)))
        tree = cKDTree(dense[:, :2])
        local_top = np.array([dense[nb, 2].max() if nb else -np.inf
                              for nb in tree.query_ball_point(pts[:, :2], reach_radius)])
        reachable = pts[:, 2] >= local_top - reach_tol
        pts, normals = pts[reachable], normals[reachable]
    pts = pts - inset * normals
    chosen = []
    for i in rng.permutation(len(pts)):
        if all(np.linalg.norm(pts[i] - pts[j]) >= min_spacing for j in chosen):
            chosen.append(i)
            if len(chosen) == n:
                break
    if len(chosen) < n:
        raise ValueError(f"Placed {len(chosen)}/{n} points at min_spacing={min_spacing} mm; reduce one of them.")
    return make_point_records(pts[chosen], rng.uniform(k_range[0], k_range[1], n))


def save_points_csv(records, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "x", "y", "z", "stiffness"])
        w.writeheader()
        w.writerows(records)


def load_points_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [{"id": i, "x": float(r["x"]), "y": float(r["y"]), "z": float(r["z"]), "stiffness": float(r["stiffness"])}
            for i, r in enumerate(rows)]


def points_xyz(records):
    return np.array([[r["x"], r["y"], r["z"]] for r in records])


# ============================================================================ STAGE B
def stiffness_to_infill(k, k_range=K_RANGE, infill_range=INFILL_RANGE):
    """Placeholder linear map, stiffness (N/mm) -> infill (%). Swap in the measured calibration curve."""
    t = np.clip((np.asarray(k, dtype=float) - k_range[0]) / (k_range[1] - k_range[0]), 0.0, 1.0)
    return infill_range[0] + t * (infill_range[1] - infill_range[0])


def voxelize(mesh, pitch, pad=2):
    """Filled voxel grid of the mesh, embedded in a grid padded by `pad` voxels on every side."""
    vg = mesh.voxelized(pitch=pitch).fill()
    occupied = np.pad(vg.matrix, pad)
    origin = vg.transform[:3, 3] - pad * pitch  # world centre of padded voxel (0, 0, 0)
    print(f"  voxelised at {pitch:g} mm: {int(vg.matrix.sum())} occupied voxels, grid {occupied.shape}")
    return {"occupied": occupied, "origin": origin, "pitch": float(pitch)}


def voxel_centers(grid):
    idx = np.indices(grid["occupied"].shape).reshape(3, -1).T
    return grid["origin"] + idx * grid["pitch"]


def label_voxels(grid, xyz):
    """Nearest stiffness point (3D) for every voxel of the padded grid.

    Labelling the empty padding too (not just occupied voxels) lets each region mesh extend a couple of
    voxels past the scan surface, so the slicer clips modifiers to the true surface instead of leaving
    an unmodified skin where voxels only partly overlap the part.
    """
    tree = cKDTree(xyz)
    _, labels = tree.query(voxel_centers(grid))
    return labels.reshape(grid["occupied"].shape)


def top_contact_mask(grid, xyz):
    """True for points on the highest surface of their (x, y) column, i.e. reachable by a top-down probe."""
    occ, pitch, origin = grid["occupied"], grid["pitch"], grid["origin"]
    ij = np.clip(np.round((xyz[:, :2] - origin[:2]) / pitch).astype(int), 0, np.array(occ.shape[:2]) - 1)
    col = occ[ij[:, 0], ij[:, 1], :]
    top_z = origin[2] + (occ.shape[2] - 1 - np.argmax(col[:, ::-1], axis=1)) * pitch
    return col.any(axis=1) & (np.abs(xyz[:, 2] - top_z) <= 2.0 * pitch)


def partition_points(xyz, mode, top_mask):
    """Coordinates the partition is computed from.

    columns: every point governs the vertical column beneath it (points flattened to one height), so a
    reading taken on the top surface never decides the infill of an unrelated region deep below it.
    3d: true volumetric Voronoi, needed once side contacts exist. auto: columns if every point is a
    top contact, otherwise 3d.
    """
    if mode == "auto":
        mode = "columns" if top_mask.all() else "3d"
    p = np.asarray(xyz, dtype=float).copy()
    if mode == "columns":
        p[:, 2] = 0.0
    return p, mode


def mask_to_mesh(mask, grid):
    """Closed surface of a boolean voxel mask (marching cubes at the 0.5 iso-level), in world mm."""
    padded = np.pad(mask.astype(np.float32), 1)
    verts, faces, _, _ = marching_cubes(padded, level=0.5)
    verts = grid["origin"] + (verts - 1.0) * grid["pitch"]
    m = trimesh.Trimesh(verts, faces, process=True)
    if m.volume < 0:
        m.invert()
    return m


def region_meshes(grid, labels, n_regions):
    """Per label: a modifier mesh (label over the padded grid) and a display mesh (label inside the part).

    Adjacent labels partition the grid exactly, so their iso-surfaces coincide: no gaps or overlaps.
    """
    modifiers, display = [], []
    for i in range(n_regions):
        sel = labels == i
        mod = mask_to_mesh(sel, grid)
        if not mod.is_watertight or mod.volume <= 0:
            raise RuntimeError(f"Region {i} modifier mesh is not a closed solid")
        modifiers.append(mod)
        inside = sel & grid["occupied"]
        display.append(mask_to_mesh(inside, grid) if inside.any() else None)
    return modifiers, display


def prism_mesh(poly, z0, z1):
    """Watertight extrusion of a convex polygon between z0 and z1, outward-facing normals."""
    ring = np.asarray(shapely.geometry.polygon.orient(poly, sign=1.0).exterior.coords[:-1])
    n = len(ring)
    verts = np.vstack([np.column_stack([ring, np.full(n, z0)]), np.column_stack([ring, np.full(n, z1)])])
    faces = [[0, i + 1, i] for i in range(1, n - 1)] + [[n, n + i, n + i + 1] for i in range(1, n - 1)]
    for i in range(n):
        j = (i + 1) % n
        faces += [[i, j, n + j], [i, n + j, n + i]]
    m = trimesh.Trimesh(verts, np.array(faces), process=False)
    if not m.is_watertight or m.volume <= 0:
        raise RuntimeError("Column modifier prism is not a closed solid")
    return m


def column_modifiers(mesh, part_xyz, overhang=2.0, seam_overlap=0.05):
    """Exact 2D Voronoi cells of the probe points' x/y, extruded through the part height.

    Each cell is convex, so the modifier is a handful of planar faces rather than a voxel staircase.
    Cells overhang the part outline (the slicer clips modifiers to the part) and overlap their
    neighbours by seam_overlap so no sliver on a shared face falls back to the base infill.
    """
    (x0, y0, z0), (x1, y1, z1) = mesh.bounds
    frame = box(x0 - overhang, y0 - overhang, x1 + overhang, y1 + overhang)
    xy = part_xyz[:, :2]
    cells = shapely.voronoi_polygons(MultiPoint(xy), extend_to=frame, ordered=True)
    cells2d, mods = [], []
    for p, cell in zip(xy, cells.geoms):
        c = cell.intersection(frame).simplify(1e-9)
        if not c.contains(Point(p)):
            raise RuntimeError("Voronoi cell order does not match the probe points")
        cells2d.append(c)
        grown = c.buffer(seam_overlap, join_style="mitre").intersection(frame).simplify(1e-9)
        mods.append(prism_mesh(grown, z0, z1))
    return cells2d, mods


def make_region_lookup(mode, part_xyz, grid, labels):
    """Region index for arbitrary part-frame points, matching the modifier geometry of each mode."""
    if mode == "columns":
        tree = cKDTree(part_xyz[:, :2])
        return lambda pts: tree.query(np.asarray(pts)[:, :2])[1]
    return lambda pts: voxel_label_at(grid, labels, pts)


def plot_regions_3d(grid, labels, infill, xyz, path, mode="3d", max_points=60000, seed=0):
    """Isometric shell of the labelled part, a top-down map, and a vertical slice through the part centre."""
    occ, pitch, origin = grid["occupied"], grid["pitch"], grid["origin"]
    infill = np.asarray(infill)
    norm = Normalize(0, 100)
    shell = occ & ~(np.roll(occ, 1, 0) & np.roll(occ, -1, 0) & np.roll(occ, 1, 1) & np.roll(occ, -1, 1)
                    & np.roll(occ, 1, 2) & np.roll(occ, -1, 2))
    centers = voxel_centers(grid)[shell.ravel()]
    vals = infill[labels[shell]]
    if len(centers) > max_points:
        keep = np.random.default_rng(seed).choice(len(centers), max_points, replace=False)
        centers, vals = centers[keep], vals[keep]

    fig = plt.figure(figsize=(19, 6.2), dpi=170)
    fig.subplots_adjust(wspace=0.32)
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    ax.computed_zorder = False
    ax.scatter(*centers.T, c=CMAP(norm(vals)), marker="s", s=5, depthshade=False, linewidths=0, zorder=1)
    ax.scatter(*xyz.T, c="red", marker="*", s=180, edgecolors="black", depthshade=False, zorder=10,
               label="probe points")
    lo, hi = centers.min(0), centers.max(0)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(np.maximum(hi - lo, 1e-6)))
    ax.view_init(elev=28, azim=-55)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    ax.set_title("Isometric (part surface)")
    ax.legend(loc="upper left")

    # Top-down: infill of the highest occupied voxel in each column.
    nx, ny, nz = occ.shape
    any_col = occ.any(axis=2)
    top_k = nz - 1 - np.argmax(occ[:, :, ::-1], axis=2)
    top = np.where(any_col, infill[labels[np.arange(nx)[:, None], np.arange(ny)[None, :], top_k]], np.nan)
    ext_xy = [origin[0] - pitch / 2, origin[0] + (nx - 0.5) * pitch, origin[1] - pitch / 2, origin[1] + (ny - 0.5) * pitch]
    a2 = fig.add_subplot(1, 3, 2)
    a2.imshow(top.T, origin="lower", extent=ext_xy, cmap=CMAP, norm=norm, interpolation="nearest")
    a2.scatter(xyz[:, 0], xyz[:, 1], c="red", marker="*", s=160, edgecolors="black")
    a2.set_title("Top-down (topmost voxel per column)")
    a2.set_xlabel("x (mm)")
    a2.set_ylabel("y (mm)")

    # Vertical x-z slice through the part centre.
    j = int(round(np.nonzero(occ)[1].mean()))
    sl = np.where(occ[:, j, :], infill[labels[:, j, :]], np.nan)
    y_c = origin[1] + j * pitch
    ext_xz = [ext_xy[0], ext_xy[1], origin[2] - pitch / 2, origin[2] + (nz - 0.5) * pitch]
    a3 = fig.add_subplot(1, 3, 3)
    a3.imshow(sl.T, origin="lower", extent=ext_xz, cmap=CMAP, norm=norm, interpolation="nearest")
    near = np.abs(xyz[:, 1] - y_c) < 3 * pitch
    a3.scatter(xyz[near, 0], xyz[near, 2], c="red", marker="*", s=160, edgecolors="black")
    a3.set_title(f"Side slice at y = {y_c:.1f} mm")
    a3.set_xlabel("x (mm)")
    a3.set_ylabel("z (mm)")
    a3.set_aspect("equal")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    fig.colorbar(sm, ax=fig.axes, fraction=0.02, pad=0.03, label="Assigned infill density (%)")
    kind = "Column (top-down probe)" if mode == "columns" else "Volumetric Voronoi"
    fig.suptitle(f"{kind} regions, coloured by assigned infill")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def export_colored_regions(display_meshes, infill, path):
    """Part split into coloured region solids; opens in Windows 3D Viewer / any glTF viewer."""
    norm = Normalize(0, 100)
    scene = trimesh.Scene()
    for i, (m, f) in enumerate(zip(display_meshes, infill)):
        if m is None:
            continue
        m = m.copy()
        m.visual.face_colors = (np.array(CMAP(norm(f))) * 255).astype(np.uint8)
        scene.add_geometry(m, node_name=f"region_{i:02d}_{f:.0f}pct")
    scene.export(path)


# ============================================================================ STAGE C
# Target: PrusaSlicer 2.5.0. Its 3MF stores every volume of an object in one merged mesh, and
# Metadata/Slic3r_PE_model.config (PrusaSlicer's counterpart of Bambu's model_settings.config) maps
# triangle-index ranges to volumes. Any volume metadata key that is not reserved is applied as a print
# setting for that volume, which is how each modifier's fill_density override reaches the slicer.
PRUSA_PINNED = "2.5.0"
DEFAULT_PROFILE = ROOT / "profiles" / "prusaslicer_2.5.0_tpu_generic_flex.ini"
_IDENTITY_4x4 = "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"


def find_prusaslicer(explicit=None):
    candidates = [explicit, os.environ.get("PRUSASLICER_CONSOLE"),
                  Path(os.environ.get("LOCALAPPDATA", "")) / "Programs/Prusa3D/PrusaSlicer/prusa-slicer-console.exe",
                  Path(os.environ.get("PROGRAMFILES", "C:/Program Files")) / "Prusa3D/PrusaSlicer/prusa-slicer-console.exe",
                  shutil.which("prusa-slicer-console")]
    for c in candidates:
        if c and Path(c).is_file():
            return str(c)
    raise FileNotFoundError("prusa-slicer-console not found; pass --slicer or set PRUSASLICER_CONSOLE")


def check_slicer_version(exe, allow_mismatch=False):
    proc = subprocess.run([exe, "--help"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    m = re.search(r"PrusaSlicer-(\d+\.\d+\.\d+)", proc.stdout + proc.stderr)
    version = m.group(1) if m else "unknown"
    if version != PRUSA_PINNED and not allow_mismatch:
        raise RuntimeError(f"PrusaSlicer {version} found; pinned to {PRUSA_PINNED} (use --allow-version-mismatch)")
    return version


def read_profile(path):
    cfg = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg


def bed_center(cfg):
    pts = np.array([[float(c) for c in p.split("x")] for p in cfg["bed_shape"].split(",")])
    return (pts.min(0) + pts.max(0)) / 2


def _fmt(v):
    return f"{v:.5f}".rstrip("0").rstrip(".") or "0"


def write_modifier_3mf(path, base_mesh, region_meshes_, infill, placement, as_parts=False, name="scanned_part",
                       extra_settings=None, names=None):
    """One volume per region, each with its own fill_density (plus any extra_settings).

    as_parts=False: base part + one ParameterModifier per region.
    as_parts=True:  no base; each region mesh is itself a ModelPart (the part already split into pieces).
    """
    extra = extra_settings or {}
    names = names or [f"region_{i:02d}" for i in range(len(region_meshes_))]
    if as_parts:
        volumes = [(nm, m, "ModelPart", {"fill_density": f"{f:g}%", **extra})
                   for nm, m, f in zip(names, region_meshes_, infill) if m is not None]
    else:
        volumes = [("base", base_mesh, "ModelPart", {})]
        volumes += [(nm, m, "ParameterModifier", {"fill_density": f"{f:g}%", **extra})
                    for nm, m, f in zip(names, region_meshes_, infill)]
    verts, tris, cfg = [], [], []
    v_off = t_off = 0
    for vname, mesh, vtype, settings in volumes:
        v = np.asarray(mesh.vertices, dtype=float)
        f = np.asarray(mesh.faces, dtype=np.int64) + v_off
        verts.extend(f'     <vertex x="{_fmt(a)}" y="{_fmt(b)}" z="{_fmt(c)}"/>' for a, b, c in v)
        tris.extend(f'     <triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in f)
        meta = [("name", vname), ("volume_type", vtype), ("matrix", _IDENTITY_4x4)]
        if vtype == "ParameterModifier":
            meta.insert(1, ("modifier", "1"))
        meta += list(settings.items())
        cfg.append(f'  <volume firstid="{t_off}" lastid="{t_off + len(f) - 1}">')
        cfg += [f'   <metadata type="volume" key={quoteattr(k)} value={quoteattr(str(val))}/>' for k, val in meta]
        cfg += ['   <mesh edges_fixed="0" degenerate_facets="0" facets_removed="0" facets_reversed="0" backwards_edges="0"/>',
                "  </volume>"]
        v_off += len(v)
        t_off += len(f)
    tx, ty, tz = placement
    model = "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" '
        'xmlns:slic3rpe="http://schemas.slic3r.org/3mf/2017/06">',
        ' <metadata name="slic3rpe:Version3mf">1</metadata>',
        f' <metadata name="Title">{escape(name)}</metadata>',
        ' <resources>', '  <object id="1" type="model">', '   <mesh>',
        '    <vertices>', *verts, '    </vertices>', '    <triangles>', *tris, '    </triangles>',
        '   </mesh>', '  </object>', ' </resources>',
        ' <build>', f'  <item objectid="1" transform="1 0 0 0 1 0 0 0 1 {_fmt(tx)} {_fmt(ty)} {_fmt(tz)}" printable="1"/>',
        ' </build>', '</model>', ''])
    config = "\n".join(['<?xml version="1.0" encoding="UTF-8"?>', '<config>', ' <object id="1" instances_count="1">',
                        f'  <metadata type="object" key="name" value={quoteattr(name)}/>', *cfg,
                        ' </object>', '</config>', ''])
    content_types = ('<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
                     ' <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n'
                     ' <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>\n</Types>\n')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            ' <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n'
            '</Relationships>\n')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", model)
        z.writestr("Metadata/Slic3r_PE_model.config", config)


def _run_slicer(exe, args, log_path):
    proc = subprocess.run([exe, "--loglevel", "2", *map(str, args)], capture_output=True, text=True,
                          encoding="utf-8", errors="replace")
    Path(log_path).write_text(f"$ {exe} {' '.join(map(str, args))}\n\n{proc.stdout}\n{proc.stderr}", encoding="utf-8")
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-12:]
        raise RuntimeError(f"PrusaSlicer exited {proc.returncode} (see {log_path}):\n" + "\n".join(tail))


MODIFIER_LIMIT = 12


def split_into_parts(mesh, region_solids, grow=0.05):
    """The part cut into one watertight piece per region (part mesh intersected with each region solid).

    Used instead of modifiers above MODIFIER_LIMIT regions: PrusaSlicer 2.5 spawns a combined region
    for every modifier whose bounding box overlaps an existing region, so slicing time explodes with
    Voronoi cells (10 regions: 1.6 s, 30: >4 min, 40: ~13 h then a crash), whereas separate parts stay
    linear. Solids are grown by `grow` mm so neighbouring pieces overlap slightly; PrusaSlicer trims
    overlapping parts and still walls only the outer surface (checked: no internal perimeters).
    """
    parts = []
    for solid in region_solids:
        s = solid.copy()
        s.vertices += s.vertex_normals * grow
        piece = trimesh.boolean.intersection([mesh, s], engine="manifold")
        parts.append(piece if piece is not None and len(piece.faces) else None)
    total = sum(p.volume for p in parts if p is not None)
    if abs(total - mesh.volume) > 0.02 * mesh.volume:
        raise RuntimeError(f"Split parts cover {total:.0f} mm^3 of a {mesh.volume:.0f} mm^3 part")
    return parts


def roundtrip_check(exe, model_3mf, run_dir, expected_infill, vtype="ParameterModifier"):
    """PrusaSlicer reloads and re-exports the 3MF; every region volume and its density must survive."""
    out = run_dir / "stageC_roundtrip_by_prusaslicer.3mf"
    _run_slicer(exe, ["--export-3mf", "--dont-arrange", "--output", out, model_3mf], run_dir / "stageC_roundtrip.log")
    with zipfile.ZipFile(out) as z:
        root = ET.fromstring(z.read("Metadata/Slic3r_PE_model.config"))
    vols = [{m.get("key"): m.get("value") for m in v.findall("metadata")} for v in root.iter("volume")]
    got = sorted(float(v["fill_density"].rstrip("%")) for v in vols
                 if v.get("volume_type") == vtype and "fill_density" in v)
    if len(got) != len(expected_infill) or not np.allclose(got, sorted(expected_infill)):
        raise RuntimeError(f"PrusaSlicer did not preserve the region volumes: {got}")
    return len(got)


def slice_with_cli(exe, model_3mf, profile, out_gcode, log_path):
    _run_slicer(exe, ["--export-gcode", "--load", profile, "--dont-arrange", "--output", out_gcode, model_3mf], log_path)


# ============================================================================ STAGE D
_WORD = re.compile(r"([A-Z])(-?\d*\.?\d+)")
SPARSE_TYPES = {"Internal infill", "Sparse infill"}  # PrusaSlicer, Bambu/Orca


def read_gcode_text(path):
    """Plain .gcode, or the G-code embedded in a Bambu/Orca .gcode.3mf."""
    path = Path(path)
    if path.name.lower().endswith(".gcode.3mf"):
        with zipfile.ZipFile(path) as z:
            names = sorted(n for n in z.namelist() if n.lower().endswith(".gcode"))
            if not names:
                raise ValueError(f"No G-code inside {path}")
            return z.read(names[0]).decode("utf-8", "replace")
    return path.read_text(encoding="utf-8", errors="replace")


def slicer_setting(text, keys, default=None):
    """First matching '; key = value' from the slicer's config dump at the end of the G-code."""
    for key in keys:
        m = re.search(rf"^;\s*{re.escape(key)}\s*=\s*([^\n,]+)", text, re.M)
        if m:
            return float(m.group(1).strip().rstrip("%"))
    return default


def parse_extrusions(text):
    """{feature type: array[(x0, y0, x1, y1, z, e)]} for every extruding XY move. Relative or absolute E."""
    x = y = z = e_last = 0.0
    relative_e = False
    cur = None
    out = {}
    for raw in text.splitlines():
        s = raw.strip()
        if s.startswith(";TYPE:"):
            cur = s[6:].strip()
            continue
        if s.startswith("; FEATURE:"):
            cur = s[10:].strip()
            continue
        code = s.split(";", 1)[0].strip().upper()
        if not code:
            continue
        w = dict(_WORD.findall(code))
        g, mcode = w.get("G"), w.get("M")
        if mcode == "83":
            relative_e = True
        elif mcode == "82":
            relative_e = False
        elif g == "92" and "E" in w:
            e_last = float(w["E"])
        elif g in ("0", "1"):
            nx, ny, nz = float(w.get("X", x)), float(w.get("Y", y)), float(w.get("Z", z))
            e = 0.0
            if "E" in w:
                ev = float(w["E"])
                if relative_e:
                    e = ev
                else:
                    e, e_last = ev - e_last, ev
            if e > 0 and (nx != x or ny != y) and cur:
                out.setdefault(cur, []).append((x, y, nx, ny, nz, e))
            x, y, z = nx, ny, nz
    return {t: np.array(v, dtype=float) for t, v in out.items()}


def bed_offset(segs_by_type, mesh):
    """Bed -> part-frame XY offset from the external perimeter's bounding-box centre."""
    ext = next((segs_by_type[t] for t in ("External perimeter", "Outer wall") if t in segs_by_type), None)
    if ext is None:
        raise ValueError("No external perimeter found in G-code")
    xs, ys = np.r_[ext[:, 0], ext[:, 2]], np.r_[ext[:, 1], ext[:, 3]]
    centre = np.array([(xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2])
    return centre - (mesh.bounds[0, :2] + mesh.bounds[1, :2]) / 2


def voxel_label_at(grid, labels, pts):
    """Step 6 label of the voxel containing each point (N x 3, part frame)."""
    idx = np.round((np.asarray(pts) - grid["origin"]) / grid["pitch"]).astype(int)
    idx = np.clip(idx, 0, np.array(labels.shape) - 1)
    return labels[idx[:, 0], idx[:, 1], idx[:, 2]]


def mesh_section(mesh, z):
    sec = mesh.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
    if sec is None:
        return Polygon()
    planar, _ = sec.to_2D(to_2D=trimesh.transformations.translation_matrix([0, 0, -z]))
    return shapely.unary_union(planar.polygons_full)


def _sample_segments(segs, step):
    """Split each segment into pieces <= step mm; return piece midpoints and their share of E."""
    lengths = np.hypot(segs[:, 2] - segs[:, 0], segs[:, 3] - segs[:, 1])
    n = np.maximum(1, np.ceil(lengths / step).astype(int))
    rep = np.repeat(np.arange(len(segs)), n)
    k = np.concatenate([np.arange(c) for c in n])
    t = (k + 0.5) / n[rep]
    s = segs[rep]
    pts = np.column_stack([s[:, 0] + t * (s[:, 2] - s[:, 0]), s[:, 1] + t * (s[:, 3] - s[:, 1])])
    return pts, s[:, 5] / n[rep], s[:, 4]


def measure_region_infill(text, mesh, region_of, n_regions, sample_step=0.25, area_step=0.3, opening=0.15):
    """Realised sparse-infill density (%) per region from sliced G-code.

    Both terms are binned with region_of (Step 6 voxel labels in 3D mode, the exact column cells in
    column mode), i.e. the same region boundaries the modifier meshes were built from. Numerator: deposited sparse-infill area, from E corrected for the extrusion
    multiplier. Denominator: per layer, the part cross-section minus the footprint of every non-sparse
    extrusion, area-sampled on a fine grid. The footprint buffer distance is measured from the G-code
    itself (how far sparse lines stop from neighbouring perimeters/solid infill).
    """
    layer_h = slicer_setting(text, ["layer_height"])
    fil_d = slicer_setting(text, ["filament_diameter"], 1.75)
    em = slicer_setting(text, ["extrusion_multiplier", "filament_flow_ratio"], 1.0)
    segs = parse_extrusions(text)
    sparse_t = next(t for t in SPARSE_TYPES if t in segs)
    offset = bed_offset(segs, mesh)

    def shift(a):
        a = a.copy()
        a[:, [0, 2]] -= offset[0]
        a[:, [1, 3]] -= offset[1]
        return a

    sparse = shift(segs[sparse_t])
    other = np.vstack([shift(v) for t, v in segs.items() if t != sparse_t])
    area_per_e = math.pi * (fil_d / 2) ** 2 / (layer_h * em)

    # Numerator.
    pts, e_share, z_top = _sample_segments(sparse, sample_step)
    lab = region_of(np.column_stack([pts, z_top - layer_h / 2]))
    deposited = np.bincount(lab, weights=e_share * area_per_e, minlength=n_regions)

    # Buffer distance: median gap between sparse-line endpoints and the nearest non-sparse centreline.
    layers = np.unique(np.round(sparse[:, 4], 4))
    gaps = []
    for z in layers[:: max(1, len(layers) // 15)]:
        o = other[np.isclose(other[:, 4], z)]
        s = sparse[np.isclose(sparse[:, 4], z)]
        if len(o) == 0 or len(s) == 0:
            continue
        op, _, _ = _sample_segments(o, 0.05)
        d, _ = cKDTree(op).query(np.r_[s[:, :2], s[:, 2:4]])
        gaps.extend(d[d < 1.5])
    buffer_d = float(np.median(gaps))

    fillable = np.zeros(n_regions)
    for z in layers:
        zc = z - layer_h / 2
        o = other[np.isclose(other[:, 4], z)]
        zone = mesh_section(mesh, zc)
        if len(o):
            zone = zone.difference(MultiLineString([[(a, b), (c, d)] for a, b, c, d in o[:, :4]]).buffer(buffer_d, quad_segs=4))
        zone = zone.buffer(-opening).buffer(opening)  # drop slivers outside the outer perimeter
        if zone.is_empty:
            continue
        x0, y0, x1, y1 = zone.bounds
        gx = np.arange(np.floor(x0 / area_step), np.ceil(x1 / area_step)) * area_step + area_step / 2
        gy = np.arange(np.floor(y0 / area_step), np.ceil(y1 / area_step)) * area_step + area_step / 2
        X, Y = (a.ravel() for a in np.meshgrid(gx, gy))
        shapely.prepare(zone)
        inside = shapely.contains_xy(zone, X, Y)
        if inside.any():
            lab = region_of(np.column_stack([X[inside], Y[inside], np.full(inside.sum(), zc)]))
            fillable += np.bincount(lab, minlength=n_regions) * area_step ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        realised = np.where(fillable > 0, 100.0 * deposited / fillable, np.nan)
    info = {"sparse_type": sparse_t, "layer_height": layer_h, "filament_diameter": fil_d,
            "extrusion_multiplier": em, "bed_offset_xy": offset.round(3).tolist(),
            "n_sparse_layers": int(len(layers)), "sparse_footprint_buffer_mm": round(buffer_d, 3),
            "fillable_area_mm2": fillable.round(1).tolist()}
    return realised, info, sparse


def plot_verification(target, realised, tolerance, path):
    idx = np.arange(len(target))
    err = realised - target
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=180, gridspec_kw={"width_ratios": [1.6, 1]})
    a1.bar(idx - 0.2, target, 0.4, label="Target (from stiffness)", color="#9aa5b1")
    a1.bar(idx + 0.2, realised, 0.4, label="Realised (from G-code)", color="#1f6f8b")
    for i, e in zip(idx, err):
        a1.text(i, max(target[i], realised[i]) + 1.5, f"{e:+.1f}", ha="center", fontsize=8, fontweight="bold",
                color="black" if abs(e) <= tolerance else "#c0392b")
    a1.set_xticks(idx, [f"R{i}" for i in idx])
    a1.set_ylabel("Infill density (%)")
    a1.set_ylim(0, np.nanmax([target.max(), np.nanmax(realised)]) + 20)
    a1.set_title("Per-region infill: target vs. realised (labels = error, pp)")
    a1.legend(loc="upper center", ncol=2)
    hi = np.nanmax([target.max(), np.nanmax(realised)]) + 5
    a2.fill_between([0, hi], [-tolerance, hi - tolerance], [tolerance, hi + tolerance], color="#2ecc71", alpha=0.15,
                    label=f"+/-{tolerance:g} pp")
    a2.plot([0, hi], [0, hi], "k--", linewidth=1, label="y = x")
    a2.scatter(target, realised, color="#1f6f8b", zorder=3)
    for i in idx:
        a2.annotate(f"R{i}", (target[i], realised[i]), textcoords="offset points", xytext=(4, 4), fontsize=8)
    a2.set_xlim(0, hi)
    a2.set_ylim(0, hi)
    a2.set_aspect("equal")
    a2.set_xlabel("Target infill (%)")
    a2.set_ylabel("Realised infill (%)")
    a2.set_title("Agreement")
    a2.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _segment_regions(sparse, region_of, layer_h):
    mids = np.column_stack([(sparse[:, 0] + sparse[:, 2]) / 2, (sparse[:, 1] + sparse[:, 3]) / 2,
                            sparse[:, 4] - layer_h / 2])
    return region_of(mids)


def plot_toolpath_layers(sparse, region_of, grid, labels, infill, layer_h, xyz, path, cells2d=None):
    """Two consecutive mid-height layers: sliced infill lines coloured by their region, region boundaries in black."""
    from matplotlib.collections import LineCollection
    infill = np.asarray(infill)
    norm = Normalize(0, 100)
    zs = np.unique(np.round(sparse[:, 4], 4))
    picks = zs[len(zs) // 2: len(zs) // 2 + 2]
    pitch, origin = grid["pitch"], grid["origin"]
    nx, ny, _ = labels.shape
    extent = [origin[0] - pitch / 2, origin[0] + (nx - 0.5) * pitch, origin[1] - pitch / 2, origin[1] + (ny - 0.5) * pitch]
    fig, axes = plt.subplots(1, len(picks), figsize=(8.5 * len(picks), 8.5), dpi=170, squeeze=False)
    for ax, z in zip(axes[0], picks):
        if cells2d is not None:
            walls = [np.asarray(c.exterior.coords) for c in cells2d]
        else:
            # Region boundaries: the voxel faces between differently labelled voxels inside the part.
            k = int(np.clip(round((z - layer_h / 2 - origin[2]) / pitch), 0, labels.shape[2] - 1))
            lab2d, occ2d = labels[:, :, k], grid["occupied"][:, :, k]
            bx = np.argwhere((lab2d[1:, :] != lab2d[:-1, :]) & (occ2d[1:, :] | occ2d[:-1, :]))
            by = np.argwhere((lab2d[:, 1:] != lab2d[:, :-1]) & (occ2d[:, 1:] | occ2d[:, :-1]))
            cx = lambda i: origin[0] + i * pitch
            cy = lambda j: origin[1] + j * pitch
            walls = [[(cx(i + 0.5), cy(j - 0.5)), (cx(i + 0.5), cy(j + 0.5))] for i, j in bx]
            walls += [[(cx(i - 0.5), cy(j + 0.5)), (cx(i + 0.5), cy(j + 0.5))] for i, j in by]
        ax.add_collection(LineCollection(walls, colors="black", linewidths=1.3, zorder=3))
        seg = sparse[np.isclose(sparse[:, 4], z)]
        reg = _segment_regions(seg, region_of, layer_h)
        ax.add_collection(LineCollection(seg[:, :4].reshape(-1, 2, 2), colors=CMAP(norm(infill[reg])), linewidths=0.7, zorder=2))
        top = xyz[:, :2]
        ax.scatter(top[:, 0], top[:, 1], c="red", marker="*", s=90, edgecolors="black", zorder=4)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"Sliced infill, layer z = {z:.2f} mm")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    fig.colorbar(sm, ax=axes[0].tolist(), fraction=0.025, pad=0.02, label="Region's assigned infill (%)")
    fig.suptitle("PrusaSlicer toolpaths coloured by region; line spacing changes at each region boundary (black)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def toolpaths_3d_html(sparse, region_of, infill, layer_h, xyz, path, max_segments=60000):
    """Rotatable 3D view (opens in a browser) of the sliced infill, coloured by region."""
    import plotly.graph_objects as go
    zs = np.unique(np.round(sparse[:, 4], 4))
    stride = max(1, int(np.ceil(len(sparse) / max_segments)))
    keep = np.isin(np.round(sparse[:, 4], 4), zs[::stride])
    s = sparse[keep]
    vals = np.asarray(infill)[_segment_regions(s, region_of, layer_h)]
    nan = np.full(len(s), np.nan)
    x = np.column_stack([s[:, 0], s[:, 2], nan]).ravel()
    y = np.column_stack([s[:, 1], s[:, 3], nan]).ravel()
    z = np.column_stack([s[:, 4], s[:, 4], nan]).ravel()
    c = np.repeat(vals, 3)
    fig = go.Figure(go.Scatter3d(x=x, y=y, z=z, mode="lines", hoverinfo="skip",
                                 line=dict(color=c, colorscale="Viridis", cmin=0, cmax=100, width=2,
                                           colorbar=dict(title="Infill %"))))
    fig.add_trace(go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="markers+text",
                               marker=dict(size=5, color="red"), text=[f"R{i}" for i in range(len(xyz))],
                               name="probe points"))
    fig.update_layout(title=f"Sliced infill coloured by region (every {stride} layer{'s' if stride > 1 else ''})",
                      scene=dict(aspectmode="data", xaxis_title="x (mm)", yaxis_title="y (mm)", zaxis_title="z (mm)"),
                      showlegend=False, margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(path, include_plotlyjs=True)


def write_verification(records, target, realised, info, tolerance, run_dir):
    rows = []
    for r, t, a, fa in zip(records, target, realised, info["fillable_area_mm2"]):
        err = a - t
        rows.append({"region": f"R{r['id']}", "x": round(r["x"], 3), "y": round(r["y"], 3), "z": round(r["z"], 3),
                     "stiffness": round(r["stiffness"], 4), "target_infill_pct": round(float(t), 2),
                     "realised_infill_pct": round(float(a), 2), "error_pp": round(float(err), 2),
                     "fillable_area_mm2": fa, "within_tolerance": bool(abs(err) <= tolerance)})
    with open(run_dir / "stageD_verification.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    plot_verification(np.asarray(target), np.asarray(realised), tolerance, run_dir / "stageD_verification.png")
    print(f"\n  {'region':>6} {'k':>6} {'target%':>8} {'realised%':>10} {'err pp':>7}")
    for r in rows:
        flag = "" if r["within_tolerance"] else "  <-- outside tolerance"
        print(f"  {r['region']:>6} {r['stiffness']:>6.2f} {r['target_infill_pct']:>8.1f} "
              f"{r['realised_infill_pct']:>10.2f} {r['error_pp']:>+7.2f}{flag}")
    errs = np.array([r["error_pp"] for r in rows])
    print(f"  MAE {np.nanmean(np.abs(errs)):.2f} pp, max {np.nanmax(np.abs(errs)):.2f} pp, "
          f"{sum(r['within_tolerance'] for r in rows)}/{len(rows)} within +/-{tolerance:g} pp")
    return rows


# ============================================================================ DRIVER
def run_stage_a(args, rng, run_dir):
    print("[Stage A] geometry + stiffness points")
    if args.dummy_cube:
        mesh = trimesh.creation.box(extents=DUMMY_CUBE_SIZE)
        mesh.apply_translation(np.array(DUMMY_CUBE_SIZE) / 2)
        records = make_point_records(np.array([p[:3] for p in DUMMY_POINTS]), [p[3] for p in DUMMY_POINTS])
        print(f"  dummy cube {DUMMY_CUBE_SIZE} mm with {len(records)} hardcoded points")
    else:
        mesh = repair_mesh(load_scan(args.scan, args.units, args.scale))
        mesh.apply_transform(registration_transform())
        if args.points:
            records = load_points_csv(args.points)
            print(f"  loaded {len(records)} stiffness points from {args.points}")
            mesh, records, shift = to_part_frame(mesh, records)
        else:
            mesh, _, shift = to_part_frame(mesh, [])
            records = sample_stiffness_points(mesh, args.n_points, rng, min_spacing=args.min_spacing,
                                              top_only=not args.all_surfaces)
            print(f"  sampled {len(records)} random stiffness points on {'all' if args.all_surfaces else 'upward-facing'} surface")
    mesh.export(run_dir / "stageA_mesh.stl")
    save_points_csv(records, run_dir / "stageA_points.csv")
    return mesh, records


def run_stage_b(args, mesh, records, run_dir, to_infill=None):
    """to_infill: stiffness array (N/mm) -> infill array (%); default is the placeholder linear map."""
    print("[Stage B] region partition")
    pitch = args.voxel_pitch or round(float(mesh.extents.max()) / 60, 2)
    grid = voxelize(mesh, pitch)
    xyz = points_xyz(records)
    top = top_contact_mask(grid, xyz)
    part_xyz, mode = partition_points(xyz, args.partition, top)
    print(f"  {int(top.sum())}/{len(xyz)} points are top contacts -> partition mode: {mode}")
    labels = label_voxels(grid, part_xyz)
    target = np.round((to_infill or stiffness_to_infill)(np.array([r["stiffness"] for r in records])), 1)
    cells2d = None
    if mode == "columns":
        cells2d, modifiers = column_modifiers(mesh, part_xyz)
        display = [mask_to_mesh((labels == i) & grid["occupied"], grid) if ((labels == i) & grid["occupied"]).any()
                   else None for i in range(len(records))]
    else:
        modifiers, display = region_meshes(grid, labels, len(records))
    regions = {"mode": mode, "grid": grid, "labels": labels, "cells2d": cells2d,
               "region_of": make_region_lookup(mode, part_xyz, grid, labels)}
    occ_counts = np.bincount(labels[grid["occupied"]], minlength=len(records))
    regions_dir = run_dir / "stageB_region_meshes"
    regions_dir.mkdir()
    for i, m in enumerate(modifiers):
        m.export(regions_dir / f"modifier_{i:02d}.stl")
    np.savez_compressed(run_dir / "stageB_voxel_labels.npz", labels=labels, occupied=grid["occupied"],
                        origin=grid["origin"], pitch=grid["pitch"])
    plot_regions_3d(grid, labels, target, xyz, run_dir / "stageB_regions_3d.png", mode=mode)
    export_colored_regions(display, target, run_dir / "stageB_regions_colored.glb")
    for i, r in enumerate(records):
        print(f"    R{i}: k={r['stiffness']:.2f} N/mm -> {target[i]:5.1f}% infill, "
              f"{occ_counts[i] * pitch ** 3:8.1f} mm^3 of part, modifier {len(modifiers[i].faces)} faces")
    if (occ_counts == 0).any():
        print(f"    note: regions {np.flatnonzero(occ_counts == 0).tolist()} own no voxels inside the part")
    return regions, target, modifiers


def run_stage_c(args, mesh, modifiers, target, run_dir):
    print("[Stage C] modifier 3MF + slicing")
    cfg = read_profile(args.profile)
    shutil.copy(args.profile, run_dir / Path(args.profile).name)
    placement = (*(bed_center(cfg) - (mesh.bounds[0, :2] + mesh.bounds[1, :2]) / 2), 0.0)
    model_3mf = run_dir / "stageC_graded_infill.3mf"
    strategy = args.region_volumes
    if strategy == "auto":
        strategy = "modifiers" if len(modifiers) <= MODIFIER_LIMIT else "parts"
    if strategy == "parts":
        parts = split_into_parts(mesh, modifiers)
        expected = [f for p, f in zip(parts, target) if p is not None]
        write_modifier_3mf(model_3mf, mesh, parts, target, placement, as_parts=True)
        vtype = "ModelPart"
        print(f"  {len(modifiers)} regions > {MODIFIER_LIMIT}: part split into {len(expected)} pieces "
              f"(PrusaSlicer modifiers scale combinatorially)")
    else:
        expected = list(target)
        write_modifier_3mf(model_3mf, mesh, modifiers, target, placement)
        vtype = "ParameterModifier"
    exe = find_prusaslicer(args.slicer)
    version = check_slicer_version(exe, args.allow_version_mismatch)
    n = roundtrip_check(exe, model_3mf, run_dir, expected, vtype)
    print(f"  PrusaSlicer {version} reloaded the 3MF and kept all {n} region volumes ({strategy}) with their fill_density")
    gcode = run_dir / "stageC_graded_infill.gcode"
    slice_with_cli(exe, model_3mf, args.profile, gcode, run_dir / "stageC_slice.log")
    print(f"  sliced -> {gcode.name} ({gcode.stat().st_size / 1e6:.1f} MB)")
    return gcode, version, placement


def run_stage_d(args, mesh, records, regions, target, gcode, run_dir):
    print("[Stage D] verification")
    text = read_gcode_text(gcode)
    region_of = regions["region_of"]
    realised, info, sparse = measure_region_infill(text, mesh, region_of, len(records))
    print(f"  {info['n_sparse_layers']} sparse-infill layers, extrusion multiplier {info['extrusion_multiplier']:g} "
          f"(divided out), bed offset {info['bed_offset_xy']}")
    rows = write_verification(records, target, realised, info, args.tolerance, run_dir)
    xyz, lh = points_xyz(records), info["layer_height"]
    plot_toolpath_layers(sparse, region_of, regions["grid"], regions["labels"], target, lh, xyz,
                         run_dir / "stageD_toolpaths_by_region.png", cells2d=regions["cells2d"])
    toolpaths_3d_html(sparse, region_of, target, lh, xyz, run_dir / "stageD_toolpaths_3d.html")
    return realised, info, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--dummy-cube", action="store_true")
    src.add_argument("--scan", type=Path, help=".glb/.gltf/.stl/.obj scan")
    ap.add_argument("--points", type=Path, help="CSV with x,y,z,stiffness (registered frame); random if omitted")
    ap.add_argument("--n-points", type=int, default=10)
    ap.add_argument("--min-spacing", type=float, default=8.0)
    ap.add_argument("--all-surfaces", action="store_true", help="sample points on any surface, not only upward-facing")
    ap.add_argument("--units", choices=["auto", "m", "mm"], default="auto")
    ap.add_argument("--scale", type=float, default=1.0, help="print scale applied after unit conversion, e.g. 0.3")
    ap.add_argument("--voxel-pitch", type=float, default=None, help="mm; default = largest extent / 60")
    ap.add_argument("--partition", choices=["auto", "columns", "3d"], default="auto")
    ap.add_argument("--region-volumes", choices=["auto", "modifiers", "parts"], default="auto",
                    help=f"auto = modifiers up to {MODIFIER_LIMIT} regions, split parts above")
    ap.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    ap.add_argument("--slicer", default=None, help="path to prusa-slicer-console")
    ap.add_argument("--allow-version-mismatch", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--tolerance", type=float, default=5.0)
    ap.add_argument("--stop-after", choices=["A", "B", "C", "D"], default="D")
    ap.add_argument("--out", type=Path, default=ROOT / "out_scan")
    args = ap.parse_args()

    tag = "cube" if args.dummy_cube else args.scan.stem
    run_dir = args.out / f"run_{datetime.now():%Y%m%d_%H%M%S}_{tag}_seed{args.seed}"
    run_dir.mkdir(parents=True)
    rng = np.random.default_rng(args.seed)
    meta = {"args": {k: str(v) for k, v in vars(args).items()}, "k_range": K_RANGE, "infill_range": INFILL_RANGE}

    stages = "ABCD"[: "ABCD".index(args.stop_after) + 1]
    mesh, records = run_stage_a(args, rng, run_dir)
    if "B" in stages:
        regions, target, modifiers = run_stage_b(args, mesh, records, run_dir)
        meta.update(voxel_pitch_mm=regions["grid"]["pitch"], partition=regions["mode"], target_infill_pct=target.tolist())
    if "C" in stages:
        gcode, version, placement = run_stage_c(args, mesh, modifiers, target, run_dir)
        meta.update(slicer=f"PrusaSlicer {version}", profile=Path(args.profile).name, placement_mm=list(placement))
    if "D" in stages:
        realised, info, rows = run_stage_d(args, mesh, records, regions, target, gcode, run_dir)
        err = np.array([r["error_pp"] for r in rows])
        meta.update(verification={k: v for k, v in info.items() if k != "fillable_area_mm2"},
                    mean_abs_error_pp=round(float(np.nanmean(np.abs(err))), 2),
                    max_abs_error_pp=round(float(np.nanmax(np.abs(err))), 2),
                    regions_within_tolerance=int(sum(r["within_tolerance"] for r in rows)))
    (run_dir / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Artifacts: {run_dir}")


if __name__ == "__main__":
    main()
