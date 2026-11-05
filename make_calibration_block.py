"""Calibration block: one block split into patches (Voronoi cells or a grid), each at a fixed infill density.

Probe each printed patch at its marked probe point to get the real stiffness-vs-infill curve (and the
softest / stiffest stiffness the material can reach), which replaces the placeholder mapping in
scan_to_infill.py.

  python make_calibration_block.py                      # 8 Voronoi cells
  python make_calibration_block.py --layout grid        # 3 x 2 rectangles
"""
import argparse
import json
import shutil
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import quoteattr

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapely
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree
from shapely.geometry import MultiPoint, Point, Polygon, box
from shapely.ops import polylabel

import scan_to_infill as s

ROOT = Path(__file__).resolve().parent


def block_outline(size, notch):
    """Block footprint with the front-left corner chamfered as an orientation mark."""
    sx, sy, _ = size
    return box(0, 0, sx, sy).difference(Polygon([(0, 0), (notch, 0), (0, notch)]))


def grid_layout(outline, size, grid):
    sx, sy, _ = size
    cols, rows = grid
    w, h = sx / cols, sy / rows
    cells = [box(c * w, r * h, (c + 1) * w, (r + 1) * h).intersection(outline)
             for r in range(rows) for c in range(cols)]
    seeds = np.array([[c.centroid.x, c.centroid.y] for c in cells])
    return seeds, cells


def voronoi_layout_varied(outline, n, min_clear_diameter, rng, candidates=150, **kw):
    """Of `candidates` valid layouts, the one whose cells differ most in area (largest / smallest).

    A wide spread in cell size makes the area-ordered infill (largest cell = sparsest) easy to see.
    """
    best, best_ratio = None, 0.0
    for _ in range(candidates):
        seeds, cells = voronoi_layout(outline, n, min_clear_diameter, rng, **kw)
        areas = [c.area for c in cells]
        ratio = max(areas) / min(areas)
        if ratio > best_ratio:
            best, best_ratio = (seeds, cells), ratio
    print(f"  picked layout with area ratio {best_ratio:.2f} (largest / smallest cell) of {candidates} candidates")
    return best


def voronoi_layout(outline, n, min_clear_diameter, rng, edge_margin=6.0, max_layouts=5000):
    """Random seeds whose clipped Voronoi cells each fit a circle of min_clear_diameter.

    That clear zone is what makes a patch measurable: at 10 % infill the line pitch is ~4 mm, so a
    probe needs several pitches of uniform infill around it, away from the neighbouring patches.
    """
    x0, y0, x1, y1 = outline.bounds
    min_spacing = 0.9 * min_clear_diameter
    for _ in range(max_layouts):
        pts = []
        for _ in range(4000):
            p = rng.uniform((x0 + edge_margin, y0 + edge_margin), (x1 - edge_margin, y1 - edge_margin))
            if outline.contains(Point(p)) and all(np.hypot(*(p - q)) >= min_spacing for q in pts):
                pts.append(p)
                if len(pts) == n:
                    break
        if len(pts) < n:
            continue
        seeds = np.array(pts)
        vor = shapely.voronoi_polygons(MultiPoint(seeds), extend_to=outline.envelope.buffer(10), ordered=True)
        cells = [c.intersection(outline).simplify(1e-9) for c in vor.geoms]
        if min(probe_target(c)[1] for c in cells) * 2 >= min_clear_diameter:
            return seeds, cells
    raise SystemExit(f"No layout of {n} cells with {min_clear_diameter} mm clear zones; lower --n-cells or --min-clear")


def probe_target(cell):
    """Point inside the cell farthest from its edges, and that distance (clear radius)."""
    p = polylabel(cell, tolerance=0.05)
    return p, cell.exterior.distance(p)


def assign_infill(n, lo, hi, rng):
    """n distinct densities spread over [lo, hi] (rounded to 5 %), shuffled so neighbours contrast."""
    vals = np.unique(np.round(np.linspace(lo, hi, n) / 5) * 5)
    if len(vals) < n:
        vals = np.round(np.linspace(lo, hi, n), 1)
    return rng.permutation(vals).tolist()


def assign_infill_by_area(cells, lo, hi):
    """Densities spread over [lo, hi] (rounded to 5 %): the largest cell gets the lowest, the smallest the highest."""
    n = len(cells)
    vals = np.unique(np.round(np.linspace(lo, hi, n) / 5) * 5)
    if len(vals) < n:
        vals = np.round(np.linspace(lo, hi, n), 1)
    order = np.argsort([-c.area for c in cells])          # largest area first
    infill = np.empty(n)
    infill[order] = vals
    return infill.tolist()


def plot_block_map(outline, cells, infill, targets, path, title):
    """Printable map: each patch coloured by infill, with its number, density and probe point."""
    norm = Normalize(0, 100)
    cmap = plt.get_cmap("viridis")
    x0, y0, x1, y1 = outline.bounds
    fig, ax = plt.subplots(figsize=(11, 11 * (y1 - y0) / (x1 - x0) + 1.2), dpi=200)
    for i, (c, f, (p, r)) in enumerate(zip(cells, infill, targets)):
        ax.fill(*c.exterior.xy, color=cmap(norm(f)), edgecolor="black", linewidth=1.8)
        ax.plot(p.x, p.y, marker="+", color="red", markersize=16, markeredgewidth=2.5)
        ax.add_patch(plt.Circle((p.x, p.y), r, fill=False, linestyle="--", color="white", linewidth=1))
        txt_col = "white" if norm(f) < 0.55 else "black"
        ax.text(p.x, p.y + 3.2, f"#{i}  {f:g}%", ha="center", va="bottom", fontsize=13, fontweight="bold", color=txt_col)
        ax.text(p.x, p.y - 3.2, f"probe ({p.x:.1f}, {p.y:.1f})", ha="center", va="top", fontsize=9, color=txt_col)
    ax.plot(*outline.exterior.xy, color="black", linewidth=2.5)
    ax.annotate("chamfered corner\n(origin 0,0)", xy=(1.5, 1.5), xytext=(6, -9), fontsize=9,
                arrowprops=dict(arrowstyle="->"), annotation_clip=False)
    ax.set_xlim(x0 - 2, x1 + 2)
    ax.set_ylim(y0 - 12, y1 + 2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm), from the chamfered corner")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.02, label="Infill density (%)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_bambu_3mf(path, pieces, names, part_settings, centre_on_bed, title="calibration_block", object_settings=None,
                    project_settings=None):
    """Bambu Studio project: one object whose parts each carry their own settings.

    Layout copied from a 3MF exported by Bambu Studio 01.07.07.89: part meshes live in
    3D/Objects/object_1.model as separate <object>s, the main model references them as components,
    and Metadata/model_settings.config holds per-part settings (sparse_infill_density etc.).
    project_settings (the slicer's own --export-settings JSON) fixes printer, process and filament; without
    it the printer and filament stay whatever the user picks in the GUI.
    """
    allv = np.vstack([p.vertices for p in pieces])
    centre = (allv.min(0) + allv.max(0)) / 2
    fmt = lambda v: f"{v:.5f}".rstrip("0").rstrip(".") or "0"
    uid = lambda: str(uuid.uuid4())
    ns = ('xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" '
          'xmlns:slic3rpe="http://schemas.slic3r.org/3mf/2017/06" '
          'xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" requiredextensions="p"')
    obj_xml, comps = [], []
    for i, p in enumerate(pieces, start=1):
        v = p.vertices - centre
        obj_xml += [f'  <object id="{i}" p:UUID="{uid()}" type="model">', '   <mesh>', '    <vertices>']
        obj_xml += [f'     <vertex x="{fmt(a)}" y="{fmt(b)}" z="{fmt(c)}"/>' for a, b, c in v]
        obj_xml += ['    </vertices>', '    <triangles>']
        obj_xml += [f'     <triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in p.faces]
        obj_xml += ['    </triangles>', '   </mesh>', '  </object>']
        comps.append(f'    <component p:path="/3D/Objects/object_1.model" objectid="{i}" p:UUID="{uid()}" '
                     f'transform="1 0 0 0 1 0 0 0 1 0 0 0"/>')
    top_id = len(pieces) + 1
    sub_model = "\n".join([f'<?xml version="1.0" encoding="UTF-8"?>', f'<model unit="millimeter" xml:lang="en-US" {ns}>',
                           ' <metadata name="BambuStudio:3mfVersion">1</metadata>', ' <resources>', *obj_xml,
                           ' </resources>', ' <build/>', '</model>', ''])
    cx, cy = centre_on_bed
    main_model = "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>', f'<model unit="millimeter" xml:lang="en-US" {ns}>',
        ' <metadata name="Application">BambuStudio-01.07.07.89</metadata>',
        ' <metadata name="BambuStudio:3mfVersion">1</metadata>',
        f' <metadata name="Title">{title}</metadata>',
        ' <resources>', f'  <object id="{top_id}" p:UUID="{uid()}" type="model">', '   <components>', *comps,
        '   </components>', '  </object>', ' </resources>', f' <build p:UUID="{uid()}">',
        f'  <item objectid="{top_id}" p:UUID="{uid()}" transform="1 0 0 0 1 0 0 0 1 {fmt(cx)} {fmt(cy)} {fmt(centre[2])}" printable="1"/>',
        ' </build>', '</model>', ''])
    parts = []
    for i, (nm, st) in enumerate(zip(names, part_settings), start=1):
        parts.append(f'    <part id="{i}" subtype="normal_part">')
        meta = {"name": nm, "matrix": "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1", **st}
        parts += [f'      <metadata key={quoteattr(k)} value={quoteattr(str(val))}/>' for k, val in meta.items()]
        parts.append('    </part>')
    settings = "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>', '<config>', f'  <object id="{top_id}">',
        f'    <metadata key="name" value={quoteattr(title)}/>', '    <metadata key="extruder" value="1"/>',
        *[f'    <metadata key={quoteattr(k)} value={quoteattr(str(v))}/>' for k, v in (object_settings or {}).items()],
        *parts,
        '  </object>', '  <plate>', '    <metadata key="plater_id" value="1"/>', '    <metadata key="plater_name" value=""/>',
        '    <metadata key="locked" value="false"/>', '    <model_instance>',
        f'      <metadata key="object_id" value="{top_id}"/>', '      <metadata key="instance_id" value="0"/>',
        '    </model_instance>', '  </plate>', '  <assemble>', '  </assemble>', '</config>', ''])
    content_types = ('<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
                     ' <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n'
                     ' <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>\n'
                     ' <Default Extension="png" ContentType="image/png"/>\n'
                     ' <Default Extension="gcode" ContentType="text/x.gcode"/>\n</Types>\n')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            ' <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n'
            '</Relationships>\n')
    model_rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
                  ' <Relationship Target="/3D/Objects/object_1.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n'
                  '</Relationships>\n')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("3D/3dmodel.model", main_model)
        z.writestr("3D/_rels/3dmodel.model.rels", model_rels)
        z.writestr("3D/Objects/object_1.model", sub_model)
        z.writestr("Metadata/model_settings.config", settings)
        if project_settings is not None:
            # Printer / process / filament for the whole plate. The Bambu CLI segfaults without it.
            z.writestr("Metadata/project_settings.config", project_settings)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--size", type=float, nargs=3, default=(100.0, 60.0, 20.0), metavar=("X", "Y", "Z"))
    ap.add_argument("--layout", choices=["voronoi", "grid"], default="voronoi")
    ap.add_argument("--n-cells", type=int, default=8, help="Voronoi layout: number of patches")
    ap.add_argument("--min-clear", type=float, default=16.0,
                    help="Voronoi layout: every patch must fit a circle this wide (mm)")
    ap.add_argument("--grid", type=int, nargs=2, default=(3, 2), metavar=("COLS", "ROWS"))
    ap.add_argument("--infill", type=float, nargs="+", default=None,
                    help="one value per patch; default spreads 10-85 %% across the patches")
    ap.add_argument("--order", choices=["area", "shuffle"], default="area",
                    help="area: largest patch = lowest infill (easy to see); shuffle: random order")
    ap.add_argument("--notch", type=float, default=4.0, help="orientation chamfer at the front-left corner (mm)")
    ap.add_argument("--seed", type=int, default=4)
    ap.add_argument("--profile", type=Path, default=s.DEFAULT_PROFILE)
    ap.add_argument("--out", type=Path, default=ROOT / "calibration_block")
    args = ap.parse_args()

    sx, sy, sz = args.size
    rng = np.random.default_rng(args.seed)
    run_dir = args.out / f"run_{datetime.now():%Y%m%d_%H%M%S}_{args.layout}"
    run_dir.mkdir(parents=True)

    outline = block_outline(args.size, args.notch)
    if args.layout == "grid":
        seeds, cells = grid_layout(outline, args.size, args.grid)
    else:
        layout = voronoi_layout_varied if args.order == "area" else voronoi_layout
        seeds, cells = layout(outline, args.n_cells, args.min_clear, rng)
    if args.infill:
        infill = args.infill
    elif args.order == "area":
        infill = assign_infill_by_area(cells, 10, 85)
    else:
        infill = assign_infill(len(cells), 10, 85, rng)
    if len(infill) != len(cells):
        raise SystemExit(f"--infill needs {len(cells)} values")
    args.infill = infill
    targets = [probe_target(c) for c in cells]
    block = s.prism_mesh(outline, 0.0, sz)
    # Pieces overlap their neighbours by 0.05 mm (PrusaSlicer trims overlaps, no internal walls) but are
    # clipped to the outline so the outside dimensions stay exact.
    pieces = [s.prism_mesh(c.buffer(0.05, join_style="mitre").intersection(outline).simplify(1e-9), 0.0, sz)
              for c in cells]
    names = [f"patch_{i}_{f:g}pct" for i, f in enumerate(args.infill)]
    block.export(run_dir / "calibration_block.stl")

    cfg = s.read_profile(args.profile)
    placement = (*(s.bed_center(cfg) - (sx / 2, sy / 2)), 0.0)
    model_3mf = run_dir / "calibration_block_prusaslicer.3mf"
    s.write_modifier_3mf(model_3mf, block, pieces, args.infill, placement, as_parts=True, name="calibration_block",
                         extra_settings={"fill_pattern": "rectilinear"}, names=names)

    exe = s.find_prusaslicer()
    version = s.check_slicer_version(exe)
    n = s.roundtrip_check(exe, model_3mf, run_dir, list(args.infill), "ModelPart")
    print(f"PrusaSlicer {version} kept all {n} patches with their fill_density")
    gcode = run_dir / "calibration_block_tpu_verified.gcode"
    s.slice_with_cli(exe, model_3mf, args.profile, gcode, run_dir / "slice.log")
    shutil.copy(args.profile, run_dir / args.profile.name)

    # Verification: every patch is the (clipped) Voronoi cell of its seed, so binning by nearest seed
    # uses exactly the patch boundaries the pieces were cut along.
    seeds3 = np.column_stack([seeds, np.zeros(len(seeds))])
    tree = cKDTree(seeds)
    region_of = lambda pts: tree.query(np.asarray(pts)[:, :2])[1]
    text = s.read_gcode_text(gcode)
    realised, info, sparse = s.measure_region_infill(text, block, region_of, len(cells))
    records = [{"id": i, "x": p.x, "y": p.y, "z": sz, "stiffness": float("nan")} for i, (p, _) in enumerate(targets)]
    target = np.asarray(args.infill, dtype=float)
    s.write_verification(records, target, realised, info, 5.0, run_dir)

    grid = s.voxelize(block, 1.0)
    labels = s.label_voxels(grid, seeds3)
    probe_xyz = np.array([[p.x, p.y, sz] for p, _ in targets])
    s.plot_toolpath_layers(sparse, region_of, grid, labels, target, info["layer_height"], probe_xyz,
                           run_dir / "calibration_block_toolpaths.png", cells2d=cells)
    plot_block_map(outline, cells, args.infill, targets, run_dir / "calibration_block_map.png",
                   f"Calibration block {sx:g} x {sy:g} x {sz:g} mm: probe each patch at the red cross")
    with open(run_dir / "probe_points.csv", "w", encoding="utf-8") as f:
        f.write("patch,infill_pct,probe_x_mm,probe_y_mm,clear_radius_mm\n")
        for i, (fi, (p, r)) in enumerate(zip(args.infill, targets)):
            f.write(f"{i},{fi:g},{p.x:.2f},{p.y:.2f},{r:.2f}\n")
    meta = {"size_mm": [sx, sy, sz], "layout": args.layout, "seed": args.seed, "infill_pct": args.infill,
            "pattern": "rectilinear", "infill_order": args.order, "patch_names": names, "seeds_mm": seeds.round(2).tolist(),
            "probe_points_mm": [[round(p.x, 2), round(p.y, 2)] for p, _ in targets],
            "min_clear_radius_mm": round(min(r for _, r in targets), 2),
            "orientation_mark": f"{args.notch:g} mm chamfer at the front-left corner (origin)",
            "slicer": f"PrusaSlicer {version}", "profile": args.profile.name,
            "verification": {k: v for k, v in info.items() if k != "fillable_area_mm2"},
            "realised_infill_pct": [round(float(r), 2) for r in realised]}
    (run_dir / "calibration_block_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # Bambu Studio project: same pieces, Bambu setting names (zig-zag = rectilinear).
    bambu_settings = [{"sparse_infill_density": f"{f:g}%", "sparse_infill_pattern": "zig-zag"} for f in args.infill]
    write_bambu_3mf(run_dir / "calibration_block_bambu.3mf", pieces, names, bambu_settings, (128.0, 128.0))
    # Display copy: no top skin, so the infill pattern of every patch is visible on the printed top face.
    write_bambu_3mf(run_dir / "calibration_block_bambu_open_top.3mf", pieces, names, bambu_settings, (128.0, 128.0),
                    title="calibration_block_open_top",
                    object_settings={"top_shell_layers": "0", "top_shell_thickness": "0"})
    print(f"Artifacts: {run_dir}")


if __name__ == "__main__":
    main()
