"""Reading scans and probe-point files in many formats, and bringing both into one frame (Z up, mm).

Meshes:  .obj .glb .gltf .ply .stl .off .3mf .dae, or a .zip containing one of them (Sketchfab / RealityScan
         downloads). .fbx .usdz .usd .abc .blend are recognised and explained (convert to .obj first).
Points:  .csv .tsv .txt (any delimiter), .json, .xlsx; flexible column names (x/X/x_mm, stiffness/k/...).
"""
import csv
import json
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import trimesh

MESH_TYPES = [".glb", ".gltf", ".obj", ".ply", ".stl", ".off", ".3mf", ".dae"]   # preference order inside a zip
UNSUPPORTED = {
    ".fbx": "FBX", ".usdz": "USDZ (Apple)", ".usd": "USD", ".usda": "USD", ".usdc": "USD",
    ".abc": "Alembic", ".blend": "Blender", ".max": "3ds Max", ".ma": "Maya", ".mb": "Maya",
}
Y_UP_TO_Z_UP = trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0])

# Accepted column names (lower-case, spaces/units stripped) for each field of a points file.
COLUMN_ALIASES = {
    "x": ["x", "x_mm", "xmm", "pos_x", "px"],
    "y": ["y", "y_mm", "ymm", "pos_y", "py"],
    "z": ["z", "z_mm", "zmm", "pos_z", "pz"],
    "stiffness": ["stiffness", "k", "k_n_mm", "k_nmm", "k_n_per_mm", "stiffness_n_mm", "stiffness_n_per_mm",
                  "stiffness_nmm", "kn_mm"],
}


class ScanError(ValueError):
    pass


# ---------------------------------------------------------------------------------------------- meshes
def _pick_from_zip(path, tmp):
    with zipfile.ZipFile(path) as z:
        z.extractall(tmp)
    files = [p for p in Path(tmp).rglob("*") if p.is_file()]
    for ext in MESH_TYPES:
        hits = sorted((p for p in files if p.suffix.lower() == ext), key=lambda p: -p.stat().st_size)
        if hits:
            return hits[0]
    found = sorted({p.suffix.lower() for p in files})
    raise ScanError(f"{path.name} contains no mesh this tool reads (found {found}). Supported: {MESH_TYPES}")


def read_mesh(path):
    """Geometry of a scan file as one triangle mesh in the file's own units and axes.

    Textures are dropped and the mesh is rebuilt from vertices + faces, which merges the duplicate vertices
    texture seams create (a textured scan otherwise loads as thousands of disconnected patches).
    """
    path = Path(path)
    if not path.exists():
        raise ScanError(f"File not found: {path}")
    ext = path.suffix.lower()
    if ext in UNSUPPORTED:
        raise ScanError(f"{UNSUPPORTED[ext]} files can't be read directly. Re-export the scan as .obj or .glb "
                        f"(RealityScan / Sketchfab download), or open it in Blender (free) and use "
                        f"File > Export > Wavefront (.obj).")
    if ext == ".zip":
        tmp = tempfile.mkdtemp(prefix="scan_zip_")
        inner = _pick_from_zip(path, tmp)
        print(f"  using {inner.name} from {path.name}")
        return read_mesh(inner)[0], inner.suffix.lower()
    if ext == ".dae":
        try:
            import collada  # noqa: F401
        except ImportError:
            raise ScanError("Reading .dae needs pycollada:  pip install pycollada  (or re-export as .obj)")
    if ext not in MESH_TYPES:
        raise ScanError(f"Unknown file type '{ext}'. Supported: {MESH_TYPES + ['.zip']}")
    loaded = trimesh.load(path, force="mesh", process=False)
    if isinstance(loaded, trimesh.PointCloud) or not hasattr(loaded, "faces") or len(loaded.faces) == 0:
        raise ScanError(f"{path.name} is a point cloud (no triangles). Export a *mesh* from the scanning app, "
                        f"not the point cloud.")
    mesh = trimesh.Trimesh(np.asarray(loaded.vertices), np.asarray(loaded.faces), process=True)
    return mesh, ext


def up_rotation(ext, up_axis):
    """4x4 rotation that makes the file's up axis +Z. auto: glTF is Y-up by spec, other formats Z-up."""
    up = ("y" if ext in (".glb", ".gltf") else "z") if up_axis == "auto" else up_axis
    if up not in ("y", "z"):
        raise ScanError(f"up_axis must be auto, y or z (got {up_axis})")
    return Y_UP_TO_Z_UP if up == "y" else np.eye(4)


def unit_factor(mesh, units):
    if units == "auto":
        return 1000.0 if float(mesh.extents.max()) < 5.0 else 1.0
    try:
        return {"m": 1000.0, "cm": 10.0, "mm": 1.0, "in": 25.4}[units]
    except KeyError:
        raise ScanError(f"units must be auto, m, cm, mm or in (got {units})")


# ---------------------------------------------------------------------------------------------- points
def _norm(name):
    n = name.strip().lower()
    for junk in ("(", ")", "[", "]", "/", "-", " ", "."):
        n = n.replace(junk, "_")
    return "_".join(p for p in n.split("_") if p)


def _rows_to_records(rows, source):
    if not rows:
        raise ScanError(f"{source}: no rows")
    cols = {_norm(c): c for c in rows[0].keys() if c is not None}
    pick = {}
    for field, aliases in COLUMN_ALIASES.items():
        hit = next((cols[a] for a in aliases if a in cols), None)
        if hit is None and field != "stiffness":
            raise ScanError(f"{source}: no '{field}' column (have {list(rows[0].keys())})")
        pick[field] = hit
    records = []
    for i, r in enumerate(rows):
        try:
            x, y, z = (float(r[pick[a]]) for a in "xyz")
        except (TypeError, ValueError):
            if all(not str(r.get(pick[a], "") or "").strip() for a in "xyz"):
                continue                                   # blank line
            raise ScanError(f"{source}: row {i + 2} has a non-numeric x/y/z")
        raw_k = str(r.get(pick["stiffness"], "") or "").strip() if pick["stiffness"] else ""
        k = float(raw_k) if raw_k else float("nan")
        records.append({"id": len(records), "x": x, "y": y, "z": z, "stiffness": k})
    return records


def read_points(path):
    """Probe points as records {id, x, y, z, stiffness}; stiffness is NaN where the file leaves it blank."""
    path = Path(path)
    if not path.exists():
        raise ScanError(f"Points file not found: {path}")
    ext = path.suffix.lower()
    if ext == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        rows = data.get("points", data) if isinstance(data, dict) else data
        if rows and isinstance(rows[0], (list, tuple)):
            rows = [dict(zip(["x", "y", "z", "stiffness"], r)) for r in rows]
        return _rows_to_records([{k: ("" if v is None else v) for k, v in r.items()} for r in rows], path.name)
    if ext in (".xlsx", ".xls"):
        try:
            import pandas as pd
        except ImportError:
            raise ScanError("Reading Excel files needs pandas + openpyxl:  pip install pandas openpyxl  "
                            "(or save the sheet as .csv)")
        try:
            df = pd.read_excel(path).fillna("")
        except ImportError as e:
            raise ScanError(f"Reading Excel needs a working openpyxl ({e}). Run  pip install -U openpyxl  "
                            f"or save the sheet as .csv")
        return _rows_to_records(df.astype(str).to_dict("records"), path.name)
    text = path.read_text(encoding="utf-8-sig")
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t |")
    except csv.Error:
        dialect = csv.excel
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if lines and all(c.replace(".", "").replace("-", "").replace("e", "").isdigit() or not c.strip()
                     for c in lines[0].replace(dialect.delimiter, " ").split()):
        lines.insert(0, dialect.delimiter.join(["x", "y", "z", "stiffness"]))     # no header: assume x y z k
    reader = csv.DictReader(lines, dialect=dialect, skipinitialspace=True)
    return _rows_to_records(list(reader), path.name)


def xyz_of(records):
    return np.array([[r["x"], r["y"], r["z"]] for r in records], dtype=float)


def place_points(records, rotation, mm_factor, scale, mesh_mm, units="auto"):
    """Move points exactly like the mesh (rotation, then to mm, then print scale).

    units: "mm" (picker output), "m" / file units, or "auto": try both and keep whichever lands on the mesh.
    Returns (records, chosen_factor, distances_mm).
    """
    raw = xyz_of(records)
    candidates = {"mm": 1.0, "m": 1000.0, "cm": 10.0, "in": 25.4, "file": mm_factor}
    tries = [units] if units != "auto" else ["mm", "file", "m"]
    best = None
    for u in tries:
        if u not in candidates:
            raise ScanError(f"points units must be auto, mm, m, cm, in or file (got {u})")
        xyz = trimesh.transform_points(raw * candidates[u], rotation) * scale
        _, d, _ = trimesh.proximity.closest_point(mesh_mm, xyz)
        if best is None or np.median(d) < np.median(best[2]):
            best = (u, xyz, d)
    u, xyz, d = best
    out = [{**r, "x": float(p[0]), "y": float(p[1]), "z": float(p[2])} for r, p in zip(records, xyz)]
    return out, u, d
