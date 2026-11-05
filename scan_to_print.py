"""Scan + probed stiffness points -> a print whose infill follows the measured stiffness.

  python scan_to_print.py scans/my_object.obj --points scans/my_object_points.csv
  python scan_to_print.py scans/my_object.glb --stl-only                    # just a clean, printable STL
  python scan_to_print.py examples/test_block.obj --points examples/test_block_points.csv   # try it

All settings live in config.toml (printer, filament, shells, infill pattern, stiffness->infill mapping...).
Command-line flags override the config for one run. See README.md for the full walkthrough.

Steps
  1. read the scan (any common mesh format or a .zip download), fix axes and units, repair it watertight
  2. read the probe points (csv/tsv/txt/json/xlsx), put them in the same frame, check they sit on the surface
  3. map stiffness -> infill and split the part into one Voronoi region per point
  4. write the outputs you ask for: Bambu project, Bambu sliced file, PrusaSlicer project, STL(s)
  5. (Bambu sliced) read the G-code back and measure the infill of every region
"""
import argparse
import csv
import json
import sys
import tomllib
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import bambu_profiles as bp
import make_calibration_block as mcb
import scan_io as io
import scan_to_infill as s

ROOT = Path(__file__).resolve().parent
FORMATS = ["bambu_project", "bambu_sliced", "prusa_project", "stl", "region_stls"]


# ------------------------------------------------------------------------------------------- config
def load_config(path, args):
    cfg = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    sec = lambda name: cfg.setdefault(name, {})
    over = {  # command-line flag -> (section, key)
        "units": ("scan", "units"), "up": ("scan", "up_axis"), "scale": ("scan", "scale"),
        "fit_largest_mm": ("scan", "fit_largest_mm"), "max_faces": ("scan", "max_faces"),
        "points_units": ("points", "units"), "surface_tol": ("points", "surface_tolerance_mm"),
        "mapping": ("mapping", "mode"), "voxel_mm": ("regions", "voxel_mm"), "partition": ("regions", "partition"),
        "out": ("output", "folder"), "printer": ("bambu", "printer"), "process": ("bambu", "process"),
        "filament": ("bambu", "filament"), "plate": ("bambu", "plate"), "bambu_exe": ("bambu", "executable"),
        "top_shells": ("print", "top_shell_layers"), "bottom_shells": ("print", "bottom_shell_layers"),
        "walls": ("print", "wall_loops"), "infill_pattern": ("print", "infill_pattern"),
        "top_pattern": ("print", "top_surface_pattern"),
    }
    for flag, (section, key) in over.items():
        v = getattr(args, flag, None)
        if v is not None:
            sec(section)[key] = v
    if args.formats:
        sec("output")["formats"] = [f.strip() for f in args.formats.split(",") if f.strip()]
    if args.stl_only:
        sec("output")["formats"] = ["stl"]
    if args.no_verify:
        sec("output")["verify"] = False
    bad = set(cfg["output"]["formats"]) - set(FORMATS)
    if bad:
        raise SystemExit(f"Unknown output format(s) {sorted(bad)}; choose from {FORMATS}")
    return cfg


def make_mapping(m):
    """Stiffness (N/mm) -> infill (%) from the [mapping] section."""
    if m["mode"] == "linear":
        k0, k1, f0, f1 = float(m["k_min"]), float(m["k_max"]), float(m["infill_min"]), float(m["infill_max"])
        return lambda k: f0 + np.clip((np.asarray(k, float) - k0) / (k1 - k0), 0, 1) * (f1 - f0)
    if m["mode"] == "table":
        tk, tf = np.asarray(m["table_k"], float), np.asarray(m["table_infill"], float)
        if len(tk) != len(tf) or len(tk) < 2 or np.any(np.diff(tk) <= 0):
            raise SystemExit("[mapping] table_k must be increasing and the same length as table_infill")
        return lambda k: np.interp(np.asarray(k, float), tk, tf)
    raise SystemExit(f"[mapping] mode must be 'linear' or 'table' (got {m['mode']})")


def bambu_process_overrides(cfg):
    p = cfg["print"]
    out = {"top_shell_layers": p["top_shell_layers"], "bottom_shell_layers": p["bottom_shell_layers"],
           # the *_shell_thickness keys are minimums that override the layer counts above, so zero them
           "top_shell_thickness": 0, "bottom_shell_thickness": 0, "wall_loops": p["wall_loops"],
           "top_surface_pattern": bp.bambu_pattern(p["top_surface_pattern"]),
           "bottom_surface_pattern": bp.bambu_pattern(p["bottom_surface_pattern"]),
           "sparse_infill_pattern": bp.bambu_pattern(p["infill_pattern"])}
    out.update(cfg["bambu"].get("process_overrides", {}))
    return out


# ------------------------------------------------------------------------------------------- steps
def prepare_mesh(path, sc):
    raw, ext = io.read_mesh(path)
    rot = io.up_rotation(ext, sc["up_axis"])
    mm = io.unit_factor(raw, sc["units"])
    mesh = raw.copy()
    mesh.apply_transform(rot)
    mesh.apply_scale(mm)
    scale = float(sc["scale"])
    if float(sc.get("fit_largest_mm", 0) or 0) > 0:
        scale = float(sc["fit_largest_mm"]) / float(mesh.extents.max())
    mesh.apply_scale(scale)
    print(f"  {Path(path).name}: {len(raw.faces)} triangles, units x{mm:g} -> mm, print scale x{scale:.4g}, "
          f"up axis {'y' if rot[1, 2] else 'z'} -> {np.round(mesh.extents, 1)} mm")
    if len(mesh.faces) > int(sc["max_faces"]):
        n0 = len(mesh.faces)
        mesh = mesh.simplify_quadric_decimation(face_count=int(sc["max_faces"]))
        print(f"  simplified {n0} -> {len(mesh.faces)} triangles")
    mesh = s.repair_mesh(mesh)
    return mesh, {"rotation": rot, "mm_factor": mm, "scale": scale, "ext": ext}


def prepare_points(path, frame, mesh, pc, allow_missing_k=False):
    records = io.read_points(path)
    records, units, dist = io.place_points(records, frame["rotation"], frame["mm_factor"], frame["scale"], mesh,
                                           pc["units"])
    print(f"  {len(records)} points ({units} units): distance to surface median {np.median(dist):.2f} mm, "
          f"max {dist.max():.2f} mm")
    far = np.flatnonzero(dist > float(pc["surface_tolerance_mm"]))
    if far.size:
        raise SystemExit(f"Points {far.tolist()} are {np.round(dist[far], 1).tolist()} mm off the mesh surface. "
                         f"Was the points file picked on this exact mesh file? Try --points-units mm or m.")
    missing = [r["id"] for r in records if not np.isfinite(r["stiffness"])]
    if missing and not allow_missing_k:
        raise SystemExit(f"Points {missing} have no stiffness value. Fill in the stiffness column (N/mm) after "
                         f"probing, or use --stl-only for an uncalibrated print.")
    return records


def write_region_stls(run_dir, pieces, names, keep, target, records):
    folder = run_dir / "region_stls"
    folder.mkdir()
    rows = []
    for piece, name, i in zip(pieces, names, keep):
        f = folder / f"{name}.stl"
        piece.export(f)
        r = records[i]
        rows.append({"file": f.name, "region": i, "infill_pct": target[i], "stiffness_n_per_mm": r["stiffness"],
                     "probe_x_mm": round(r["x"], 2), "probe_y_mm": round(r["y"], 2), "probe_z_mm": round(r["z"], 2)})
    with open(folder / "infill_map.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    (folder / "HOW_TO_IMPORT.txt").write_text(
        "These STLs are the pieces of one part, already in the right place relative to each other.\n\n"
        "Bambu Studio / OrcaSlicer: select ALL the .stl files and drag them onto the plate together. When asked\n"
        "'Load these files as a single object with multiple parts?' choose Yes. Then click each part in the\n"
        "object list, open its settings and set 'Sparse infill density' to the value in infill_map.csv.\n\n"
        "PrusaSlicer: add the first STL, then right-click it > Add part > Load... for the others. Set each part's\n"
        "infill density (right-click > Add settings > Infill) from infill_map.csv.\n\n"
        "Cura: load them all, select all, Ctrl+Alt+G (merge), then use 'Per Model Settings' on each piece.\n",
        encoding="utf-8")
    return folder


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", type=Path, nargs="?", help="scan file: .obj .glb .gltf .ply .stl .off .3mf .dae or .zip")
    ap.add_argument("--points", type=Path, help="probe points: .csv .tsv .txt .json .xlsx with x, y, z, stiffness")
    ap.add_argument("--config", type=Path, default=ROOT / "config.toml")
    ap.add_argument("--formats", help=f"comma list from {FORMATS} (overrides [output] formats)")
    ap.add_argument("--stl-only", action="store_true", help="just write the cleaned watertight STL")
    ap.add_argument("--demo-points", type=int, default=0,
                    help="no probe data yet: put N reachable points with random stiffness on the scan (testing)")
    ap.add_argument("--no-verify", action="store_true", help="skip reading the sliced G-code back")
    ap.add_argument("--list-bambu-profiles", action="store_true", help="print Bambu printer/process/filament names")
    g = ap.add_argument_group("overrides for config.toml")
    g.add_argument("--units", choices=["auto", "m", "cm", "mm", "in"])
    g.add_argument("--up", choices=["auto", "y", "z"])
    g.add_argument("--scale", type=float)
    g.add_argument("--fit-largest-mm", type=float)
    g.add_argument("--max-faces", type=int)
    g.add_argument("--points-units", choices=["auto", "mm", "m", "cm", "in", "file"])
    g.add_argument("--surface-tol", type=float)
    g.add_argument("--mapping", choices=["linear", "table"])
    g.add_argument("--voxel-mm", type=float)
    g.add_argument("--partition", choices=["auto", "columns", "3d"])
    g.add_argument("--out")
    g.add_argument("--printer")
    g.add_argument("--process")
    g.add_argument("--filament")
    g.add_argument("--plate")
    g.add_argument("--bambu-exe")
    g.add_argument("--top-shells", type=int)
    g.add_argument("--bottom-shells", type=int)
    g.add_argument("--walls", type=int)
    g.add_argument("--infill-pattern")
    g.add_argument("--top-pattern")
    args = ap.parse_args(argv)
    cfg = load_config(args.config, args)

    if args.list_bambu_profiles:
        exe = bp.find_executable(cfg["bambu"].get("executable", ""))
        profs = bp.list_profiles(exe, cfg["bambu"]["printer"])
        print(f"Bambu Studio {bp.version(exe)} at {exe}\n\nPrinters:")
        print("\n".join(f"  {n}" for n in profs["machine"] if "nozzle" in n))
        print(f"\nProcesses for {cfg['bambu']['printer']}:\n" + "\n".join(f"  {n}" for n in profs["process"]))
        print(f"\nFilaments for {cfg['bambu']['printer']}:\n" + "\n".join(f"  {n}" for n in profs["filament"]))
        return
    if args.scan is None:
        ap.error("give a scan file (or --list-bambu-profiles)")

    formats = cfg["output"]["formats"]
    needs_regions = any(f in formats for f in ("bambu_project", "bambu_sliced", "prusa_project", "region_stls"))
    if needs_regions and not args.points and not args.demo_points:
        print("No --points given: writing the cleaned STL only (add --points for stiffness-matched infill).")
        formats, needs_regions = ["stl"], False

    stem = args.scan.stem
    summary = {"scan": str(args.scan), "points": str(args.points) if args.points else None, "formats": formats,
               "outputs": {}}

    print("[1] scan")
    try:
        mesh, frame = prepare_mesh(args.scan, cfg["scan"])
    except io.ScanError as e:
        raise SystemExit(str(e))

    records = []
    if needs_regions:
        print("[2] probe points")
        if args.points:
            try:
                records = prepare_points(args.points, frame, mesh, cfg["points"])
            except io.ScanError as e:
                raise SystemExit(str(e))
        else:
            records = s.sample_stiffness_points(mesh, args.demo_points, np.random.default_rng(0),
                                                k_range=(cfg["mapping"]["k_min"], cfg["mapping"]["k_max"]),
                                                min_spacing=max(5.0, float(mesh.extents[:2].min()) / 6))
            print(f"  DEMO: {len(records)} random reachable points with random stiffness")
    mesh, records, shift = s.to_part_frame(mesh, records)
    summary["part_size_mm"] = np.round(mesh.extents, 2).tolist()
    # Output folder only once the inputs are known good, so failed runs leave nothing behind.
    out_root = Path(cfg["output"]["folder"])
    run_dir = (out_root if out_root.is_absolute() else ROOT / out_root) / f"{stem}_{datetime.now():%Y%m%d_%H%M%S}"
    run_dir.mkdir(parents=True)
    if "stl" in formats:
        f = run_dir / f"{stem}_clean.stl"
        mesh.export(f)
        summary["outputs"]["stl"] = f.name
    if not needs_regions:
        _finish(run_dir, summary)
        return

    print("[3] stiffness -> infill regions")
    s.save_points_csv(records, run_dir / "points_in_part_frame.csv")
    reg_args = SimpleNamespace(voxel_pitch=float(cfg["regions"]["voxel_mm"]) or None,
                               partition=cfg["regions"]["partition"])
    regions, target, region_solids = s.run_stage_b(reg_args, mesh, records, run_dir,
                                                   to_infill=make_mapping(cfg["mapping"]))
    parts = s.split_into_parts(mesh, region_solids)
    keep = [i for i, p in enumerate(parts) if p is not None]
    pieces = [parts[i] for i in keep]
    names = [f"R{i:02d}_{target[i]:g}pct" for i in keep]
    summary.update(partition=regions["mode"], regions=len(keep), target_infill_pct=target.tolist())

    if "region_stls" in formats:
        summary["outputs"]["region_stls"] = write_region_stls(run_dir, pieces, names, keep, target, records).name

    if "prusa_project" in formats:
        prof = ROOT / cfg["prusa"]["profile"]
        placement = (*(s.bed_center(s.read_profile(prof)) - (mesh.bounds[0, :2] + mesh.bounds[1, :2]) / 2), 0.0)
        f = run_dir / f"{stem}_prusaslicer.3mf"
        s.write_modifier_3mf(f, mesh, pieces, [target[i] for i in keep], placement, as_parts=True, name=stem,
                             extra_settings={"fill_pattern": cfg["print"]["infill_pattern"]}, names=names)
        summary["outputs"]["prusa_project"] = f.name

    if "bambu_project" in formats or "bambu_sliced" in formats:
        print("[4] Bambu Studio")
        b = cfg["bambu"]
        exe = bp.find_executable(b.get("executable", ""))
        overrides = bambu_process_overrides(cfg)
        profiles = bp.flat_profiles(exe, run_dir / "bambu_profiles", b["printer"], b["process"], b["filament"],
                                    b["plate"], process_overrides=overrides,
                                    filament_overrides=b.get("filament_overrides", {}))
        settings = bp.project_settings(exe, profiles, run_dir / "bambu_profiles" / "export_settings.log")
        pattern = bp.bambu_pattern(cfg["print"]["infill_pattern"])
        part_settings = [{"sparse_infill_density": f"{target[i]:g}%", "sparse_infill_pattern": pattern} for i in keep]
        project = run_dir / f"{stem}_bambu.3mf"
        mcb.write_bambu_3mf(project, pieces, names, part_settings, tuple(b["bed_centre_mm"]), title=stem,
                            project_settings=settings)
        summary["outputs"]["bambu_project"] = project.name
        summary["bambu"] = {"studio_version": bp.version(exe), "printer": b["printer"], "process": b["process"],
                            "filament": b["filament"], "plate": b["plate"], "process_overrides": overrides}
        print(f"  {project.name}: {len(pieces)} parts, {b['printer']}, {b['filament']}, {b['plate']}")
        if "bambu_sliced" in formats:
            sliced = bp.slice_3mf(exe, project, profiles, run_dir, f"{stem}_bambu.gcode.3mf",
                                  run_dir / "bambu_profiles" / "slice.log")
            summary["outputs"]["bambu_sliced"] = sliced.name
            print(f"  sliced -> {sliced.name}")
            if cfg["output"].get("verify", True):
                print("[5] verify infill in the G-code")
                text = s.read_gcode_text(sliced)
                realised, info, sparse = s.measure_region_infill(text, mesh, regions["region_of"], len(records))
                rows = s.write_verification(records, target, realised, info, 5.0, run_dir)
                err = np.array([r["error_pp"] for r in rows], dtype=float)
                summary["verification"] = {"mean_abs_error_pct_points": round(float(np.nanmean(np.abs(err))), 2),
                                           "max_abs_error_pct_points": round(float(np.nanmax(np.abs(err))), 2)}
    _finish(run_dir, summary)


def _finish(run_dir, summary):
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nDone. Output folder: {run_dir}")
    labels = {"bambu_sliced": "ready to print (open in Bambu Studio > Print, or copy to the SD card)",
              "bambu_project": "open in Bambu Studio, check, Slice, Print",
              "prusa_project": "open in PrusaSlicer", "stl": "clean watertight model, any slicer",
              "region_stls": "one STL per region + infill_map.csv (see HOW_TO_IMPORT.txt)"}
    for k, v in summary["outputs"].items():
        print(f"  {v:<45} {labels[k]}")


if __name__ == "__main__":
    sys.exit(main())
