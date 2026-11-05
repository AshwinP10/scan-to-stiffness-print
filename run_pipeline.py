"""Stiffness-to-infill proof of concept: dummy part + random probe points -> graded-infill 3MF ->
PrusaSlicer CLI -> G-code -> per-region infill verification.

Usage:  python run_pipeline.py [--n-points 10] [--seed 7] [--size 40 40 10]
"""
import argparse
import csv
import json
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np

from infill_poc import mapping, regions, slicer, threemf, verify, viz

ROOT = Path(__file__).resolve().parent
DEFAULT_PROFILE = ROOT / "profiles" / "prusaslicer_2.5.0_poc.ini"


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


def flow_spacing(width, height):
    return width - height * (1 - np.pi / 4)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-points", type=int, default=10)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--size", type=float, nargs=3, default=(40.0, 40.0, 10.0), metavar=("X", "Y", "Z"))
    ap.add_argument("--k-range", type=float, nargs=2, default=(mapping.K_MIN, mapping.K_MAX))
    ap.add_argument("--min-spacing", type=float, default=8.0, help="minimum distance between probe points (mm)")
    ap.add_argument("--tolerance", type=float, default=5.0, help="allowed |measured - target| in percentage points")
    ap.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    ap.add_argument("--slicer", default=None, help="path to prusa-slicer-console")
    ap.add_argument("--allow-version-mismatch", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "out")
    args = ap.parse_args()

    sx, sy, sz = args.size
    run_dir = args.out / f"run_{datetime.now():%Y%m%d_%H%M%S}_seed{args.seed}"
    run_dir.mkdir(parents=True)
    cfg = read_profile(args.profile)
    shutil.copy(args.profile, run_dir / args.profile.name)
    rng = np.random.default_rng(args.seed)

    print("[0] Inputs: dummy part + random probe points")
    part = regions.make_box_model(sx, sy, sz)
    xy, k = regions.sample_probe_points(args.n_points, sx, sy, args.k_range, args.min_spacing, rng)

    print("[1] Stiffness -> infill mapping")
    target = np.round(mapping.stiffness_to_infill(k, *args.k_range), 1)

    print("[2-3] Voronoi partition + modifier extrusion")
    cells, modifier_meshes = regions.build_regions(xy, sx, sy, sz)

    print("[4] Visualising regions")
    viz.plot_voronoi_map(cells, xy, k, target, sx, sy, run_dir / "voronoi_infill_map.png")
    viz.plot_modifiers_3d(cells, target, sx, sy, sz, run_dir / "modifiers_3d.png")

    print("[5] Writing 3MF with per-modifier fill_density")
    placement = (*(bed_center(cfg) - (sx / 2, sy / 2)), 0.0)
    model_3mf = run_dir / "graded_infill.3mf"
    mods = [(f"region_{i:02d}", m, {"fill_density": f"{f:g}%"}) for i, (m, f) in enumerate(zip(modifier_meshes, target))]
    threemf.write_prusa_3mf(model_3mf, part, mods, placement)

    exe = slicer.find_prusaslicer(args.slicer)
    version = slicer.check_pinned_version(exe, args.allow_version_mismatch)
    print(f"    PrusaSlicer {version}: round-trip check of the 3MF")
    rt_3mf = run_dir / "roundtrip_by_prusaslicer.3mf"
    slicer.roundtrip_3mf(exe, model_3mf, rt_3mf, run_dir / "prusaslicer_roundtrip.log")
    rt_mods = [v for v in threemf.read_volume_settings(rt_3mf) if v.get("volume_type") == "ParameterModifier"]
    rt_density = sorted(float(v.get("fill_density", "nan").rstrip("%")) for v in rt_mods)
    if len(rt_mods) != len(mods) or not np.allclose(rt_density, sorted(target)):
        raise RuntimeError(f"PrusaSlicer did not preserve the modifiers: got {len(rt_mods)} with densities {rt_density}")
    print(f"    OK: PrusaSlicer re-exported all {len(rt_mods)} modifiers with their fill_density overrides")

    print("[6] Slicing via CLI")
    gcode = run_dir / "graded_infill.gcode"
    slicer.slice_to_gcode(exe, model_3mf, args.profile, gcode, run_dir / "prusaslicer_slice.log")

    print("[7] Verifying infill per region from G-code")
    segs = verify.parse_extrusions(gcode)
    expected_offset = np.array(placement[:2])
    offset = verify.measured_offset(segs[verify.EXT_PERIMETER], sx, sy)
    offset_err = float(np.hypot(*(offset - expected_offset)))
    if offset_err > 0.5:
        print(f"    WARNING: part sits {offset_err:.2f} mm from the requested placement; using measured offset")
    layer_h = float(cfg["layer_height"])
    measured, info = verify.region_densities(segs[verify.INFILL], cells, offset, layer_h,
                                             float(cfg["filament_diameter"]))
    err = measured - target

    spacing = flow_spacing(float(cfg["infill_extrusion_width"]), layer_h)
    rows = []
    for i in range(len(target)):
        pitch = spacing / (target[i] / 100)
        rows.append({
            "region": f"R{i}", "x_mm": round(xy[i, 0], 3), "y_mm": round(xy[i, 1], 3),
            "stiffness_N_per_mm": round(k[i], 4), "target_infill_pct": target[i],
            "measured_infill_pct": round(measured[i], 2), "error_pp": round(err[i], 2),
            "cell_area_mm2": round(cells[i].area, 2), "fillable_area_mm2": round(info["fill_areas"][i], 2),
            "expected_line_pitch_mm": round(pitch, 2),
            "cell_width_in_pitches": round(np.sqrt(info["fill_areas"][i]) / pitch, 1),
            "within_tolerance": bool(abs(err[i]) <= args.tolerance),
        })
    with open(run_dir / "verification_report.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    viz.plot_verification(target, measured, args.tolerance, run_dir / "verification.png")
    viz.plot_layer_toolpaths(cells, target, segs[verify.INFILL], offset, sx, sy, run_dir / "gcode_toolpaths.png")

    summary = {
        "seed": args.seed, "n_points": args.n_points, "part_size_mm": [sx, sy, sz],
        "k_range_N_per_mm": list(args.k_range),
        "mapping": f"linear {mapping.INFILL_MIN:g}-{mapping.INFILL_MAX:g}% over k_range (placeholder)",
        "slicer": f"PrusaSlicer {version}", "slicer_exe": exe, "profile": args.profile.name,
        "placement_mm": list(placement), "measured_offset_mm": offset.round(3).tolist(),
        "placement_error_mm": round(offset_err, 3),
        "n_infill_layers": info["n_infill_layers"], "n_infill_segments": info["n_segments"],
        "n_segments_split_across_regions": info["n_split_segments"],
        "fill_rect_mm": [round(v, 3) for v in info["fill_rect"]],
        "mean_abs_error_pp": round(float(np.abs(err).mean()), 2),
        "max_abs_error_pp": round(float(np.abs(err).max()), 2),
        "tolerance_pp": args.tolerance,
        "regions_within_tolerance": int(sum(r["within_tolerance"] for r in rows)),
    }
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\n{'region':>6} {'k':>6} {'target%':>8} {'measured%':>10} {'err pp':>7} {'cell mm2':>9} {'pitches':>8}")
    for r in rows:
        flag = "" if r["within_tolerance"] else "  <-- outside tolerance"
        print(f"{r['region']:>6} {r['stiffness_N_per_mm']:>6.2f} {r['target_infill_pct']:>8.1f} "
              f"{r['measured_infill_pct']:>10.2f} {r['error_pp']:>+7.2f} {r['cell_area_mm2']:>9.1f} "
              f"{r['cell_width_in_pitches']:>8.1f}{flag}")
    print(f"\nMAE {summary['mean_abs_error_pp']} pp, max {summary['max_abs_error_pp']} pp, "
          f"{summary['regions_within_tolerance']}/{len(rows)} regions within +/-{args.tolerance:g} pp")
    print(f"Artifacts: {run_dir}")


if __name__ == "__main__":
    main()
