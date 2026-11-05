"""Pick probe points on the scan (Open3D window) and save them for scan_to_print.py.

  python pick_points.py scans/my_object.obj          # writes scans/my_object_points.csv

Shift + left-click on the surface to pick a point, then close the window (extra points can be deleted
from the CSV). Coordinates are saved in mm, in the scan file's own axes (x1000 if the file is in metres, as
RealityScan and glTF exports are). The stiffness column is left blank: fill it in with each point's measured
k (N/mm) after probing, then run scan_to_print.py with --points.

Accepts every scan type scan_to_print.py does (.obj .glb .gltf .ply .stl .off .3mf .zip ...).
Based on the team's 3d_check.py vertex picker.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

import scan_io as io


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", type=Path)
    ap.add_argument("--out", type=Path, default=None, help="default: <scan name>_points.csv next to the scan")
    ap.add_argument("--units", choices=["auto", "m", "cm", "mm", "in"], default="auto",
                    help="scan file units (auto: under 5 across = metres)")
    args = ap.parse_args()
    try:
        import open3d as o3d
    except ImportError:
        raise SystemExit("The point picker needs Open3D:  pip install open3d")

    try:
        mesh, _ = io.read_mesh(args.scan)
    except io.ScanError as e:
        raise SystemExit(str(e))
    to_mm = io.unit_factor(mesh, args.units)
    om = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(np.asarray(mesh.vertices)),
                                   o3d.utility.Vector3iVector(np.asarray(mesh.faces)))
    om.compute_vertex_normals()
    om.paint_uniform_color([0.75, 0.75, 0.78])

    vis = o3d.visualization.VisualizerWithVertexSelection()
    vis.create_window(window_name="Tactile Setup - Point Picker (Shift+click to pick, close when done)",
                      width=1280, height=860)
    vis.add_geometry(om)
    print("Shift + left-click on the surface to pick probe points; close the window when finished.")
    vis.run()
    vis.destroy_window()

    picked = vis.get_picked_points()
    if not picked:
        raise SystemExit("No points picked.")
    out = args.out or args.scan.with_name(f"{args.scan.stem}_points.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "x", "y", "z", "stiffness"])
        for i, p in enumerate(picked):
            x, y, z = np.asarray(p.coord) * to_mm
            w.writerow([i, f"{x:.4f}", f"{y:.4f}", f"{z:.4f}", ""])
            print(f"Point {i} -> vertex #{p.index} | X={x:.3f} Y={y:.3f} Z={z:.3f} mm")
    print(f"\nSaved {len(picked)} points to {out}\nNext: probe each point, type its stiffness (N/mm) into the "
          f"'stiffness' column, then run  python scan_to_print.py {args.scan} --points {out}")


if __name__ == "__main__":
    main()
