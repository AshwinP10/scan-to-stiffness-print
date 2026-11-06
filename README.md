# Scan → Stiffness → Print

Turn a 3D scan of a real object into a 3D print whose **squishiness matches the original, spot by spot**.

You scan the object (RealityScan), press a stiffness probe into it at chosen points, and this code builds a
print file where every region of the part gets its own infill density: stiff spots print dense, soft spots
print sparse. The default target is a **Bambu Lab P1S (0.4 mm nozzle) printing Generic TPU on the Textured PEI
plate**, with thin skins (2 top / 2 bottom layers) so the infill dominates how the print feels.

```
 phone scan          laptop                                                     printer
┌───────────┐   ┌──────────────┐   ┌─────────────┐   ┌──────────────────────┐   ┌──────────────┐
│RealityScan│──►│ scan file    │──►│pick_points  │──►│ scan_to_print.py     │──►│ Bambu P1S    │
│  capture  │   │ .obj / .glb  │   │ .csv points │   │ regions + infill     │   │ TPU, Textured│
└───────────┘   └──────────────┘   └──────┬──────┘   │ -> .3mf / .gcode.3mf │   │ PEI plate    │
                                          │ probe    │    / .stl            │   └──────────────┘
                                          ▼ each pt  └──────────────────────┘
                                   stiffness (N/mm) typed into the .csv
```

**Two kinds of output, pick either:**

| You want | Use this output | What you do with it |
|---|---|---|
| A finished file, ready to print | `<name>_bambu.gcode.3mf` | Open in Bambu Studio → **Print** (or copy to the SD card) |
| A file you slice and tweak yourself | `<name>_bambu.3mf` (Bambu Studio), `<name>_prusaslicer.3mf` (PrusaSlicer) | Open → check → **Slice** → Print |
| A plain STL | `<name>_clean.stl` | Any slicer, one infill everywhere |
| STLs for any other slicer, with varying infill | `region_stls/` (one STL per region + `infill_map.csv`) | Import all as one object, set each part's infill (see `HOW_TO_IMPORT.txt` inside) |

> **Why isn't "varying infill" just one STL?** An STL only stores the shape. It has no place to say "this
> part is 20 % infill, that part is 70 %". The `.3mf` files carry the shape *and* each region's infill; the
> `region_stls/` folder carries the pieces plus a table of their infill for slicers that need you to type it.

---

## Contents
1. [Install (once)](#1-install-once)
2. [Try it on the example (2 minutes)](#2-try-it-on-the-example-2-minutes)
3. [Scan the object with RealityScan](#3-scan-the-object-with-realityscan)
4. [Move the scan to the laptop](#4-move-the-scan-to-the-laptop)
5. [Pick the probe points](#5-pick-the-probe-points)
6. [Probe and fill in stiffness](#6-probe-and-fill-in-stiffness)
7. [Build the print](#7-build-the-print)
8. [Print it](#8-print-it)
9. [Changing settings](#9-changing-settings)
10. [Supported file types](#10-supported-file-types)
11. [Calibration block](#11-calibration-block-measure-the-real-stiffness-vs-infill-curve)
12. [How it works / where to change the code](#12-how-it-works--where-to-change-the-code)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Install (once)

**You need:** Python 3.11 or newer, and Bambu Studio (free, from bambulab.com) if you want the Bambu files.

```bash
git clone <this repo's URL>
cd <repo folder>
python -m venv .venv
```
Activate the environment, then install:
```bash
# Windows (PowerShell)
.venv\Scripts\Activate.ps1
# macOS / Linux
source .venv/bin/activate

pip install -r requirements.txt
```

Check it works:
```bash
python scan_to_print.py --help
python scan_to_print.py --list-bambu-profiles
```
The second command should print your Bambu Studio version and the printer / process / filament names it found.
If it says Bambu Studio is not installed, either install it, or set its path in `config.toml` under
`[bambu] executable`, or skip the Bambu outputs (`--formats stl,region_stls,prusa_project`).

| Platform | Where the code looks for Bambu Studio |
|---|---|
| Windows | `C:\Program Files\Bambu Studio\bambu-studio.exe` |
| macOS | `/Applications/BambuStudio.app/Contents/MacOS/BambuStudio` |
| Linux | `bambu-studio` on your PATH |

> Tested with Bambu Studio **01.07.07.89** on Windows. A newer Bambu Studio should work (the settings are
> exported from whatever version you have), but run the example once after upgrading.

---

## 2. Try it on the example (2 minutes)

The repo ships a small test part (`examples/test_block.obj`, a stepped block saved in metres like a
RealityScan export) and six probe points with stiffness values.

```bash
python scan_to_print.py examples/test_block.obj --points examples/test_block_points.csv
```

You get a folder like `out_print/test_block_20261105_143931/` containing:

```
test_block_clean.stl            the repaired, watertight part
region_stls/                    R00_13.2pct.stl ... one STL per region, infill_map.csv, HOW_TO_IMPORT.txt
test_block_bambu.3mf            Bambu Studio project (slice it yourself)
test_block_bambu.gcode.3mf      sliced for the P1S, ready to print
stageD_verification.csv/.png    infill measured back out of the G-code, per region
stageB_regions_3d.png           picture of the regions, coloured by infill
stageB_regions_colored.glb      3D view of the regions (opens in Windows 3D Viewer / any glTF viewer)
summary.json                    everything that was done, with the exact settings
```

The terminal ends with a verification table: target infill vs. the infill actually in the G-code, for every
region. On the example every region lands within about 0.6 % points of its target.

---

## 3. Scan the object with RealityScan

The **RealityScan** app (Epic Games, free, iOS and Android) turns photos into a textured 3D model.

**Capture tips** (from Epic's guide, plus what worked for us):
- Walk around the object in **circles and arches at several heights**; keep the object filling at least half
  of each photo. Every spot you want in the model must appear in at least two photos.
- **Small steps:** neighbouring photos should overlap about 70 %; don't swing more than ~30° between shots.
- Sharp, steady photos; even, diffuse light (no hard shadows); a plain, non-shiny background.
- Shiny, transparent or totally plain surfaces scan badly. A light dusting of matte spray or painter's tape
  with a few marks helps.
- The app needs at least 20 photos to start and allows up to 300.
- **Turntable rig:** if the object sits on our turntable, put markers (the turntable crosshair / AprilTags) in
  view. They give the scan a fixed origin, and on the desktop version they can give it real-world scale.

**Process and export (phone app):**
1. Finish capturing and tap **Process project** (choose *High* mesh resolution for small parts). The model
   downloads to the phone when processing finishes.
2. Open the project (tap its card, or the three-dot button → **View/Edit Scene**). Use the **crop** tool to
   cut away the table and background so only the object is left.
3. Tap **Export 3D model** and save the file to the phone. **Pick OBJ or GLB if the app offers a choice**
   (both work directly). If it only offers FBX or USDZ, see [file types](#10-supported-file-types) for a
   one-minute conversion.
4. *Alternative:* **Share to Sketchfab**, then on sketchfab.com open the model → **Download 3D Model** and pick
   glTF/GLB or the original format. The `.zip` you get can be used as-is.

**Desktop (RealityScan for PC, formerly RealityCapture):** align the photos, build the model, then export the
mesh with the **Export model** dialog as **OBJ**. If you defined markers / AprilTags as ground-control points,
the OBJ is already in your workspace frame and in metres. This is what `mesh_Model_7.obj` from our rig is.

> **Scale check.** Phone photogrammetry has no reliable real-world size unless markers fix it. Measure the
> longest side of the real object with calipers. If the cleaned STL is the wrong size, rerun with
> `--fit-largest-mm <that length>` (or set `fit_largest_mm` in `config.toml`).

---

## 4. Move the scan to the laptop

Any of these works. Put the file in a folder called `scans/` inside the repo (it is git-ignored, so your
scans stay private):

| From | How |
|---|---|
| iPhone → Mac | AirDrop the exported file |
| iPhone → Windows | iCloud Drive / Google Drive / OneDrive app, or a USB cable + the Windows Photos/Files app |
| Android → anything | Google Drive, or USB cable (choose *File transfer* on the phone) |
| Sketchfab | Download on the laptop directly; keep the `.zip` |
| RealityScan desktop | It's already on the PC; just copy the `.obj` (and its `.mtl`/textures, not required) |

Email works for small files but often blocks anything over ~25 MB. Scans are usually bigger.

Check the scan opens before going further:
```bash
python scan_to_print.py scans/my_object.obj --stl-only
```
This writes `out_print/my_object_<time>/my_object_clean.stl`. Open it in Bambu Studio (or Windows 3D Viewer)
and check its size and that it is right side up. Use `--up y` if it lies on its side, and `--fit-largest-mm` if
the size is wrong.

> **Just want a printable STL of the scan?** That command is all you need: `my_object_clean.stl` is watertight,
> in millimetres and Z-up, and prints in any slicer with whatever single infill you choose.

---

## 5. Pick the probe points

```bash
python pick_points.py scans/my_object.obj
```
A 3D window opens.
- **Shift + left-click** on the surface to drop a point (in most Open3D versions **Shift + right-click** undoes
  the last one; otherwise just delete extra rows from the CSV afterwards).
- Rotate with the left mouse button, pan with the right, zoom with the wheel.
- Close the window when done.

The points are saved as `scans/my_object_points.csv` (x, y, z in **mm**, in the scan's own axes) with an empty
`stiffness` column. Tips:
- One point = one region of the print. 6–40 points is typical. More points mean smaller regions.
- Keep points at least ~15 mm apart. At 10 % infill the lines are ~4 mm apart, so a region needs room for a
  few of them.
- The probe comes straight down, so pick points on **upward-facing** surfaces with nothing above them.
- Put points on each distinct material (soft pad, hard rim, ...), not exactly on the boundary between them.

Already have coordinates from the gantry or another tool? Skip the picker and write the CSV yourself (see
`examples/points_template.csv`). Any of these column names work: `x y z` / `x_mm` / `X (mm)`, and
`stiffness` / `k` / `k (N/mm)`.

---

## 6. Probe and fill in stiffness

Press the probe into the real object at each point and record the stiffness **k = force ÷ indentation (N/mm)**.
Open the CSV in Excel, Google Sheets or any text editor and type each value into the `stiffness` column:

```
id,x,y,z,stiffness
0,12.0,12.0,15.0,0.6
1,12.0,38.0,15.0,0.9
...
```
Save it as **CSV** (`.xlsx` also works if `pandas` and `openpyxl` are installed). The program stops and tells
you which points are missing a stiffness value.

> No probe data yet? `--demo-points 12` puts 12 random points with random stiffness on the scan so you can
> test the whole chain. Don't print the result as if it meant anything.

---

## 7. Build the print

```bash
python scan_to_print.py scans/my_object.obj --points scans/my_object_points.csv
```

What happens:
1. **Scan:** reads the file, makes Z point up, converts to mm, simplifies very dense scans, and repairs it into
   a closed (watertight) solid. Phone scans always have holes, e.g. the underside the camera never saw.
2. **Points:** reads your CSV and moves the points exactly like the mesh. It refuses to continue if any point
   is more than 3 mm off the surface, which usually means the wrong scan file or units.
3. **Regions:** converts stiffness to infill (see [mapping](#how-stiffness-becomes-infill)) and splits the part
   into one region per point: each bit of the part belongs to its nearest probe point (a Voronoi partition).
   If every point is on a top surface, a point owns the whole column of material beneath it.
4. **Outputs:** writes the files you asked for (default: clean STL, region STLs, Bambu project, Bambu sliced).
5. **Check:** reads the sliced G-code back and measures the infill in each region, so you know the slicer did
   what was asked before you spend hours printing.

Choose outputs with `--formats` (comma-separated) or in `config.toml`:

| Format name | File | Needs |
|---|---|---|
| `bambu_sliced` | `<name>_bambu.gcode.3mf`: **finished, ready to print** | Bambu Studio installed |
| `bambu_project` | `<name>_bambu.3mf`: slice it yourself, all settings pre-loaded | Bambu Studio installed |
| `prusa_project` | `<name>_prusaslicer.3mf`: PrusaSlicer project, infill per part | – |
| `stl` | `<name>_clean.stl`: one watertight STL | – |
| `region_stls` | `region_stls/`: one STL per region + `infill_map.csv` | – |

Examples:
```bash
# only the ready-to-print file
python scan_to_print.py scans/box.obj --points scans/box_points.csv --formats bambu_sliced

# a project to slice yourself in Bambu Studio, and the STLs
python scan_to_print.py scans/box.obj --points scans/box_points.csv --formats bambu_project,stl,region_stls

# print at half size, with 3 top layers for this run only
python scan_to_print.py scans/box.obj --points scans/box_points.csv --scale 0.5 --top-shells 3

# a Sketchfab download, straight from the zip
python scan_to_print.py scans/box_sketchfab.zip --points scans/box_points.csv
```

---

## 8. Print it

**Ready-to-print file (`_bambu.gcode.3mf`):** open it in Bambu Studio. It opens on the *Preview* tab, already
sliced. Click **Print plate** and send it to the printer, or copy the file to the printer's SD card.

**Project to slice yourself (`_bambu.3mf`):**
1. Bambu Studio → **File → Open Project** → the `.3mf`.
2. Printer, filament and plate come in from the file: *Bambu Lab P1S 0.4 nozzle*, *Generic TPU*, *Textured PEI
   Plate*. Make sure the plate type matches the plate actually on the printer.
3. In the object list, expand the object: there is one part per region, named like `R03_58pct`. Click a part's
   settings icon to see its **Sparse infill density**.
4. **Slice plate** → **Preview**. Set the legend to *Line Type* and show only *Sparse infill*, then drag the
   layer slider to mid-height to see the density change from region to region.
5. **Print plate.**

**TPU on the P1S:**
- Feed TPU from the **external spool holder** on the back, not the AMS (soft TPU jams in it). Turn off
  *Use AMS* in the print dialog.
- Dry the spool if it has been open a while. Wet TPU strings and bubbles.
- Clean the Textured PEI plate (dish soap + water). Don't use glue stick for TPU on PEI.

---

## 9. Changing settings

**Everything adjustable is in [`config.toml`](config.toml)**, with a comment on every line. Edit it, or copy it
(`cp config.toml my_settings.toml`) and run with `--config my_settings.toml`. Command-line flags override the
file for a single run (`python scan_to_print.py --help` lists them).

### Common changes

| I want to... | In `config.toml` | Or on the command line |
|---|---|---|
| Use a different printer | `[bambu] printer = "Bambu Lab X1 Carbon 0.4 nozzle"` (and a matching `process`) | `--printer "..." --process "..."` |
| Use Bambu's own TPU | `[bambu] filament = "Bambu TPU 95A @BBL X1C"` | `--filament "..."` |
| Use another plate | `[bambu] plate = "Cool Plate"` | `--plate "Cool Plate"` |
| Thicker / thinner skins | `[print] top_shell_layers`, `bottom_shell_layers` | `--top-shells 3 --bottom-shells 3` |
| More / fewer walls | `[print] wall_loops` | `--walls 3` |
| Different infill pattern | `[print] infill_pattern = "gyroid"` | `--infill-pattern gyroid` |
| Different top surface | `[print] top_surface_pattern = "concentric"` | `--top-pattern concentric` |
| Change the stiffness → infill range | `[mapping] k_min, k_max, infill_min, infill_max` | – |
| Use measured calibration data | `[mapping] mode = "table"` + `table_k`, `table_infill` | `--mapping table` |
| Print smaller / bigger | `[scan] scale = 0.5` or `fit_largest_mm = 120` | `--scale 0.5` / `--fit-largest-mm 120` |
| Finer / coarser region map | `[regions] voxel_mm = 1.0` | `--voxel-mm 1.0` |
| Change any other Bambu setting | `[bambu.process_overrides]` / `[bambu.filament_overrides]` | – |
| Save outputs somewhere else | `[output] folder` | `--out D:/prints` |

Exact printer / process / filament names for your Bambu Studio install:
```bash
python scan_to_print.py --list-bambu-profiles
```

**Any other Bambu setting.** Use its internal key (open any Bambu `.3mf` as a zip and look in
`Metadata/project_settings.config` for names) under `[bambu.process_overrides]`, e.g.
```toml
[bambu.process_overrides]
outer_wall_speed = "40"
brim_type = "no_brim"
```

> **Why shells need two settings.** Bambu has both `top_shell_layers` (a count) and `top_shell_thickness` (a
> minimum in mm). The stock profile's 0.8 mm minimum would force 4 top layers even if you ask for 2, so the
> code sets both thicknesses to 0 and the layer counts decide. Same for the bottom.

### How stiffness becomes infill
- `mode = "linear"` (default): `k_min` → `infill_min`, `k_max` → `infill_max`, straight line in between,
  clipped outside. The defaults (0.5–3.0 N/mm → 10–90 %) are a **placeholder** until you calibrate.
- `mode = "table"`: piecewise-linear through your measured pairs. Print the calibration block
  ([section 11](#11-calibration-block-measure-the-real-stiffness-vs-infill-curve)), probe each patch, and put
  the results in `table_k` (stiffness you measured) and `table_infill` (that patch's infill).

### Why rectilinear infill by default
Rectilinear lines sit `0.41 mm / density` apart (≈ 4 mm at 10 %), so even small, sparse regions get several
full lines. Gyroid and honeycomb need a whole repeating unit inside a region to behave as designed, which small
low-density regions don't have room for. Switch with `infill_pattern` if your regions are large.

---

## 10. Supported file types

**Scans**

| Type | Read? | Notes |
|---|---|---|
| `.obj` | ✅ | Best choice. Textures (`.mtl`, images) and vertex colours are ignored and not needed. Large files use a fast reader. |
| `.glb` / `.gltf` | ✅ | Y-up by the glTF standard; the code turns it Z-up automatically. |
| `.ply` | ✅ | Must be a **mesh**. A point-cloud PLY is rejected with a message. |
| `.stl` | ✅ | |
| `.off`, `.3mf` | ✅ | |
| `.dae` (Collada) | ✅ | Needs `pip install pycollada`. |
| `.zip` | ✅ | Sketchfab / RealityScan downloads: the largest readable mesh inside is used. |
| `.fbx`, `.usdz`, `.usd`, `.abc`, `.blend` | ❌ → convert | Open in **Blender** (free): File → Import → (that format), then File → Export → Wavefront (.obj). Or re-export from the app as OBJ/GLB. |

Units are detected (a model under 5 units across is taken to be in metres) and can be forced with
`--units m|cm|mm|in`. The up axis is detected per format and can be forced with `--up y|z`.

**Probe points**

| Type | Notes |
|---|---|
| `.csv`, `.tsv`, `.txt` | Any delimiter (comma, tab, semicolon, space). A header row is optional: without one, columns are x, y, z, stiffness. Lines starting with `#` are ignored. |
| `.json` | `[{"x":..,"y":..,"z":..,"stiffness":..}, ...]`, or `{"points": [...]}`, or `[[x, y, z, k], ...]` |
| `.xlsx` | Needs `pandas` + `openpyxl`. |

Point units are detected (mm from `pick_points.py`, or the scan's own units such as metres): the code tries
both and keeps whichever lands on the surface. Force it with `--points-units mm|m|cm|in|file`.

**Outputs**: `.gcode.3mf` (Bambu, sliced), `.3mf` (Bambu project), `.3mf` (PrusaSlicer project), `.stl`, plus
`.csv` / `.png` / `.glb` reports.

---

## 11. Calibration block (measure the real stiffness-vs-infill curve)

The linear stiffness → infill mapping is a guess until you measure how stiff printed TPU actually is at each
infill. `make_calibration_block.py` makes a 100 × 60 × 20 mm block split into 8 Voronoi patches from 10 % to
85 % infill. The largest patch is the sparsest, so the gradient is easy to see.

```bash
python make_calibration_block.py
```
It writes (in `calibration_block/run_<time>/`) the Bambu project, an *open-top* copy that shows the infill
on the printed top face (for display, not probing), a map PNG with each patch's probe point, and
`probe_points.csv`. Print it, probe each patch at its red cross, and put the results into `[mapping] table_k` /
`table_infill`.

---

## 12. How it works / where to change the code

| File | What it does | Change it to... |
|---|---|---|
| `scan_to_print.py` | **Main program**: config, steps 1–5, outputs | add an output format, change the order of steps |
| `config.toml` | all settings | change behaviour without touching code |
| `scan_io.py` | reading scans + points in every format, units/axes, putting points on the mesh | support a new file type, new column names (`COLUMN_ALIASES`) |
| `bambu_profiles.py` | finds Bambu Studio, merges its profiles, exports settings, slices via the CLI | new printer quirks, friendly pattern names (`PATTERN_NAMES`) |
| `scan_to_infill.py` | the geometry engine: mesh repair, voxels, Voronoi regions, 3MF writer, G-code reader and infill verification | region logic, verification maths |
| `make_calibration_block.py` | calibration block + the Bambu 3MF writer (`write_bambu_3mf`) | block size, patch count, Bambu 3MF layout |
| `pick_points.py` | Open3D point picker → points CSV | picker behaviour |
| `examples/` | test part + points in CSV and JSON, blank template | – |
| `profiles/` | PrusaSlicer 2.5 TPU profile (for `prusa_project` and the older script) | – |
| `run_pipeline.py`, `infill_poc/` | the first proof of concept (cube + PrusaSlicer) | kept for reference |

**The pipeline in one paragraph.** The scan is rebuilt from vertices and faces only, which merges the duplicate
vertices that texture seams create (a textured scan otherwise loads as thousands of loose pieces). It is turned
Z-up, scaled to mm and repaired watertight with trimesh, falling back to MeshFix. Probe points get exactly the
same transform. The part is voxelised and each voxel is labelled with its nearest probe point; with top-only
contacts the regions are exact 2-D Voronoi cells extruded vertically. The part is then cut into one solid per
region with manifold3d booleans. Separate parts slice in seconds, where slicer "modifier" volumes blew up at
40 regions (13.5 h and a crash). Each part carries its own infill in the `.3mf`. For Bambu, the full printer /
process / filament config is exported by Bambu Studio's own CLI and embedded in the project, which is then
sliced by the same CLI. Finally the G-code's sparse-infill moves are binned by region and turned into a
measured infill % to compare with the target.

**Adding an output format:** write a function in `scan_to_print.py` that takes `pieces`, `names`, `target`,
add its name to `FORMATS`, and call it in `main()` next to the others.

**Adding a slicer setting:** add it to `config.toml` under `[print]`, then map it to the Bambu key in
`bambu_process_overrides()` in `scan_to_print.py`. Or just put the Bambu key in `[bambu.process_overrides]`.

---

## 13. Troubleshooting

| Message / problem | Fix |
|---|---|
| `... is a point cloud (no triangles)` | Export the **mesh** from the scanning app, not the point cloud. |
| `FBX files can't be read directly` | Re-export as OBJ/GLB, or convert in Blender (see [file types](#10-supported-file-types)). |
| `Points [...] are ... mm off the mesh surface` | The points were picked on a different file, or the units are wrong: try `--points-units mm` or `m`, and pick on the same file you pass to `scan_to_print.py`. |
| `Points [...] have no stiffness value` | Fill in the `stiffness` column, or use `--stl-only` / `--demo-points`. |
| Print is the wrong size | `--fit-largest-mm <measured length>` or `--scale`. |
| Model lies on its side | `--up y` (or `--up z` if it was detected as y). |
| `Bambu Studio is not installed...` | Install it, set `[bambu] executable`, or drop the `bambu_*` formats. |
| `No Bambu process profile named ...` | Copy the exact name from `--list-bambu-profiles`. |
| `Bambu Studio CLI exited ...` | Open the `.log` named in the message (in `bambu_profiles/` of the run folder); `process not compatible with printer` means the process name doesn't match the printer. |
| Verification shows +2 to +4 % points everywhere | Infill anchor lines are back on. `config.toml` turns them off by default (`sparse_infill_anchor` and `sparse_infill_anchor_max` = `"0"` under `[bambu.process_overrides]`); check those lines are still there. |
| Very slow on a huge scan | Big RealityScan OBJs (hundreds of MB, millions of triangles) are read by a fast built-in reader (a 456 MB, 6.5 M-triangle scan loads in ~30 s) and simplified to `[scan] max_faces`. If it is still slow, lower `max_faces` or raise `[regions] voxel_mm`. Keep the laptop awake for long runs. |
| `Reading Excel needs ...` | `pip install -U pandas openpyxl`, or save the sheet as CSV. |
