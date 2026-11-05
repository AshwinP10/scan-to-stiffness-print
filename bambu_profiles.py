"""Bambu Studio integration: find the install, resolve its system profiles, export settings, slice.

The installed profiles inherit from each other ("inherits": parent name). The CLI wants every key, so a
profile is resolved by merging child over parent. A Bambu project .3mf must also carry
Metadata/project_settings.config (the full merged plate config) or the CLI crashes on load; we get that
file from the CLI itself with --export-settings so it always matches the installed version.
"""
import json
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

import trimesh

# Friendly infill / surface names -> Bambu's internal values.
PATTERN_NAMES = {
    "rectilinear": "zig-zag", "zigzag": "zig-zag", "aligned_rectilinear": "alignedrectilinear",
    "monotonic_line": "monotonicline", "3d_honeycomb": "3dhoneycomb", "adaptive_cubic": "adaptivecubic",
    "tri_hexagon": "tri-hexagon", "hilbert": "hilbertcurve", "archimedean": "archimedeanchords",
    "octagram": "octagramspiral",
}

CANDIDATE_EXES = {
    "Windows": [r"C:\Program Files\Bambu Studio\bambu-studio.exe", r"C:\Program Files (x86)\Bambu Studio\bambu-studio.exe"],
    "Darwin": ["/Applications/BambuStudio.app/Contents/MacOS/BambuStudio"],
    "Linux": ["/usr/bin/bambu-studio", "/opt/bambu-studio/bambu-studio"],
}


def bambu_pattern(name):
    key = str(name).strip().lower().replace(" ", "_").replace("-", "_")
    return PATTERN_NAMES.get(key, str(name).strip().lower())


def find_executable(configured=""):
    if configured:
        p = Path(configured).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"Bambu Studio not found at {p} (set [bambu] executable in config.toml)")
    for c in CANDIDATE_EXES.get(platform.system(), []):
        if Path(c).exists():
            return Path(c)
    w = shutil.which("bambu-studio") or shutil.which("BambuStudio")
    if w:
        return Path(w)
    raise FileNotFoundError("Bambu Studio is not installed in a standard location. Install it from bambulab.com, "
                            "or set [bambu] executable in config.toml, or drop 'bambu_sliced' from [output] formats.")


def profile_dir(exe):
    exe = Path(exe)
    for p in (exe.parent / "resources" / "profiles" / "BBL",                     # Windows / Linux
              exe.parent.parent / "Resources" / "profiles" / "BBL"):              # macOS app bundle
        if p.is_dir():
            return p
    raise FileNotFoundError(f"No Bambu profiles folder next to {exe}")


def _index(pdir):
    idx = {}
    for kind in ("machine", "process", "filament"):
        for f in (pdir / kind).glob("*.json"):
            try:
                d = json.loads(f.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            if "name" in d:
                idx[(kind, d["name"])] = d
    return idx


def list_profiles(exe, printer_filter=""):
    """Printable profiles (no abstract bases), optionally only those compatible with one printer."""
    idx = _index(profile_dir(exe))
    out = {"machine": [], "process": [], "filament": []}
    for (kind, name), d in sorted(idx.items()):
        if d.get("instantiation") == "false":
            continue
        if printer_filter and kind != "machine" and printer_filter not in d.get("compatible_printers", [printer_filter]):
            continue
        out[kind].append(name)
    return out


def resolve(idx, kind, name):
    if (kind, name) not in idx:
        close = [n for (k, n) in idx if k == kind and name.split()[0].lower() in n.lower()][:8]
        raise KeyError(f"No Bambu {kind} profile named '{name}'. Similar: {close}  "
                       f"(python scan_to_print.py --list-bambu-profiles)")
    d = idx[(kind, name)]
    merged = resolve(idx, kind, d["inherits"]) if d.get("inherits") else {}
    merged.update(d)
    return merged


def _as_setting(kind, value):
    """Bambu stores settings as strings; filament settings are lists (one entry per filament slot)."""
    if isinstance(value, list):
        return [str(v) for v in value]
    if isinstance(value, bool):
        value = int(value)
    return [str(value)] if kind == "filament" else str(value)


def flat_profiles(exe, out_dir, printer, process, filament, plate, process_overrides=None, filament_overrides=None,
                  machine_overrides=None):
    """Write resolved machine.json, process.json, filament.json into out_dir; return their paths."""
    idx = _index(profile_dir(exe))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for kind, name, extra in [("machine", printer, {"curr_bed_type": plate, **(machine_overrides or {})}),
                              ("process", process, process_overrides or {}),
                              ("filament", filament, filament_overrides or {})]:
        d = resolve(idx, kind, name)
        d.update({k: _as_setting(kind, v) for k, v in extra.items()})
        for k in ("inherits", "instantiation"):
            d.pop(k, None)
        # The CLI checks process/filament compatibility against the printer by name, and only accepts a
        # "system" printer, so the machine profile keeps that tag even with the plate type added.
        d["name"], d["from"] = name, "User" if extra and kind != "machine" else "system"
        paths[kind] = out_dir / f"{kind}.json"
        paths[kind].write_text(json.dumps(d, indent=2), encoding="utf-8")
    return paths


def _run(exe, args, log_path, timeout=3600):
    proc = subprocess.run([str(exe), "--debug", "2", *map(str, args)], capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout)
    Path(log_path).write_text(f"$ bambu-studio {' '.join(map(str, args))}\n\n{proc.stdout}\n{proc.stderr}",
                              encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"Bambu Studio CLI exited {proc.returncode} (see {log_path})")


def version(exe):
    out = subprocess.run([str(exe), "--help"], capture_output=True, text=True, errors="replace").stdout
    return next((ln.split("-", 1)[1].rstrip(":") for ln in out.splitlines() if ln.startswith("BambuStudio-")), "unknown")


def project_settings(exe, profiles, log_path):
    """The full merged plate config (Metadata/project_settings.config), from the CLI's --export-settings."""
    with tempfile.TemporaryDirectory() as tmp:
        stub = Path(tmp) / "stub.stl"
        trimesh.creation.box(extents=(10, 10, 10)).export(stub)
        out = Path(tmp) / "settings.json"
        _run(exe, ["--load-settings", f"{profiles['machine']};{profiles['process']}", "--load-filaments",
                   profiles["filament"], "--export-settings", out, stub], log_path, timeout=300)
        return out.read_text(encoding="utf-8")


def slice_3mf(exe, project_3mf, profiles, out_dir, sliced_name, log_path):
    """Slice every plate of a Bambu project 3MF; returns the .gcode.3mf (G-code + project, ready to print)."""
    out_dir = Path(out_dir)
    _run(exe, ["--slice", "0", "--load-settings", f"{profiles['machine']};{profiles['process']}",
               "--load-filaments", profiles["filament"], "--outputdir", out_dir, "--export-3mf", sliced_name,
               project_3mf], log_path)
    sliced = out_dir / sliced_name
    if not sliced.exists():
        raise RuntimeError(f"Bambu Studio CLI finished but wrote no {sliced_name} (see {log_path})")
    return sliced
