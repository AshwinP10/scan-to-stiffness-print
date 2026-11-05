import os
import re
import shutil
import subprocess
from pathlib import Path

PINNED_VERSION = "2.5.0"

_CANDIDATES = [
    Path(os.environ.get("LOCALAPPDATA", "")) / "Programs/Prusa3D/PrusaSlicer/prusa-slicer-console.exe",
    Path(os.environ.get("PROGRAMFILES", "C:/Program Files")) / "Prusa3D/PrusaSlicer/prusa-slicer-console.exe",
]


def find_prusaslicer(explicit=None):
    for c in [explicit, os.environ.get("PRUSASLICER_CONSOLE"), *_CANDIDATES,
              shutil.which("prusa-slicer-console"), shutil.which("prusa-slicer")]:
        if c and Path(c).is_file():
            return str(c)
    raise FileNotFoundError("PrusaSlicer console executable not found; pass --slicer or set PRUSASLICER_CONSOLE.")


def _run(exe, args, log_path):
    proc = subprocess.run([exe, "--loglevel", "2", *args], capture_output=True, text=True,
                          encoding="utf-8", errors="replace")
    Path(log_path).write_text(f"$ {exe} {' '.join(args)}\n\n[stdout]\n{proc.stdout}\n[stderr]\n{proc.stderr}",
                              encoding="utf-8")
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-15:]
        raise RuntimeError(f"PrusaSlicer exited {proc.returncode} (log: {log_path}):\n" + "\n".join(tail))
    return proc


def slicer_version(exe):
    proc = subprocess.run([exe, "--help"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    m = re.search(r"PrusaSlicer-(\d+\.\d+\.\d+)", proc.stdout + proc.stderr)
    if not m:
        raise RuntimeError("Could not determine PrusaSlicer version from --help output")
    return m.group(1)


def check_pinned_version(exe, allow_mismatch=False):
    version = slicer_version(exe)
    if version != PINNED_VERSION and not allow_mismatch:
        raise RuntimeError(f"PrusaSlicer {version} found, pipeline is pinned to {PINNED_VERSION}. "
                           "Install the pinned version or pass --allow-version-mismatch.")
    return version


def roundtrip_3mf(exe, model_3mf, out_3mf, log_path):
    """Have PrusaSlicer load and re-export the 3MF, proving it parsed every volume."""
    _run(exe, ["--export-3mf", "--dont-arrange", "--output", str(out_3mf), str(model_3mf)], log_path)


def slice_to_gcode(exe, model_3mf, profile_ini, out_gcode, log_path):
    _run(exe, ["--export-gcode", "--load", str(profile_ini), "--dont-arrange",
               "--output", str(out_gcode), str(model_3mf)], log_path)
