"""Write the exact Python and library versions to documentation/ENVIRONMENT.txt.

Run this in the same environment where you run frozen_refit.py, then paste the
versions into requirements.txt.
"""
import platform
import sys
from importlib import metadata
from pathlib import Path

PKGS = ["numpy", "pandas", "scikit-learn", "scikit-survival"]
lines = [f"python=={platform.python_version()}", f"platform: {platform.platform()}"]
for p in PKGS:
    try:
        lines.append(f"{p}=={metadata.version(p)}")
    except metadata.PackageNotFoundError:
        lines.append(f"{p}: NOT INSTALLED")
out = Path(__file__).resolve().parent.parent / "documentation" / "ENVIRONMENT.txt"
out.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"written to {out}", file=sys.stderr)
