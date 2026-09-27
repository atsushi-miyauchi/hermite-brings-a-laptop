#!/usr/bin/env python3
"""Run the exact Section 4 application certificate checks used in the manuscript."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable


def run(*args: str) -> None:
    cmd = [PYTHON, *args]
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> int:
    run(str(ROOT / "maxagree" / "verify.py"), "--quiet")

    maxagree_k = {
        3: "0.8151", 4: "0.8025", 5: "0.7964", 6: "0.7934",
        7: "0.7914", 8: "0.7902", 9: "0.7893", 10: "0.7885",
        11: "0.7878", 12: "0.7871", 13: "0.7866", 14: "0.7863",
        15: "0.7859", 16: "0.7856",
    }
    for K, gamma in maxagree_k.items():
        run(
            str(ROOT / "maxagree-k" / "verify.py"),
            "--K", str(K), "--target-ratio", gamma, "--quiet",
        )

    max_k_cut = {
        3: ("0.836008", "0.836009"),
        4: ("0.857487", "0.857488"),
        5: ("0.876609", "0.876610"),
        6: ("0.891534", "0.891535"),
        7: ("0.903259", "0.903260"),
        8: ("0.912664", "0.912665"),
        9: ("0.920366", "0.920367"),
        10: ("0.926786", "0.926787"),
        11: ("0.932222", "0.932223"),
        12: ("0.936885", "0.936886"),
        13: ("0.940931", "0.940932"),
        14: ("0.944476", "0.944477"),
        15: ("0.947609", "0.947610"),
        16: ("0.950398", "0.950399"),
    }
    for K, (lower, upper) in max_k_cut.items():
        run(
            str(ROOT / "max-k-cut" / "verify.py"),
            "--K", str(K), "--lower", lower, "--upper", upper, "--quiet",
        )

    run(str(ROOT / "modularity" / "verify.py"), "--quiet")
    run(str(ROOT / "modularity_limitation" / "verify.py"), "--quiet")
    print("All Section 4 application certificate checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
