#!/usr/bin/env python3
"""Exact verifier for the Section 4.6 modularity near-optimality witness."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import Q, collision_bounds_at, frac_str, load_coefficients

WITNESS_X = Fraction(4, 5)
DEFAULT_TARGET = Fraction(3785, 10000)
DEFAULT_SPOKES_MAX = 32
DEFAULT_COEFFS = ROOT / "hermite_coefficients" / "res" / "coeffs_k2-32_M24.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coeff-json", default=str(DEFAULT_COEFFS))
    ap.add_argument("--M", type=int, default=24)
    ap.add_argument("--target-value", "--target", dest="target", type=Q, default=DEFAULT_TARGET)
    ap.add_argument("--spokes-max", type=int, default=DEFAULT_SPOKES_MAX)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    target = Fraction(args.target)
    table = load_coefficients(args.coeff_json)
    required = list(range(2, args.spokes_max + 1))
    M = table.require(required, args.M) if required else args.M

    rows: list[tuple[int, Fraction, Fraction]] = []
    ok = True
    for k in range(1, args.spokes_max + 1):
        if k == 1:
            p_upper = Fraction(1)
        else:
            _, p_upper = collision_bounds_at(
                WITNESS_X, table.c_lower[k], table.c_upper[k], M
            )
        d_lower = WITNESS_X - p_upper + Fraction(1, k)
        slack = d_lower - target
        rows.append((k, d_lower, slack))
        ok = ok and slack >= 0

    if not args.quiet:
        print("modularity near-optimality witness:", "PASS" if ok else "FAIL")
        print(f"  coefficient file = {table.path}")
        print(f"  M = {M}")
        print(f"  witness x = {frac_str(WITNESS_X)}")
        print(f"  target = {frac_str(target)}")
        print(f"  spokes family = 1,...,{args.spokes_max}")
        for k, d, slack in rows:
            print(f"    k={k:2d}: D_lower={float(d):.12g} slack={float(slack):+.3e} {'PASS' if slack >= 0 else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
