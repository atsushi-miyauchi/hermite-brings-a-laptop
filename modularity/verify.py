#!/usr/bin/env python3
"""Exact verifier for the Section 4.6 modularity certificate (additive error 0.3790)."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import Q, frac_str, load_coefficients, verify_nonnegative_bernstein

K = 7
DEFAULT_EPSILON = Fraction(379, 1000)
DEFAULT_COEFFS = ROOT / "hermite_coefficients" / "res" / "coeffs_k2-16_M32.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coeff-json", default=str(DEFAULT_COEFFS))
    ap.add_argument("--M", type=int, default=32)
    ap.add_argument("--target-error", "--epsilon", dest="epsilon", type=Q, default=DEFAULT_EPSILON)
    ap.add_argument("--max-depth", type=int, default=32)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    epsilon = Fraction(args.epsilon)
    table = load_coefficients(args.coeff_json)
    M = table.require([K], args.M)

    # epsilon - x + P^-_{7,M}(x) - 1/7 >= 0 on [0,1].
    coeffs = [table.c_lower[K][m] for m in range(M + 1)]
    coeffs[0] += epsilon - Fraction(1, K)
    if len(coeffs) < 2:
        coeffs.append(Fraction(-1))
    else:
        coeffs[1] -= 1

    result = verify_nonnegative_bernstein(
        coeffs, lo=Fraction(0), hi=Fraction(1), max_depth=args.max_depth
    )
    ok = result.ok

    if not args.quiet:
        print("modularity certificate:", "PASS" if ok else "FAIL")
        print(f"  coefficient file = {table.path}")
        print(f"  rounding = {K}-spokes")
        print(f"  M = {M}")
        print(f"  target additive error = {frac_str(epsilon)}")
        print(f"  Bernstein intervals checked = {result.intervals_checked}")
        print(f"  max depth reached = {result.max_depth_reached}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
