#!/usr/bin/env python3
"""Exact verifier for the Section 4.3 MaxAgree certificate (ratio 0.7818)."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import (
    Q,
    frac_str,
    load_coefficients,
    verify_nonnegative_bernstein,
)

BETA = Fraction(3909, 5000)
MIXTURE = {4: Fraction(91, 250), 5: Fraction(159, 250)}
DEFAULT_COEFFS = ROOT / "hermite_coefficients" / "res" / "coeffs_k2-32_M24.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coeff-json", default=str(DEFAULT_COEFFS))
    ap.add_argument("--M", type=int, default=24)
    ap.add_argument("--target-ratio", type=Q, default=BETA)
    ap.add_argument("--max-depth", type=int, default=32)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    table = load_coefficients(args.coeff_json)
    M = table.require(MIXTURE.keys(), args.M)
    beta = Fraction(args.target_ratio)

    coeffs = [
        sum((w * table.c_lower[k][m] for k, w in MIXTURE.items()), Fraction(0))
        for m in range(M + 1)
    ]
    if len(coeffs) < 2:
        coeffs.extend([Fraction(0)] * (2 - len(coeffs)))
    coeffs[1] -= beta
    result = verify_nonnegative_bernstein(
        coeffs, lo=Fraction(0), hi=Fraction(1), max_depth=args.max_depth
    )

    # The negative-edge condition at a=0 is exact because P_k(0)=1/k.
    p0 = sum((w * Fraction(1, k) for k, w in MIXTURE.items()), Fraction(0))
    negative_slack = Fraction(1) - p0 - beta
    ok = result.ok and negative_slack >= 0

    if not args.quiet:
        print("MaxAgree certificate:", "PASS" if ok else "FAIL")
        print(f"  coefficient file = {table.path}")
        print(f"  M = {M}")
        print(f"  target ratio = {frac_str(beta)}")
        print("  mixture = " + ", ".join(f"lambda_{k}={frac_str(w)}" for k, w in MIXTURE.items()))
        print("  positive-edge Bernstein check:", "PASS" if result.ok else "FAIL")
        print(f"    intervals checked = {result.intervals_checked}")
        print(f"    max depth reached = {result.max_depth_reached}")
        print(f"  negative-edge endpoint slack = {frac_str(negative_slack)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
