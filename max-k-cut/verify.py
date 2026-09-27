#!/usr/bin/env python3
"""Exact verifier for the Section 4.5 lower/upper bounds on the FJ Max K-Cut ratio."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import Q, collision_bounds_at, decimal_str, frac_str, load_coefficients

DEFAULT_COEFFS = ROOT / "hermite_coefficients" / "res" / "coeffs_k2-16_M32.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--lower", "--lower-gamma", dest="lower", type=Q, default=None,
                    help="candidate rigorous lower bound on alpha_K^FJ")
    ap.add_argument("--upper", "--upper-gamma", dest="upper", type=Q, default=None,
                    help="candidate rigorous upper bound on alpha_K^FJ")
    ap.add_argument("--coeff-json", default=str(DEFAULT_COEFFS))
    ap.add_argument("--M", type=int, default=32)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if args.K < 3:
        raise ValueError("K must be at least 3")
    if args.lower is None and args.upper is None:
        raise ValueError("provide --lower and/or --upper")

    table = load_coefficients(args.coeff_json)
    M = table.require([args.K], args.M)
    a = Fraction(-1, args.K - 1)
    p_lower, p_upper = collision_bounds_at(a, table.c_lower[args.K], table.c_upper[args.K], M)
    alpha_lower = Fraction(1) - p_upper
    alpha_upper = Fraction(1) - p_lower

    ok = True
    if args.lower is not None:
        ok = ok and alpha_lower >= Fraction(args.lower)
    if args.upper is not None:
        ok = ok and alpha_upper <= Fraction(args.upper)

    if not args.quiet:
        print("Max K-Cut certificate:", "PASS" if ok else "FAIL")
        print(f"  K = {args.K}")
        print(f"  coefficient file = {table.path}")
        print(f"  M = {M}")
        print(f"  a_K = {frac_str(a)}")
        print(f"  P_K(a_K) in [{frac_str(p_lower)}, {frac_str(p_upper)}]")
        print(
            f"  alpha_K^FJ in [{frac_str(alpha_lower)}, {frac_str(alpha_upper)}] "
            f"~= [{decimal_str(alpha_lower)}, {decimal_str(alpha_upper)}]"
        )
        if args.lower is not None:
            print(f"  requested lower bound {frac_str(args.lower)}: {'PASS' if alpha_lower >= args.lower else 'FAIL'}")
        if args.upper is not None:
            print(f"  requested upper bound {frac_str(args.upper)}: {'PASS' if alpha_upper <= args.upper else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
