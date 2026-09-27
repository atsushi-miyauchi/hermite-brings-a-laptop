#!/usr/bin/env python3
"""Generate reusable certified Hermite coefficient bounds for FJ k-spokes rounding."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import parse_int_set
from hermite_coefficients.generator import compute_coefficients_acb, write_coefficient_json


def compact_tag(ks: list[int]) -> str:
    if not ks:
        return "none"
    parts: list[str] = []
    start = prev = ks[0]
    for k in ks[1:]:
        if k == prev + 1:
            prev = k
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = k
    parts.append(str(start) if start == prev else f"{start}-{prev}")
    return "_".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ks", required=True, help="spokes parameters, e.g. 2-32 or 2,4,7-10")
    ap.add_argument("--M", type=int, required=True, help="largest Hermite degree to certify")
    ap.add_argument("--R", type=Fraction, default=Fraction(8), help="finite integration radius (default: 8)")
    ap.add_argument("--prec", type=int, default=256, help="Arb/ACB precision in bits (default: 256)")
    ap.add_argument("--digits", type=int, default=40, help="decimal outward-rounding digits (default: 40)")
    ap.add_argument("--rel-tol", default=None)
    ap.add_argument("--abs-tol", default=None)
    ap.add_argument("--deg-limit", type=int, default=None)
    ap.add_argument("--eval-limit", type=int, default=None)
    ap.add_argument("--depth-limit", type=int, default=None)
    ap.add_argument("--details-csv", default=None, help="optional per-integral certificate CSV")
    ap.add_argument("--output", default=None, help="output JSON; default is res/coeffs_k<range>_M<M>.json")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if args.M < 0:
        raise ValueError("M must be nonnegative")
    ks = parse_int_set(args.ks, minimum=2)
    if not ks:
        raise ValueError("--ks must contain at least one k>=2")

    output = Path(args.output) if args.output else Path(__file__).resolve().parent / "res" / f"coeffs_k{compact_tag(ks)}_M{args.M}.json"
    c_lower, c_upper = compute_coefficients_acb(
        M=args.M,
        ks=ks,
        R=args.R,
        prec=args.prec,
        digits=args.digits,
        rel_tol=args.rel_tol,
        abs_tol=args.abs_tol,
        deg_limit=args.deg_limit,
        eval_limit=args.eval_limit,
        depth_limit=args.depth_limit,
        progress=not args.quiet,
        details_csv=args.details_csv,
    )
    write_coefficient_json(
        output,
        M=args.M,
        ks=ks,
        R=args.R,
        prec=args.prec,
        digits=args.digits,
        rel_tol=args.rel_tol,
        abs_tol=args.abs_tol,
        deg_limit=args.deg_limit,
        eval_limit=args.eval_limit,
        depth_limit=args.depth_limit,
        c_lower=c_lower,
        c_upper=c_upper,
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
