#!/usr/bin/env python3
"""Exact verifier for the Section 4.4 MaxAgree[K] certificates in Table 2."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from certlib.exact import Q, decimal_str, frac_str, load_coefficients, verify_nonnegative_bernstein

DEFAULT_COEFFS = ROOT / "hermite_coefficients" / "res" / "coeffs_k2-16_M32.json"

TABLE2_MIXTURES: dict[int, dict[int, Fraction]] = {
    3: {2: Fraction(123, 1000), 3: Fraction(877, 1000)},
    4: {3: Fraction(720, 1000), 4: Fraction(280, 1000)},
    5: {3: Fraction(453, 1000), 4: Fraction(547, 1000)},
    6: {3: Fraction(291, 1000), 4: Fraction(709, 1000)},
    7: {3: Fraction(182, 1000), 4: Fraction(818, 1000)},
    8: {3: Fraction(100, 1000), 4: Fraction(900, 1000)},
    9: {3: Fraction(40, 1000), 4: Fraction(960, 1000)},
    10: {4: Fraction(986, 1000), 5: Fraction(14, 1000)},
    11: {4: Fraction(925, 1000), 5: Fraction(75, 1000)},
    12: {4: Fraction(878, 1000), 5: Fraction(122, 1000)},
    13: {4: Fraction(832, 1000), 5: Fraction(168, 1000)},
    14: {4: Fraction(800, 1000), 5: Fraction(200, 1000)},
    15: {4: Fraction(769, 1000), 5: Fraction(231, 1000)},
    16: {4: Fraction(742, 1000), 5: Fraction(258, 1000)},
}


def splus_coefficients(K: int, gamma: Fraction, degree: int) -> list[Fraction]:
    n = max(1, degree)
    out = [Fraction(0) for _ in range(n + 1)]
    out[0] = gamma / K
    out[1] = gamma * Fraction(K - 1, K)
    return out


def subtract_polynomials(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    n = max(len(p), len(q))
    out = [Fraction(0) for _ in range(n)]
    for i in range(n):
        out[i] = (p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def lower_positive_polynomial(table, mixture: dict[int, Fraction], M: int) -> list[Fraction]:
    return [
        sum((w * table.c_lower[k][m] for k, w in mixture.items()), Fraction(0))
        for m in range(M + 1)
    ]


def lower_negative_polynomial(table, mixture: dict[int, Fraction], M: int) -> list[Fraction]:
    coeffs = [Fraction(0) for _ in range(M + 2)]
    for k, w in mixture.items():
        for m in range(M + 1):
            c = table.c_lower[k][m] if m % 2 == 0 else table.c_upper[k][m]
            coeffs[m] += w * c
        residual = Fraction(1) - sum(table.c_lower[k][: M + 1], Fraction(0))
        coeffs[M + 1] += w * ((-1) ** M) * residual
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def upper_collision_at_left_endpoint(table, K: int, mixture: dict[int, Fraction], M: int) -> Fraction:
    a = Fraction(-1, K - 1)
    total = Fraction(0)
    for k, w in mixture.items():
        partial = Fraction(0)
        for m in range(M + 1):
            c = table.c_upper[k][m] if m % 2 == 0 else table.c_lower[k][m]
            partial += c * a**m
        residual = Fraction(1) - sum(table.c_lower[k][: M + 1], Fraction(0))
        total += w * (partial + abs(a) ** (M + 1) * residual)
    return total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--target-ratio", "--gamma", dest="gamma", type=Q, required=True)
    ap.add_argument("--coeff-json", default=str(DEFAULT_COEFFS))
    ap.add_argument("--M", type=int, default=32)
    ap.add_argument("--max-depth", type=int, default=40)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if args.K not in TABLE2_MIXTURES:
        raise ValueError("K must lie in 3,...,16")
    gamma = Fraction(args.gamma)
    mixture = TABLE2_MIXTURES[args.K]
    table = load_coefficients(args.coeff_json)
    M = table.require(mixture.keys(), args.M)
    a = Fraction(-1, args.K - 1)

    lower_pos = lower_positive_polynomial(table, mixture, M)
    q_pos = subtract_polynomials(lower_pos, splus_coefficients(args.K, gamma, len(lower_pos) - 1))
    pos = verify_nonnegative_bernstein(q_pos, lo=Fraction(0), hi=Fraction(1), max_depth=args.max_depth)

    lower_neg = lower_negative_polynomial(table, mixture, M)
    q_neg = subtract_polynomials(lower_neg, splus_coefficients(args.K, gamma, len(lower_neg) - 1))
    neg_pos = verify_nonnegative_bernstein(q_neg, lo=a, hi=Fraction(0), max_depth=args.max_depth)

    upper_Pa = upper_collision_at_left_endpoint(table, args.K, mixture, M)
    endpoint_slack = Fraction(1) - upper_Pa - gamma
    endpoint_ok = endpoint_slack >= 0
    ok = pos.ok and neg_pos.ok and endpoint_ok

    if not args.quiet:
        print("MaxAgree[K] certificate:", "PASS" if ok else "FAIL")
        print(f"  K = {args.K}")
        print(f"  target ratio = {frac_str(gamma)} ~ {decimal_str(gamma)}")
        print(f"  coefficient file = {table.path}")
        print(f"  M = {M}")
        print("  mixture = " + ", ".join(f"lambda_{k}={frac_str(w)}" for k, w in sorted(mixture.items())))
        print(f"  positive-edge [0,1] = {'PASS' if pos.ok else 'FAIL'}")
        print(f"  positive-edge [{frac_str(a)},0] = {'PASS' if neg_pos.ok else 'FAIL'}")
        print(f"  negative-edge endpoint = {'PASS' if endpoint_ok else 'FAIL'}")
        print(f"    upper P(a_K) = {frac_str(upper_Pa)}")
        print(f"    slack = {frac_str(endpoint_slack)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
