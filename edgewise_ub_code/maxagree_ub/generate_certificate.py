#!/usr/bin/env python3
"""Generate a MaxAgree edge-wise upper-bound certificate (Sections 6.1--6.2).

The generator loads reusable certified Hermite coefficient bounds and then solves
the finite dual LP in floating point only to locate a candidate support, rounds
that support to rational numbers, normalizes it exactly, and constructs U_cert
from certified rational collision-probability bounds.

All generated files are written to ./res/ next to this script.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from scipy import sparse
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from certlib.exact import load_coefficients

from common import (
    f2float,
    fj_lower_nonnegative,
    fj_upper_nonnegative,
    frac_str,
    gw_bounds_arb,
    parse_frac,
    parse_int_set,
    rounding_name,
    sha256_file,
)


def compact_int_set(values: Sequence[int]) -> str:
    """Return a compact deterministic range string, e.g. [1,2,3,5,6] -> '1-3_5-6'."""
    vals = sorted(set(int(x) for x in values))
    if not vals:
        return "none"
    parts = []
    start = prev = vals[0]
    for x in vals[1:]:
        if x == prev + 1:
            prev = x
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = x
    parts.append(str(start) if start == prev else f"{start}-{prev}")
    return "_".join(parts)


def rationalize_nonnegative(x: float, max_den: int, zero_tol: float) -> Fraction:
    if x <= zero_tol:
        return Fraction(0)
    q = Fraction(str(float(x))).limit_denominator(max_den)
    if q < 0:
        raise ValueError(f"negative rationalized dual variable: {q}")
    return q


def fj_poly_float(coeff: Sequence[Fraction], rho: np.ndarray) -> np.ndarray:
    a = np.array([f2float(x) for x in coeff], dtype=np.float64)
    out = np.zeros_like(rho)
    for c in a[::-1]:
        out = out * rho + c
    return out


def fj_upper_float(k, rho, M, c_lower, c_upper):
    if k == 1:
        return np.ones_like(rho)
    out = fj_poly_float(c_upper[k], rho)
    tail = f2float(Fraction(1) - sum(c_lower[k]))
    out += tail * np.power(rho, M + 1)
    return out


def fj_lower_float(k, rho, c_lower):
    if k == 1:
        return np.ones_like(rho)
    return fj_poly_float(c_lower[k], rho)


def gw_float(k: int, rho: np.ndarray) -> np.ndarray:
    return np.power(1.0 - np.arccos(np.clip(rho, 0.0, 1.0)) / np.pi, k)


def build_and_solve_float_dual(kfj, kgw, N, M, c_lower, c_upper):
    rho_plus = np.arange(1, N + 1, dtype=np.float64) / float(N)
    rho_minus = np.arange(0, N, dtype=np.float64) / float(N)
    roundings = [("gw", k) for k in kgw] + [("fj", k) for k in kfj]
    nr = len(roundings)

    pplus = np.empty((nr, N), dtype=np.float64)
    qminus = np.empty((nr, N), dtype=np.float64)
    for r, (kind, k) in enumerate(roundings):
        if kind == "fj":
            pplus[r, :] = fj_upper_float(k, rho_plus, M, c_lower, c_upper)
            qminus[r, :] = 1.0 - fj_lower_float(k, rho_minus, c_lower)
        else:
            pplus[r, :] = gw_float(k, rho_plus)
            qminus[r, :] = 1.0 - gw_float(k, rho_minus)

    nvars = 2 * N + 1
    Uidx = nvars - 1
    c = np.zeros(nvars, dtype=np.float64)
    c[Uidx] = 1.0

    eq = np.zeros(nvars, dtype=np.float64)
    eq[:N] = rho_plus
    eq[N : 2 * N] = 1.0 - rho_minus
    Aeq = sparse.csr_matrix(eq.reshape(1, -1))
    beq = np.array([1.0], dtype=np.float64)

    Aub = sparse.hstack(
        [
            sparse.csr_matrix(pplus),
            sparse.csr_matrix(qminus),
            sparse.csr_matrix(-np.ones((nr, 1), dtype=np.float64)),
        ],
        format="csr",
    )
    bub = np.zeros(nr, dtype=np.float64)

    res = linprog(
        c,
        A_ub=Aub,
        b_ub=bub,
        A_eq=Aeq,
        b_eq=beq,
        bounds=(0.0, None),
        method="highs",
    )
    if not res.success:
        raise RuntimeError(res.message)
    return roundings, res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kfj", default="1-32", help="spokes k values, e.g. '2,3,4-8' or '1-32'")
    ap.add_argument("--kgw", default="1-5", help="hyperplane k values, e.g. '1-5'; empty means none")
    ap.add_argument("--N", "--grid", dest="N", type=int, default=100000, help="grid denominator")
    ap.add_argument("--M", type=int, default=24, help="Hermite truncation degree")
    ap.add_argument("--coeff-json", default=None, help="reusable Hermite coefficient JSON")
    ap.add_argument("--max-den", type=int, default=10**12, help="max denominator when rationalizing dual variables")
    ap.add_argument("--zero-tol", type=float, default=1e-12, help="dual values at most this are set to zero")
    ap.add_argument("--gw-prec", type=int, default=256, help="Arb precision in bits for certified GW bounds")
    args = ap.parse_args(argv)

    if args.N <= 0:
        raise ValueError("N must be positive")
    if args.M < 0:
        raise ValueError("M must be nonnegative")
    if args.max_den <= 0:
        raise ValueError("--max-den must be positive")

    kfj = parse_int_set(args.kfj)
    kgw = parse_int_set(args.kgw)
    if not kfj and not kgw:
        raise ValueError("specify at least one rounding")

    script_dir = Path(__file__).resolve().parent
    res_dir = script_dir / "res"
    res_dir.mkdir(parents=True, exist_ok=True)
    fj_tag = compact_int_set(kfj)
    gw_tag = compact_int_set(kgw)
    coeff_path = Path(args.coeff_json).resolve() if args.coeff_json else (ROOT / "hermite_coefficients" / "res" / "coeffs_k2-32_M24.json")
    table = load_coefficients(coeff_path)
    nontrivial_kfj = [k for k in kfj if k >= 2]
    table.require(nontrivial_kfj, args.M)
    c_lower = {k: table.c_lower[k][: args.M + 1] for k in nontrivial_kfj}
    c_upper = {k: table.c_upper[k][: args.M + 1] for k in nontrivial_kfj}

    print(
        f"solving MaxAgree finite dual: N={args.N}, M={args.M}, FJ={kfj}, GW={kgw}",
        file=sys.stderr,
        flush=True,
    )
    roundings, res = build_and_solve_float_dual(kfj, kgw, args.N, args.M, c_lower, c_upper)

    active_u_raw = []
    active_v_raw = []
    for j in range(args.N):
        q = rationalize_nonnegative(float(res.x[j]), args.max_den, args.zero_tol)
        if q:
            active_u_raw.append((j + 1, Fraction(j + 1, args.N), q))
    for j in range(args.N):
        q = rationalize_nonnegative(float(res.x[args.N + j]), args.max_den, args.zero_tol)
        if q:
            active_v_raw.append((j + 1, Fraction(j, args.N), q))

    Z = sum((q * rho for _, rho, q in active_u_raw), Fraction(0))
    Z += sum((q * (1 - rho) for _, rho, q in active_v_raw), Fraction(0))
    if Z <= 0:
        raise RuntimeError("rationalized dual support has nonpositive normalization factor Z")

    active_u = [(idx, rho, q, q / Z) for idx, rho, q in active_u_raw]
    active_v = [(idx, rho, q, q / Z) for idx, rho, q in active_v_raw]

    bounds: Dict[str, Dict[str, List[str]]] = {}
    dual_values: Dict[str, Fraction] = {}
    for kind, k in roundings:
        name = rounding_name(kind, k)
        plus_upper: List[Fraction] = []
        minus_lower: List[Fraction] = []

        for _, rho, _, _ in active_u:
            if kind == "fj":
                pu = fj_upper_nonnegative(k, rho, args.M, c_lower, c_upper)
            else:
                _, pu = gw_bounds_arb(k, rho, args.gw_prec)
            plus_upper.append(pu)

        for _, rho, _, _ in active_v:
            if kind == "fj":
                pl = fj_lower_nonnegative(k, rho, args.M, c_lower)
            else:
                pl, _ = gw_bounds_arb(k, rho, args.gw_prec)
            minus_lower.append(pl)

        val = Fraction(0)
        for pu, (_, _, _, u) in zip(plus_upper, active_u):
            val += pu * u
        for pl, (_, _, _, v) in zip(minus_lower, active_v):
            val += (1 - pl) * v

        bounds[name] = {
            "plus_upper": [frac_str(x) for x in plus_upper],
            "minus_lower": [frac_str(x) for x in minus_lower],
        }
        dual_values[name] = val

    U_cert = max(dual_values.values())

    cert = {
        "schema": "edgewise-upper-bound-certificate-v1",
        "problem": "MaxAgree",
        "a": "0",
        "N": args.N,
        "M": args.M,
        "kfj": kfj,
        "kgw": kgw,
        "roundings": [{"kind": kind, "k": k, "name": rounding_name(kind, k)} for kind, k in roundings],
        "coefficient_source": {
            "sha256": sha256_file(str(coeff_path)),
            "relative_path": os.path.relpath(coeff_path, res_dir),
            "basename": coeff_path.name,
        },
        "gw_prec_bits": args.gw_prec,
        "candidate_float_lp_objective": float(res.fun),
        "rationalization": {
            "max_den": args.max_den,
            "zero_tol": args.zero_tol,
            "Z_before_normalization": frac_str(Z),
        },
        "active_u": [
            {"grid_index": idx, "rho": frac_str(rho), "raw_value": frac_str(raw), "value": frac_str(value)}
            for idx, rho, raw, value in active_u
        ],
        "active_v": [
            {"grid_index": idx, "rho": frac_str(rho), "raw_value": frac_str(raw), "value": frac_str(value)}
            for idx, rho, raw, value in active_v
        ],
        "probability_bounds": bounds,
        "dual_values": {name: frac_str(v) for name, v in dual_values.items()},
        "U_cert": frac_str(U_cert),
    }

    output_path = res_dir / f"maxagree_ub_N{args.N}_M{args.M}_fj{fj_tag}_gw{gw_tag}.json"
    with output_path.open("w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"coefficient bounds: {coeff_path}")
    print(f"floating dual objective: {res.fun:.17g}")
    print(f"active positive variables: {len(active_u)}")
    print(f"active negative variables: {len(active_v)}")
    print(f"U_cert = {frac_str(U_cert)}")
    print(f"U_cert ~= {f2float(U_cert):.17g}")
    print(f"certificate: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
