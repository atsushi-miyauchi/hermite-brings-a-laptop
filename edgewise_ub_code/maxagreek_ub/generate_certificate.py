#!/usr/bin/env python3
"""Generate MaxAgree[K] edge-wise upper-bound certificates for K=3,...,K_max.

For each K>=3, a_K=-1/(K-1), and the rounding family is spokes k=1,...,K.
Given K_max, the generator loads a reusable Hermite coefficient table and reuses it
to solve the finite dual LP separately for every K=3,...,K_max.  Corollary 2.2
reduces the negative-edge side to the single endpoint rho=a_K.

All generated files are written to ./res/ next to this script.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from certlib.exact import load_coefficients

from common import f2float, fj_lower, fj_upper, frac_str, parse_frac, s_minus, s_plus, sha256_file


def rationalize_nonnegative(x: float, max_den: int, zero_tol: float) -> Fraction:
    if x <= zero_tol:
        return Fraction(0)
    q = Fraction(str(float(x))).limit_denominator(max_den)
    if q < 0:
        raise ValueError(f"negative rationalized dual variable: {q}")
    return q


def horner_float(coeff: np.ndarray, rho: np.ndarray) -> np.ndarray:
    out = np.zeros_like(rho)
    for c in coeff[::-1]:
        out = out * rho + c
    return out


def fj_upper_float(k, rho, M, c_lower, c_upper):
    if k == 1:
        return np.ones_like(rho)
    lo = np.array([f2float(x) for x in c_lower[k]], dtype=np.float64)
    hi = np.array([f2float(x) for x in c_upper[k]], dtype=np.float64)
    out = np.empty_like(rho)
    nonneg = rho >= 0.0
    if np.any(nonneg):
        rr = rho[nonneg]
        out[nonneg] = horner_float(hi, rr)
        out[nonneg] += (1.0 - float(sum(c_lower[k]))) * np.power(rr, M + 1)
    if np.any(~nonneg):
        rr = rho[~nonneg]
        mixed = np.array([hi[m] if m % 2 == 0 else lo[m] for m in range(M + 1)], dtype=np.float64)
        out[~nonneg] = horner_float(mixed, rr)
        out[~nonneg] += (1.0 - float(sum(c_lower[k]))) * np.power(np.abs(rr), M + 1)
    return out


def fj_lower_float_scalar(k, rho, M, c_lower, c_upper):
    if k == 1:
        return 1.0
    lo = np.array([f2float(x) for x in c_lower[k]], dtype=np.float64)
    hi = np.array([f2float(x) for x in c_upper[k]], dtype=np.float64)
    if rho >= 0:
        val = 0.0
        for c in lo[::-1]:
            val = val * rho + c
        return val
    mixed = np.array([lo[m] if m % 2 == 0 else hi[m] for m in range(M + 1)], dtype=np.float64)
    val = 0.0
    for c in mixed[::-1]:
        val = val * rho + c
    val -= (1.0 - float(sum(c_lower[k]))) * abs(rho) ** (M + 1)
    return val


def build_and_solve_float_dual(K, N, M, c_lower, c_upper):
    a = -1.0 / float(K - 1)
    t = np.arange(1, N + 1, dtype=np.float64) / float(N)
    rho_plus = a + t * (1.0 - a)
    splus = t  # exactly s_a^+(rho_i^+)=i/N

    pplus = np.empty((K, N), dtype=np.float64)
    qminus = np.empty(K, dtype=np.float64)
    for row, k in enumerate(range(1, K + 1)):
        pplus[row, :] = fj_upper_float(k, rho_plus, M, c_lower, c_upper)
        qminus[row] = 1.0 - fj_lower_float_scalar(k, a, M, c_lower, c_upper)

    nvars = N + 2
    vidx = N
    Uidx = N + 1
    c = np.zeros(nvars, dtype=np.float64)
    c[Uidx] = 1.0

    eq = np.zeros(nvars, dtype=np.float64)
    eq[:N] = splus
    eq[vidx] = 1.0  # s_a^-(a)=1
    Aeq = sparse.csr_matrix(eq.reshape(1, -1))
    beq = np.array([1.0], dtype=np.float64)

    Aub = sparse.hstack(
        [
            sparse.csr_matrix(pplus),
            sparse.csr_matrix(qminus.reshape(-1, 1)),
            sparse.csr_matrix(-np.ones((K, 1), dtype=np.float64)),
        ],
        format="csr",
    )
    bub = np.zeros(K, dtype=np.float64)

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
    return res


def generate_certificate_for_K(*, K, args, c_lower, c_upper, coeff_path: Path, res_dir: Path):
    a = Fraction(-1, K - 1)
    kfj = list(range(1, K + 1))

    print(
        f"solving MaxAgree[{K}] finite dual: a={frac_str(a)}, N={args.N}, M={args.M}, spokes=1-{K}",
        file=sys.stderr,
        flush=True,
    )
    res = build_and_solve_float_dual(K, args.N, args.M, c_lower, c_upper)

    active_u_raw = []
    for j in range(args.N):
        q = rationalize_nonnegative(float(res.x[j]), args.max_den, args.zero_tol)
        if q:
            idx = j + 1
            rho = a + Fraction(idx, args.N) * (1 - a)
            active_u_raw.append((idx, rho, q))

    v_raw = rationalize_nonnegative(float(res.x[args.N]), args.max_den, args.zero_tol)
    Z = sum((q * s_plus(a, rho) for _, rho, q in active_u_raw), Fraction(0))
    Z += v_raw * s_minus(a, a)
    if Z <= 0:
        raise RuntimeError(f"K={K}: rationalized dual support has nonpositive normalization factor Z")

    active_u = [(idx, rho, raw, raw / Z) for idx, rho, raw in active_u_raw]
    active_v = [(a, v_raw, v_raw / Z)] if v_raw else []

    probability_bounds: Dict[str, Dict[str, List[str]]] = {}
    dual_values: Dict[str, Fraction] = {}
    for k in kfj:
        name = f"fj{k}"
        plus_upper = [fj_upper(k, rho, args.M, c_lower, c_upper) for _, rho, _, _ in active_u]
        minus_lower = [fj_lower(k, a, args.M, c_lower, c_upper)] if active_v else []

        val = Fraction(0)
        for pu, (_, _, _, u) in zip(plus_upper, active_u):
            val += pu * u
        if active_v:
            val += (1 - minus_lower[0]) * active_v[0][2]

        probability_bounds[name] = {
            "plus_upper": [frac_str(x) for x in plus_upper],
            "minus_lower": [frac_str(x) for x in minus_lower],
        }
        dual_values[name] = val

    U_cert = max(dual_values.values())
    cert = {
        "schema": "edgewise-upper-bound-certificate-v1",
        "problem": "MaxAgree[K]",
        "K": K,
        "K_max_generation_run": args.K_max,
        "a": frac_str(a),
        "N": args.N,
        "M": args.M,
        "kfj": kfj,
        "roundings": [{"kind": "fj", "k": k, "name": f"fj{k}"} for k in kfj],
        "negative_edge_mode": "spokes_endpoint",
        "coefficient_source": {
            "sha256": sha256_file(str(coeff_path)),
            "relative_path": os.path.relpath(coeff_path, res_dir),
            "basename": coeff_path.name,
        },
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
            {"rho": frac_str(rho), "raw_value": frac_str(raw), "value": frac_str(value)}
            for rho, raw, value in active_v
        ],
        "probability_bounds": probability_bounds,
        "dual_values": {name: frac_str(v) for name, v in dual_values.items()},
        "U_cert": frac_str(U_cert),
    }

    output_path = res_dir / f"maxagreek_K{K}_ub_N{args.N}_M{args.M}.json"
    with output_path.open("w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
        f.write("\n")

    return {
        "K": K,
        "float_objective": float(res.fun),
        "U_cert": U_cert,
        "active_u": len(active_u),
        "active_v": len(active_v),
        "certificate": output_path,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--K-max",
        dest="K_max",
        type=int,
        required=True,
        help="generate certificates for every K=3,...,K_max",
    )
    ap.add_argument("--N", "--grid", dest="N", type=int, default=100000, help="positive-edge grid denominator")
    ap.add_argument("--M", type=int, default=32, help="Hermite truncation degree")
    ap.add_argument("--coeff-json", default=None, help="reusable Hermite coefficient JSON")
    ap.add_argument("--max-den", type=int, default=10**12, help="max denominator when rationalizing dual variables")
    ap.add_argument("--zero-tol", type=float, default=1e-12, help="dual values at most this are set to zero")
    args = ap.parse_args(argv)

    if args.K_max < 3:
        raise ValueError("K_max must be at least 3")
    if args.N <= 0:
        raise ValueError("N must be positive")
    if args.M < 0:
        raise ValueError("M must be nonnegative")
    if args.max_den <= 0:
        raise ValueError("--max-den must be positive")

    K_max = args.K_max
    nontrivial_kfj = list(range(2, K_max + 1))

    script_dir = Path(__file__).resolve().parent
    res_dir = script_dir / "res"
    res_dir.mkdir(parents=True, exist_ok=True)
    coeff_path = Path(args.coeff_json).resolve() if args.coeff_json else (ROOT / "hermite_coefficients" / "res" / "coeffs_k2-16_M32.json")
    table = load_coefficients(coeff_path)
    table.require(nontrivial_kfj, args.M)
    c_lower = {k: table.c_lower[k][: args.M + 1] for k in nontrivial_kfj}
    c_upper = {k: table.c_upper[k][: args.M + 1] for k in nontrivial_kfj}

    results = []
    for K in range(3, K_max + 1):
        results.append(
            generate_certificate_for_K(
                K=K,
                args=args,
                c_lower=c_lower,
                c_upper=c_upper,
                coeff_path=coeff_path,
                res_dir=res_dir,
            )
        )

    print(f"coefficient bounds: {coeff_path}")
    print("\ncertified upper bounds")
    for row in results:
        print(
            f"  K={row['K']}: U_cert={frac_str(row['U_cert'])} "
            f"~= {f2float(row['U_cert']):.17g}; "
            f"active_u={row['active_u']}, active_v={row['active_v']}"
        )
        print(f"       certificate: {row['certificate']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
