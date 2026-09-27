#!/usr/bin/env python3
"""Check MaxAgree[K] edge-wise upper-bound certificates in exact arithmetic.

The checker uses no NumPy, SciPy, LP solver, or floating-point arithmetic for
verification.  It can check either one certificate JSON or, with --K-max,
all certificates for K=3,...,K_max produced by one generation run.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from common import fj_lower, fj_upper, frac_str, load_coeff_json, parse_frac, s_minus, s_plus, sha256_file


def fail(msg: str):
    raise AssertionError(msg)


def check_one(cert_path: Path, *, coeff_json: str | None, require_below: str | None, verbose: bool = True):
    cert_path = cert_path.resolve()
    with cert_path.open() as f:
        cert = json.load(f)

    if cert.get("schema") != "edgewise-upper-bound-certificate-v1":
        fail("unexpected certificate schema")
    if cert.get("problem") != "MaxAgree[K]":
        fail("certificate is not for MaxAgree[K]")

    K = int(cert["K"])
    if K < 3:
        fail("K must be at least 3")
    expected_a = Fraction(-1, K - 1)
    a = parse_frac(cert["a"])
    if a != expected_a:
        fail(f"a mismatch: expected {expected_a}, got {a}")
    if cert.get("negative_edge_mode") != "spokes_endpoint":
        fail("MaxAgree[K] certificate must use the spokes endpoint negative-edge reduction")

    N = int(cert["N"])
    M = int(cert["M"])
    if N <= 0 or M < 0:
        fail("invalid N or M")

    expected_kfj = list(range(1, K + 1))
    kfj = [int(k) for k in cert["kfj"]]
    if kfj != expected_kfj:
        fail(f"rounding family must be spokes k=1,...,K; got {kfj}")

    source = cert.get("coefficient_source", {})
    if coeff_json is None:
        rel = source.get("relative_path")
        basename = source.get("basename")
        if rel:
            coeff_path = (cert_path.parent / rel).resolve()
        elif basename:
            coeff_path = cert_path.parent / basename
        else:
            fail("certificate does not specify its coefficient JSON")
    else:
        coeff_path = Path(coeff_json).resolve()
    if not coeff_path.exists():
        fail(f"coefficient JSON not found: {coeff_path}")

    expected_sha = source.get("sha256")
    if expected_sha and sha256_file(str(coeff_path)) != expected_sha:
        fail("coefficient JSON SHA-256 does not match the certificate")
    c_lower, c_upper = load_coeff_json(str(coeff_path), kfj, M)

    active_u = []
    seen_u = set()
    for row in cert["active_u"]:
        idx = int(row["grid_index"])
        rho = parse_frac(row["rho"])
        value = parse_frac(row["value"])
        if not (1 <= idx <= N):
            fail(f"positive grid index out of range: {idx}")
        expected_rho = a + Fraction(idx, N) * (1 - a)
        if rho != expected_rho:
            fail(f"positive rho does not match rational grid at index {idx}")
        if value < 0:
            fail("negative u variable")
        if idx in seen_u:
            fail(f"duplicate positive grid index {idx}")
        seen_u.add(idx)
        active_u.append((idx, rho, value))

    active_v = []
    if len(cert["active_v"]) > 1:
        fail("spokes endpoint reduction has at most one active negative dual variable")
    for row in cert["active_v"]:
        rho = parse_frac(row["rho"])
        value = parse_frac(row["value"])
        if rho != a:
            fail(f"negative test point must be rho=a={a}")
        if value < 0:
            fail("negative v variable")
        active_v.append((rho, value))

    normalization = sum((u * s_plus(a, rho) for _, rho, u in active_u), Fraction(0))
    normalization += sum((v * s_minus(a, rho) for rho, v in active_v), Fraction(0))
    if normalization != 1:
        fail(f"dual normalization is {normalization}, not 1")

    rounding_rows = cert["roundings"]
    expected_names = [f"fj{k}" for k in expected_kfj]
    actual_names = [row["name"] for row in rounding_rows]
    if actual_names != expected_names:
        fail(f"rounding list mismatch: expected {expected_names}, got {actual_names}")

    bounds = cert["probability_bounds"]
    recomputed_values = {}
    for k in expected_kfj:
        name = f"fj{k}"
        if name not in bounds:
            fail(f"missing probability bounds for {name}")
        plus_upper = [parse_frac(x) for x in bounds[name]["plus_upper"]]
        minus_lower = [parse_frac(x) for x in bounds[name]["minus_lower"]]
        if len(plus_upper) != len(active_u) or len(minus_lower) != len(active_v):
            fail(f"probability-bound length mismatch for {name}")

        for j, (_, rho, _) in enumerate(active_u):
            expected = fj_upper(k, rho, M, c_lower, c_upper)
            if plus_upper[j] != expected:
                fail(f"FJ upper bound mismatch for {name} at rho={rho}")
        for j, (rho, _) in enumerate(active_v):
            expected = fj_lower(k, rho, M, c_lower, c_upper)
            if minus_lower[j] != expected:
                fail(f"FJ lower bound mismatch for {name} at rho={rho}")

        val = Fraction(0)
        for pu, (_, _, u) in zip(plus_upper, active_u):
            val += pu * u
        for pl, (_, v) in zip(minus_lower, active_v):
            val += (1 - pl) * v
        recomputed_values[name] = val

        claimed = parse_frac(cert["dual_values"][name])
        if val != claimed:
            fail(f"dual value mismatch for {name}: recomputed {val}, claimed {claimed}")

    U = max(recomputed_values.values())
    claimed_U = parse_frac(cert["U_cert"])
    if U != claimed_U:
        fail(f"U_cert mismatch: recomputed {U}, claimed {claimed_U}")

    if require_below is not None:
        threshold = parse_frac(require_below)
        if not U < threshold:
            fail(f"required U_cert < {threshold}, but U_cert = {U}")

    if verbose:
        print("certificate OK")
        print(f"K = {K}")
        print(f"a = {frac_str(a)}")
        print(f"coefficient JSON = {coeff_path}")
        print(f"normalization = {frac_str(normalization)}")
        print(f"U_cert = {frac_str(U)}")
        print(f"U_cert ~= {float(U):.17g}")

    return K, U, coeff_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate_json", nargs="?", help="one certificate JSON to check")
    ap.add_argument("--K-max", dest="K_max", type=int, default=None, help="check K=3,...,K_max in ./res/")
    ap.add_argument("--N", "--grid", dest="N", type=int, default=100000, help="grid denominator for --K-max mode")
    ap.add_argument("--M", type=int, default=32, help="truncation degree for --K-max mode")
    ap.add_argument("--coeff-json", default=None, help="override automatically located Hermite coefficient JSON")
    ap.add_argument(
        "--require-below",
        default=None,
        help="optionally require U_cert < this value (single-certificate mode only)",
    )
    args = ap.parse_args(argv)

    if args.certificate_json is not None and args.K_max is not None:
        ap.error("specify either one certificate_json or --K-max, not both")
    if args.certificate_json is None and args.K_max is None:
        ap.error("specify certificate_json or --K-max")

    if args.certificate_json is not None:
        check_one(
            Path(args.certificate_json),
            coeff_json=args.coeff_json,
            require_below=args.require_below,
            verbose=True,
        )
        return 0

    if args.require_below is not None:
        ap.error("--require-below is only supported in single-certificate mode")
    if args.K_max < 3:
        ap.error("--K-max must be at least 3")
    if args.N <= 0 or args.M < 0:
        ap.error("N must be positive and M must be nonnegative")

    script_dir = Path(__file__).resolve().parent
    res_dir = script_dir / "res"
    rows = []
    common_coeff = None
    for K in range(3, args.K_max + 1):
        cert_path = res_dir / f"maxagreek_K{K}_ub_N{args.N}_M{args.M}.json"
        if not cert_path.exists():
            fail(f"certificate not found for K={K}: {cert_path}")
        k, U, coeff_path = check_one(
            cert_path,
            coeff_json=args.coeff_json,
            require_below=None,
            verbose=False,
        )
        if common_coeff is None:
            common_coeff = coeff_path
        elif coeff_path != common_coeff:
            fail("batch certificates do not reference the same coefficient JSON")
        rows.append((k, U))

    print(f"all certificates OK for K=3,...,{args.K_max}")
    if common_coeff is not None:
        print(f"coefficient JSON = {common_coeff}")
    for K, U in rows:
        print(f"  K={K}: U_cert = {frac_str(U)} ~= {float(U):.17g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
