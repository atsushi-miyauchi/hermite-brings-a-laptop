#!/usr/bin/env python3
"""Check a MaxAgree edge-wise upper-bound certificate.

No LP solver and no floating-point arithmetic are used for the dual-feasibility
check or for U_cert.  The checker automatically resolves the reusable Hermite coefficient JSON
recorded by the certificate, unless --coeff-json is supplied.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from common import (
    fj_lower_nonnegative,
    fj_upper_nonnegative,
    frac_str,
    gw_bounds_arb,
    load_coeff_json,
    parse_frac,
    sha256_file,
)


def fail(msg: str):
    raise AssertionError(msg)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate_json")
    ap.add_argument("--coeff-json", default=None, help="override automatically located Hermite coefficient JSON")
    ap.add_argument("--require-below", default=None, help="optionally require U_cert < this rational/decimal value")
    ap.add_argument("--verify-gw-with-arb", action="store_true", help="independently audit stored GW intervals with Arb")
    ap.add_argument("--gw-verify-prec", type=int, default=512, help="precision in bits for --verify-gw-with-arb")
    args = ap.parse_args(argv)

    cert_path = Path(args.certificate_json).resolve()
    with cert_path.open() as f:
        cert = json.load(f)

    if cert.get("schema") != "edgewise-upper-bound-certificate-v1":
        fail("unexpected certificate schema")
    if cert.get("problem") != "MaxAgree":
        fail("certificate is not for MaxAgree")
    if parse_frac(cert["a"]) != 0:
        fail("MaxAgree certificate must have a=0")

    N = int(cert["N"])
    M = int(cert["M"])
    kfj = [int(k) for k in cert["kfj"]]
    kgw = [int(k) for k in cert["kgw"]]
    if N <= 0 or M < 0:
        fail("invalid N or M")

    source = cert.get("coefficient_source", {})
    if args.coeff_json is None:
        rel = source.get("relative_path")
        basename = source.get("basename")
        if rel:
            coeff_path = (cert_path.parent / rel).resolve()
        elif basename:
            coeff_path = cert_path.parent / basename
        else:
            fail("certificate does not specify its coefficient JSON")
    else:
        coeff_path = Path(args.coeff_json).resolve()
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
        if rho != Fraction(idx, N):
            fail(f"positive rho does not match grid index {idx}")
        if value < 0:
            fail("negative u variable")
        if idx in seen_u:
            fail(f"duplicate positive grid index {idx}")
        seen_u.add(idx)
        active_u.append((idx, rho, value))

    active_v = []
    seen_v = set()
    for row in cert["active_v"]:
        idx = int(row["grid_index"])
        rho = parse_frac(row["rho"])
        value = parse_frac(row["value"])
        if not (1 <= idx <= N):
            fail(f"negative grid index out of range: {idx}")
        if rho != Fraction(idx - 1, N):
            fail(f"negative rho does not match grid index {idx}")
        if value < 0:
            fail("negative v variable")
        if idx in seen_v:
            fail(f"duplicate negative grid index {idx}")
        seen_v.add(idx)
        active_v.append((idx, rho, value))

    normalization = sum((u * rho for _, rho, u in active_u), Fraction(0))
    normalization += sum((v * (1 - rho) for _, rho, v in active_v), Fraction(0))
    if normalization != 1:
        fail(f"dual normalization is {normalization}, not 1")

    rounding_rows = cert["roundings"]
    expected_names = [f"gw{k}" for k in kgw] + [f"fj{k}" for k in kfj]
    actual_names = [row["name"] for row in rounding_rows]
    if actual_names != expected_names:
        fail(f"rounding list mismatch: expected {expected_names}, got {actual_names}")

    bounds = cert["probability_bounds"]
    recomputed_values = {}
    for rr in rounding_rows:
        kind, k, name = rr["kind"], int(rr["k"]), rr["name"]
        if name not in bounds:
            fail(f"missing probability bounds for {name}")
        plus_upper = [parse_frac(x) for x in bounds[name]["plus_upper"]]
        minus_lower = [parse_frac(x) for x in bounds[name]["minus_lower"]]
        if len(plus_upper) != len(active_u) or len(minus_lower) != len(active_v):
            fail(f"probability-bound length mismatch for {name}")

        if kind == "fj":
            for j, (_, rho, _) in enumerate(active_u):
                expected = fj_upper_nonnegative(k, rho, M, c_lower, c_upper)
                if plus_upper[j] != expected:
                    fail(f"FJ upper bound mismatch for {name} at rho={rho}")
            for j, (_, rho, _) in enumerate(active_v):
                expected = fj_lower_nonnegative(k, rho, M, c_lower)
                if minus_lower[j] != expected:
                    fail(f"FJ lower bound mismatch for {name} at rho={rho}")
        elif kind == "gw":
            if args.verify_gw_with_arb:
                prec = max(int(cert.get("gw_prec_bits", 0)), args.gw_verify_prec)
                for j, (_, rho, _) in enumerate(active_u):
                    _, fresh_hi = gw_bounds_arb(k, rho, prec)
                    if plus_upper[j] < fresh_hi:
                        fail(f"stored GW upper bound fails Arb audit for {name} at rho={rho}")
                for j, (_, rho, _) in enumerate(active_v):
                    fresh_lo, _ = gw_bounds_arb(k, rho, prec)
                    if minus_lower[j] > fresh_lo:
                        fail(f"stored GW lower bound fails Arb audit for {name} at rho={rho}")
        else:
            fail(f"unknown rounding kind {kind!r}")

        val = Fraction(0)
        for pu, (_, _, u) in zip(plus_upper, active_u):
            val += pu * u
        for pl, (_, _, v) in zip(minus_lower, active_v):
            val += (1 - pl) * v
        recomputed_values[name] = val

        claimed = parse_frac(cert["dual_values"][name])
        if val != claimed:
            fail(f"dual value mismatch for {name}: recomputed {val}, claimed {claimed}")

    U = max(recomputed_values.values())
    claimed_U = parse_frac(cert["U_cert"])
    if U != claimed_U:
        fail(f"U_cert mismatch: recomputed {U}, claimed {claimed_U}")

    if args.require_below is not None:
        threshold = parse_frac(args.require_below)
        if not U < threshold:
            fail(f"required U_cert < {threshold}, but U_cert = {U}")

    print("certificate OK")
    print(f"coefficient JSON = {coeff_path}")
    print(f"normalization = {frac_str(normalization)}")
    print(f"U_cert = {frac_str(U)}")
    print(f"U_cert ~= {float(U):.17g}")
    if kgw and not args.verify_gw_with_arb:
        print("note: stored GW rational enclosures were used as certified inputs; pass --verify-gw-with-arb to audit them")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
