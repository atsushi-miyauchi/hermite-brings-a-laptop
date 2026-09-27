#!/usr/bin/env python3
"""Exact utilities for the MaxAgree edge-wise upper-bound certificate."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def parse_frac(x) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x, 1)
    return Fraction(str(x).strip())


def frac_str(x: Fraction) -> str:
    x = Fraction(x)
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def f2float(x: Fraction) -> float:
    return x.numerator / x.denominator


def parse_int_set(spec: str) -> List[int]:
    """Parse e.g. '2,3,4-8' into [2,3,4,5,6,7,8]."""
    out = set()
    if not spec.strip():
        return []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a_s, b_s = part.split("-", 1)
            a, b = int(a_s), int(b_s)
            if a > b:
                raise ValueError(f"invalid range {part!r}")
            out.update(range(a, b + 1))
        else:
            out.add(int(part))
    if any(k < 1 for k in out):
        raise ValueError("rounding parameters k must be positive integers")
    return sorted(out)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_coeff_json(
    path: str,
    kfj: Sequence[int],
    M: int,
) -> Tuple[Dict[int, List[Fraction]], Dict[int, List[Fraction]]]:
    """Load c_lower/c_upper through degree M for the requested spokes k values.

    k=1 is the exact one-cluster rounding P_1(rho)=1 and therefore does not
    need coefficient data.
    """
    if M < 0:
        raise ValueError("M must be nonnegative")
    with open(path) as f:
        data = json.load(f)

    data_M = int(data.get("M", -1))
    if data_M >= 0 and M > data_M:
        raise ValueError(f"requested M={M}, but coefficient file only certifies through M={data_M}")

    c_lower_raw = data["c_lower"]
    c_upper_raw = data["c_upper"]
    c_lower: Dict[int, List[Fraction]] = {}
    c_upper: Dict[int, List[Fraction]] = {}

    for k in kfj:
        if k == 1:
            continue
        sk = str(k)
        if sk not in c_lower_raw or sk not in c_upper_raw:
            raise KeyError(f"coefficient file is missing k={k}")
        lo_all = [parse_frac(x) for x in c_lower_raw[sk]]
        hi_all = [parse_frac(x) for x in c_upper_raw[sk]]
        if len(lo_all) < M + 1 or len(hi_all) < M + 1:
            raise ValueError(f"k={k} has fewer than M+1={M+1} coefficient bounds")
        lo = lo_all[: M + 1]
        hi = hi_all[: M + 1]
        for m, (l, u) in enumerate(zip(lo, hi)):
            if l < 0:
                raise ValueError(f"negative c_lower at k={k}, m={m}")
            if l > u:
                raise ValueError(f"c_lower > c_upper at k={k}, m={m}")
        if sum(lo) > 1:
            raise ValueError(f"sum of c_lower exceeds 1 for k={k}")
        c_lower[k] = lo
        c_upper[k] = hi
    return c_lower, c_upper


def fj_lower_nonnegative(
    k: int,
    rho: Fraction,
    M: int,
    c_lower: Dict[int, List[Fraction]],
) -> Fraction:
    if not (0 <= rho <= 1):
        raise ValueError("MaxAgree FJ lower bound expects rho in [0,1]")
    if k == 1:
        return Fraction(1)
    out = Fraction(0)
    power = Fraction(1)
    for m in range(M + 1):
        out += c_lower[k][m] * power
        power *= rho
    return out


def fj_upper_nonnegative(
    k: int,
    rho: Fraction,
    M: int,
    c_lower: Dict[int, List[Fraction]],
    c_upper: Dict[int, List[Fraction]],
) -> Fraction:
    if not (0 <= rho <= 1):
        raise ValueError("MaxAgree FJ upper bound expects rho in [0,1]")
    if k == 1:
        return Fraction(1)
    out = Fraction(0)
    power = Fraction(1)
    for m in range(M + 1):
        out += c_upper[k][m] * power
        power *= rho
    tail_mass_upper = Fraction(1) - sum(c_lower[k])
    out += tail_mass_upper * power  # power = rho^(M+1)
    return out


def _exact_arb_endpoint_to_fraction(x) -> Fraction:
    """Convert an exact finite arb endpoint to an exact Python Fraction."""
    man, exp = x.man_exp()
    man_i, exp_i = int(man), int(exp)
    if exp_i >= 0:
        return Fraction(man_i << exp_i, 1)
    return Fraction(man_i, 1 << (-exp_i))


def gw_bounds_arb(k: int, rho: Fraction, prec: int = 256) -> Tuple[Fraction, Fraction]:
    """Certified rational enclosure of (1-acos(rho)/pi)^k using python-flint/Arb."""
    if k < 1:
        raise ValueError("GW k must be positive")
    if not (-1 <= rho <= 1):
        raise ValueError("rho must lie in [-1,1]")
    try:
        from flint import arb, ctx
    except Exception as exc:  # pragma: no cover - depends on external package
        raise RuntimeError(
            "GW bounds require python-flint. Install python-flint, or generate a certificate with --kgw ''."
        ) from exc

    old_prec = ctx.prec
    try:
        ctx.prec = int(prec)
        r = arb(f"{rho.numerator}/{rho.denominator}")
        val = (1 - r.acos() / arb.pi()) ** k
        lo = _exact_arb_endpoint_to_fraction(val.lower())
        hi = _exact_arb_endpoint_to_fraction(val.upper())
    finally:
        ctx.prec = old_prec

    lo = max(Fraction(0), lo)
    hi = min(Fraction(1), hi)
    if lo > hi:
        raise RuntimeError(f"invalid Arb enclosure for GW(k={k}, rho={rho}): [{lo},{hi}]")
    return lo, hi


def rounding_name(kind: str, k: int) -> str:
    return f"{kind}{k}"
