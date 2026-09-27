#!/usr/bin/env python3
"""Exact utilities for MaxAgree[K] edge-wise upper-bound certificates."""

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


def s_plus(a: Fraction, rho: Fraction) -> Fraction:
    return (rho - a) / (1 - a)


def s_minus(a: Fraction, rho: Fraction) -> Fraction:
    return (1 - rho) / (1 - a)


def fj_lower(
    k: int,
    rho: Fraction,
    M: int,
    c_lower: Dict[int, List[Fraction]],
    c_upper: Dict[int, List[Fraction]],
) -> Fraction:
    """Proposition 2.4 lower bound, including rho<0."""
    if not (-1 <= rho <= 1):
        raise ValueError("rho must lie in [-1,1]")
    if k == 1:
        return Fraction(1)
    if rho >= 0:
        out = Fraction(0)
        power = Fraction(1)
        for m in range(M + 1):
            out += c_lower[k][m] * power
            power *= rho
        return out

    out = Fraction(0)
    power = Fraction(1)
    for m in range(M + 1):
        coeff = c_lower[k][m] if (m % 2 == 0) else c_upper[k][m]
        out += coeff * power
        power *= rho
    tail_mass_upper = Fraction(1) - sum(c_lower[k])
    out -= (abs(rho) ** (M + 1)) * tail_mass_upper
    return out


def fj_upper(
    k: int,
    rho: Fraction,
    M: int,
    c_lower: Dict[int, List[Fraction]],
    c_upper: Dict[int, List[Fraction]],
) -> Fraction:
    """Proposition 2.4 upper bound, including rho<0."""
    if not (-1 <= rho <= 1):
        raise ValueError("rho must lie in [-1,1]")
    if k == 1:
        return Fraction(1)
    if rho >= 0:
        out = Fraction(0)
        power = Fraction(1)
        for m in range(M + 1):
            out += c_upper[k][m] * power
            power *= rho
        tail_mass_upper = Fraction(1) - sum(c_lower[k])
        out += (rho ** (M + 1)) * tail_mass_upper
        return out

    out = Fraction(0)
    power = Fraction(1)
    for m in range(M + 1):
        coeff = c_upper[k][m] if (m % 2 == 0) else c_lower[k][m]
        out += coeff * power
        power *= rho
    tail_mass_upper = Fraction(1) - sum(c_lower[k])
    out += (abs(rho) ** (M + 1)) * tail_mass_upper
    return out
