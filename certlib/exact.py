from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence


def Q(x: str | int | float | Fraction) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x, 1)
    return Fraction(str(x).strip())


def frac_str(x: Fraction) -> str:
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def decimal_str(x: Fraction, digits: int = 15) -> str:
    x = Fraction(x)
    with localcontext() as ctx:
        ctx.prec = max(50, digits + 20)
        d = Decimal(x.numerator) / Decimal(x.denominator)
        return f"{d:.{digits}f}"


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_int_set(spec: str, *, minimum: int = 1) -> list[int]:
    out: set[int] = set()
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
    if any(k < minimum for k in out):
        raise ValueError(f"all integers must be at least {minimum}")
    return sorted(out)


@dataclass(frozen=True)
class CoefficientTable:
    path: Path
    stored_M: int
    c_lower: dict[int, list[Fraction]]
    c_upper: dict[int, list[Fraction]]
    metadata: dict

    def require(self, ks: Sequence[int], M: int | None = None) -> int:
        use_M = self.stored_M if M is None else int(M)
        if use_M < 0:
            raise ValueError("M must be nonnegative")
        if use_M > self.stored_M:
            raise ValueError(
                f"requested M={use_M}, but coefficient file only certifies through M={self.stored_M}"
            )
        for k in ks:
            if k == 1:
                continue
            if k not in self.c_lower or k not in self.c_upper:
                raise ValueError(f"coefficient file is missing k={k}")
            if len(self.c_lower[k]) <= use_M or len(self.c_upper[k]) <= use_M:
                raise ValueError(f"coefficient file does not contain k={k} through M={use_M}")
        return use_M


def _validate_one(k: int, lower: list[Fraction], upper: list[Fraction], M: int) -> None:
    if len(lower) < M + 1 or len(upper) < M + 1:
        raise ValueError(f"k={k}: coefficient arrays are shorter than M+1={M+1}")
    if k >= 2:
        expected = Fraction(1, k)
        if lower[0] != expected or upper[0] != expected:
            raise ValueError(
                f"k={k}: expected exact c_{{k,0}}=1/k={expected}, got [{lower[0]}, {upper[0]}]"
            )
    for m, (lo, hi) in enumerate(zip(lower[: M + 1], upper[: M + 1])):
        if lo < 0:
            raise ValueError(f"k={k}, m={m}: negative lower bound")
        if lo > hi:
            raise ValueError(f"k={k}, m={m}: lower bound exceeds upper bound")
    if sum(lower[: M + 1], Fraction(0)) > 1:
        raise ValueError(f"k={k}: sum of lower bounds through M exceeds 1")


def load_coefficients(path: str | Path) -> CoefficientTable:
    path = Path(path).resolve()
    with path.open() as f:
        data = json.load(f)

    # Preferred common schema: one uniform M and flat c_lower/c_upper dictionaries.
    if "c_lower" in data and "c_upper" in data:
        stored_M = int(data["M"])
        lower = {int(k): [Q(x) for x in vals] for k, vals in data["c_lower"].items()}
        upper = {int(k): [Q(x) for x in vals] for k, vals in data["c_upper"].items()}
    # Backward compatibility with earlier problem-specific files.
    elif "coefficients" in data:
        entries = data["coefficients"]
        if not entries:
            raise ValueError("coefficient JSON contains no entries")
        Ms = {int(entry["M"]) for entry in entries.values()}
        if len(Ms) != 1:
            raise ValueError("problem-specific coefficient JSON has nonuniform M values")
        stored_M = next(iter(Ms))
        lower = {int(k): [Q(x) for x in entry["c_lower"]] for k, entry in entries.items()}
        upper = {int(k): [Q(x) for x in entry["c_upper"]] for k, entry in entries.items()}
    else:
        raise ValueError("unrecognized coefficient JSON schema")

    if set(lower) != set(upper):
        raise ValueError("lower/upper coefficient files contain different k sets")
    for k in lower:
        _validate_one(k, lower[k], upper[k], stored_M)
    return CoefficientTable(path, stored_M, lower, upper, data)


def poly_eval(coeffs: Sequence[Fraction], x: Fraction) -> Fraction:
    out = Fraction(0)
    for c in reversed(coeffs):
        out = out * x + c
    return out


def compose_on_interval(coeffs: Sequence[Fraction], a: Fraction, b: Fraction) -> list[Fraction]:
    h = b - a
    n = len(coeffs) - 1
    out = [Fraction(0) for _ in range(n + 1)]
    for i, ci in enumerate(coeffs):
        if ci == 0:
            continue
        for j in range(i + 1):
            out[j] += ci * Fraction(math.comb(i, j)) * (a ** (i - j)) * (h ** j)
    return out


def power_to_bernstein_unit(power: Sequence[Fraction]) -> list[Fraction]:
    n = len(power) - 1
    if n < 0:
        return []
    out: list[Fraction] = []
    for i in range(n + 1):
        bi = Fraction(0)
        for j in range(i + 1):
            bi += power[j] * Fraction(math.comb(i, j), math.comb(n, j))
        out.append(bi)
    return out


def bernstein_coeffs_on_interval(
    coeffs: Sequence[Fraction], a: Fraction, b: Fraction
) -> list[Fraction]:
    return power_to_bernstein_unit(compose_on_interval(coeffs, a, b))


@dataclass(frozen=True)
class BernsteinResult:
    ok: bool
    intervals_checked: int
    failed_interval: tuple[Fraction, Fraction] | None
    min_bernstein: Fraction
    min_accepted_bernstein: Fraction | None
    max_depth_reached: int


def verify_nonnegative_bernstein(
    coeffs: Sequence[Fraction], *, lo: Fraction, hi: Fraction, max_depth: int
) -> BernsteinResult:
    if not lo < hi:
        raise ValueError("need lo < hi")
    intervals_checked = 0
    min_seen: Fraction | None = None
    min_accepted: Fraction | None = None
    max_depth_seen = 0
    stack: list[tuple[Fraction, Fraction, int]] = [(lo, hi, 0)]

    while stack:
        a, b, depth = stack.pop()
        intervals_checked += 1
        max_depth_seen = max(max_depth_seen, depth)
        bern = bernstein_coeffs_on_interval(coeffs, a, b)
        local_min = min(bern)
        min_seen = local_min if min_seen is None else min(min_seen, local_min)
        if local_min >= 0:
            min_accepted = local_min if min_accepted is None else min(min_accepted, local_min)
            continue
        mid = (a + b) / 2
        if poly_eval(coeffs, mid) < 0 or depth >= max_depth:
            return BernsteinResult(
                False, intervals_checked, (a, b), min_seen, min_accepted, max_depth_seen
            )
        stack.append((mid, b, depth + 1))
        stack.append((a, mid, depth + 1))

    return BernsteinResult(
        True,
        intervals_checked,
        None,
        min_seen or Fraction(0),
        min_accepted,
        max_depth_seen,
    )


def collision_bounds_at(
    rho: Fraction,
    lower: Sequence[Fraction],
    upper: Sequence[Fraction],
    M: int,
) -> tuple[Fraction, Fraction]:
    """Certified lower/upper bounds on P_k(rho) from coefficients through M."""
    rho = Fraction(rho)
    if not -1 <= rho <= 1:
        raise ValueError("rho must lie in [-1,1]")
    if len(lower) <= M or len(upper) <= M:
        raise ValueError("coefficient arrays do not contain enough terms")

    residual = Fraction(1) - sum(lower[: M + 1], Fraction(0))
    if residual < 0:
        raise ValueError("sum of lower coefficient bounds exceeds 1")

    if rho >= 0:
        p_lo = sum((lower[m] * rho**m for m in range(M + 1)), Fraction(0))
        p_hi = sum((upper[m] * rho**m for m in range(M + 1)), Fraction(0))
        p_hi += rho ** (M + 1) * residual
        return p_lo, p_hi

    p_lo = Fraction(0)
    p_hi = Fraction(0)
    for m in range(M + 1):
        rm = rho**m
        if m % 2 == 0:
            p_lo += lower[m] * rm
            p_hi += upper[m] * rm
        else:
            p_lo += upper[m] * rm
            p_hi += lower[m] * rm
    tail = abs(rho) ** (M + 1) * residual
    return p_lo - tail, p_hi + tail


def mixture_collision_bounds_at(
    rho: Fraction,
    mixture: Mapping[int, Fraction],
    table: CoefficientTable,
    M: int,
) -> tuple[Fraction, Fraction]:
    lo = Fraction(0)
    hi = Fraction(0)
    for k, w in mixture.items():
        if k == 1:
            klo = khi = Fraction(1)
        else:
            klo, khi = collision_bounds_at(rho, table.c_lower[k], table.c_upper[k], M)
        lo += w * klo
        hi += w * khi
    return lo, hi
