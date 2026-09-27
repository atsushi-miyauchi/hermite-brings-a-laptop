#!/usr/bin/env python3
"""Rigorous Hermite coefficient generation for spokes/FJ roundings.

This is the coefficient-generation part of the paper's ACB/Arb pipeline.  Given
spokes parameters k and a truncation degree M, it computes certified rational
lower and upper bounds on c_{k,m}, m=0,...,M.

The final edge-wise checker does not import this module; python-flint is needed
only at certificate-generation time.
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterator, Sequence


Poly = list[int]


def frac_str(x: Fraction) -> str:
    x = Fraction(x)
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def _require_flint():
    try:
        from flint import acb, arb, ctx
    except Exception as exc:
        raise RuntimeError(
            "Hermite coefficient generation requires python-flint (Arb/ACB)."
        ) from exc
    return acb, arb, ctx


def _frac_to_acb(acb, x: Fraction | int):
    x = Fraction(x)
    if x.denominator == 1:
        return acb(str(x.numerator))
    return acb(f"{x.numerator}/{x.denominator}")


def _frac_to_arb(arb, x: Fraction | int):
    x = Fraction(x)
    if x.denominator == 1:
        return arb(str(x.numerator))
    return arb(f"{x.numerator}/{x.denominator}")


def _arb_endpoint_fraction(x) -> Fraction:
    mant, exp = x.man_exp()
    mant = int(mant)
    exp = int(exp)
    if exp >= 0:
        return Fraction(mant * (1 << exp), 1)
    return Fraction(mant, 1 << (-exp))


def _arb_lower_fraction(x) -> Fraction:
    return _arb_endpoint_fraction(x.lower())


def _arb_upper_fraction(x) -> Fraction:
    return _arb_endpoint_fraction(x.upper())


def _floor_fraction_to_decimal(x: Fraction, digits: int) -> Fraction:
    scale = 10**digits
    return Fraction((x.numerator * scale) // x.denominator, scale)


def _ceil_fraction_to_decimal(x: Fraction, digits: int) -> Fraction:
    scale = 10**digits
    return Fraction(-((-x.numerator * scale) // x.denominator), scale)


def _arb_floor_decimal(x, digits: int) -> Fraction:
    return _floor_fraction_to_decimal(_arb_lower_fraction(x), digits)


def _arb_ceil_decimal(x, digits: int) -> Fraction:
    return _ceil_fraction_to_decimal(_arb_upper_fraction(x), digits)


def _trim_poly(p: Sequence[int]) -> Poly:
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def _poly_sub(p: Sequence[int], q: Sequence[int]) -> Poly:
    n = max(len(p), len(q))
    out = [0] * n
    for i in range(n):
        out[i] = (p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0)
    return _trim_poly(out)


def _poly_scalar_mul(p: Sequence[int], c: int) -> Poly:
    return _trim_poly([c * x for x in p])


def _poly_x_mul(p: Sequence[int]) -> Poly:
    return [0] + list(p)


def _poly_mul(p: Sequence[int], q: Sequence[int]) -> Poly:
    out = [0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if qj:
                out[i + j] += pi * qj
    return _trim_poly(out)


def _hermite_polynomials(M: int) -> list[Poly]:
    if M < 0:
        raise ValueError("M must be nonnegative")
    H: list[Poly] = [[1]]
    if M == 0:
        return H
    H.append([0, 1])
    for n in range(1, M):
        H.append(_poly_sub(_poly_x_mul(H[n]), _poly_scalar_mul(H[n - 1], n)))
    return H


def _multisets_of_size_sum(size: int, total: int, minimum: int = 0) -> Iterator[tuple[int, ...]]:
    if size == 0:
        if total == 0:
            yield ()
        return
    for first in range(minimum, total + 1):
        for rest in _multisets_of_size_sum(size - 1, total - first, first):
            yield (first,) + rest


def _multiset_count(lam: Sequence[int]) -> int:
    counts: dict[int, int] = {}
    for b in lam:
        counts[b] = counts.get(b, 0) + 1
    out = math.factorial(len(lam))
    for c in counts.values():
        out //= math.factorial(c)
    return out


def _denominator_factor(a: int, lam: Sequence[int]) -> int:
    out = math.factorial(a)
    for b in lam:
        out *= math.factorial(b)
    return out


def _S_polynomial(a: int, lam: Sequence[int], H: Sequence[Poly]) -> Poly:
    out = list(H[a])
    for b in lam:
        if b >= 1:
            out = _poly_mul(out, H[b - 1])
    return out


def _G_bounds(arb, mu: int, R: Fraction, max_j: int):
    if mu <= 0 or R <= 0:
        raise ValueError("mu and R must be positive")
    R_a = _frac_to_arb(arb, R)
    mu_a = arb(mu)
    e = (-mu_a * R_a * R_a / 2).exp()
    G = [arb(0)] * (max_j + 1)
    if max_j >= 0:
        G[0] = e / (mu_a * R_a)
    if max_j >= 1:
        G[1] = e / mu_a
    for j in range(2, max_j + 1):
        G[j] = (R_a ** (j - 1)) * e / mu_a + arb(j - 1) * G[j - 2] / mu_a
    return G


def _tail_bound_arb(arb, a: int, lam: Sequence[int], *, R: Fraction, H: Sequence[Poly]):
    S = _S_polynomial(a, lam, H)
    r = sum(1 for b in lam if b >= 1)
    G = _G_bounds(arb, r + 1, R, len(S) - 1)
    weighted = arb(0)
    for j, sj in enumerate(S):
        if sj:
            weighted += arb(abs(sj)) * G[j]
    prefactor = 2 / ((2 * arb.pi()) ** (r + 1)).sqrt()
    return prefactor * weighted


def _sq_lower_bound(lo: Fraction, hi: Fraction) -> Fraction:
    if lo <= 0 <= hi:
        return Fraction(0)
    return min(lo * lo, hi * hi)


def _sq_upper_bound(lo: Fraction, hi: Fraction) -> Fraction:
    return max(lo * lo, hi * hi)


def _poly_eval_acb(acb, p: Sequence[int], x):
    acc = acb(0)
    for coeff in reversed(p):
        acc = acc * x + acb(coeff)
    return acc


@dataclass(frozen=True)
class IntegralCertificate:
    k: int
    m: int
    a: int
    lam: tuple[int, ...]
    multiplicity: int
    finite_lo: Fraction
    finite_hi: Fraction
    tail_eps: Fraction
    full_lo: Fraction
    full_hi: Fraction
    sq_lb: Fraction
    sq_ub: Fraction
    denominator: int
    contribution: Fraction
    contribution_upper: Fraction


def _finite_integral_enclosure_acb(
    acb,
    a: int,
    lam: Sequence[int],
    *,
    R: Fraction,
    H: Sequence[Poly],
    rel_tol,
    abs_tol,
    deg_limit: int | None,
    eval_limit: int | None,
    depth_limit: int | None,
):
    r = sum(1 for b in lam if b >= 1)
    s = len(lam) - r
    S = _S_polynomial(a, lam, H)
    sign = -1 if r % 2 else 1
    sqrt_two_pi = (2 * acb.pi()).sqrt()
    sqrt_two = acb(2).sqrt()

    def integrand(z, analytic: bool):
        phi = (-z * z / 2).exp() / sqrt_two_pi
        Phi = (1 + (z / sqrt_two).erf()) / 2
        return acb(sign) * (phi ** (r + 1)) * (Phi ** s) * _poly_eval_acb(acb, S, z)

    kwargs = {}
    if rel_tol is not None:
        kwargs["rel_tol"] = rel_tol
    if abs_tol is not None:
        kwargs["abs_tol"] = abs_tol
    if deg_limit is not None:
        kwargs["deg_limit"] = deg_limit
    if eval_limit is not None:
        kwargs["eval_limit"] = eval_limit
    if depth_limit is not None:
        kwargs["depth_limit"] = depth_limit

    val = acb.integral(integrand, _frac_to_acb(acb, -R), _frac_to_acb(acb, R), **kwargs)
    return val.real


def _certify_integral_acb(
    acb,
    arb,
    k: int,
    m: int,
    a: int,
    lam: Sequence[int],
    *,
    R: Fraction,
    digits: int,
    H: Sequence[Poly],
    rel_tol,
    abs_tol,
    deg_limit: int | None,
    eval_limit: int | None,
    depth_limit: int | None,
) -> IntegralCertificate:
    finite = _finite_integral_enclosure_acb(
        acb,
        a,
        lam,
        R=R,
        H=H,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        deg_limit=deg_limit,
        eval_limit=eval_limit,
        depth_limit=depth_limit,
    )
    finite_lo = _arb_floor_decimal(finite, digits)
    finite_hi = _arb_ceil_decimal(finite, digits)
    tail = _tail_bound_arb(arb, a, lam, R=R, H=H)
    tail_eps = _ceil_fraction_to_decimal(_arb_upper_fraction(tail), digits)
    full_lo = finite_lo - tail_eps
    full_hi = finite_hi + tail_eps
    sq_lb = _sq_lower_bound(full_lo, full_hi)
    sq_ub = _sq_upper_bound(full_lo, full_hi)
    mult = _multiset_count(lam)
    denom = _denominator_factor(a, lam)
    weight = Fraction(k * mult, denom)
    return IntegralCertificate(
        k=k,
        m=m,
        a=a,
        lam=tuple(lam),
        multiplicity=mult,
        finite_lo=finite_lo,
        finite_hi=finite_hi,
        tail_eps=tail_eps,
        full_lo=full_lo,
        full_hi=full_hi,
        sq_lb=sq_lb,
        sq_ub=sq_ub,
        denominator=denom,
        contribution=weight * sq_lb,
        contribution_upper=weight * sq_ub,
    )


def compute_coefficients_acb(
    *,
    M: int,
    ks: Sequence[int],
    R: Fraction = Fraction(8, 1),
    prec: int = 256,
    digits: int = 40,
    rel_tol: str | None = None,
    abs_tol: str | None = None,
    deg_limit: int | None = None,
    eval_limit: int | None = None,
    depth_limit: int | None = None,
    progress: bool = True,
    details_csv: str | None = None,
) -> tuple[dict[int, list[Fraction]], dict[int, list[Fraction]]]:
    """Compute certified c_lower/c_upper for every requested k and m<=M."""
    if M < 0:
        raise ValueError("M must be nonnegative")
    ks = sorted(set(int(k) for k in ks if int(k) >= 2))
    if any(k < 2 for k in ks):
        raise ValueError("nontrivial spokes parameters must satisfy k>=2")
    if not ks:
        return {}, {}
    if R <= 0:
        raise ValueError("R must be positive")
    if digits < 1:
        raise ValueError("digits must be positive")

    acb, arb, ctx = _require_flint()
    old_prec = ctx.prec
    ctx.prec = int(prec)
    try:
        H = _hermite_polynomials(M)
        rel_tol_arb = None if rel_tol is None else _frac_to_arb(arb, Fraction(rel_tol))
        abs_tol_arb = None if abs_tol is None else _frac_to_arb(arb, Fraction(abs_tol))
        c_lower: dict[int, list[Fraction]] = {k: [] for k in ks}
        c_upper: dict[int, list[Fraction]] = {k: [] for k in ks}

        detail_file = None
        writer = None
        if details_csv:
            detail_path = Path(details_csv)
            detail_path.parent.mkdir(parents=True, exist_ok=True)
            detail_file = detail_path.open("w", newline="")
            writer = csv.DictWriter(
                detail_file,
                fieldnames=[
                    "k", "m", "a", "lambda", "N", "finite_lo", "finite_hi",
                    "tail_eps", "full_lo", "full_hi", "sq_lb", "sq_ub",
                    "denominator", "contribution", "contribution_upper",
                ],
            )
            writer.writeheader()

        try:
            for k in ks:
                for m in range(M + 1):
                    if m == 0:
                        ckm = Fraction(1, k)
                        c_lower[k].append(ckm)
                        c_upper[k].append(ckm)
                        if progress:
                            print(
                                f"k={k} m=0: exact c_{{k,0}}={frac_str(ckm)}",
                                file=sys.stderr,
                                flush=True,
                            )
                        continue

                    tasks = [
                        (a, lam)
                        for a in range(m + 1)
                        for lam in _multisets_of_size_sum(k - 1, m - a)
                    ]
                    if progress:
                        print(
                            f"k={k} m={m}: {len(tasks)} integrals",
                            file=sys.stderr,
                            flush=True,
                        )

                    lower = Fraction(0)
                    upper = Fraction(0)
                    for idx, (a, lam) in enumerate(tasks, start=1):
                        if progress and (idx == 1 or idx == len(tasks) or idx % 25 == 0):
                            print(
                                f"  integral {idx}/{len(tasks)} a={a} lambda={lam}",
                                file=sys.stderr,
                                flush=True,
                            )
                        cert = _certify_integral_acb(
                            acb,
                            arb,
                            k,
                            m,
                            a,
                            lam,
                            R=R,
                            digits=digits,
                            H=H,
                            rel_tol=rel_tol_arb,
                            abs_tol=abs_tol_arb,
                            deg_limit=deg_limit,
                            eval_limit=eval_limit,
                            depth_limit=depth_limit,
                        )
                        lower += cert.contribution
                        upper += cert.contribution_upper
                        if writer is not None:
                            writer.writerow(
                                {
                                    "k": cert.k,
                                    "m": cert.m,
                                    "a": cert.a,
                                    "lambda": " ".join(map(str, cert.lam)),
                                    "N": cert.multiplicity,
                                    "finite_lo": frac_str(cert.finite_lo),
                                    "finite_hi": frac_str(cert.finite_hi),
                                    "tail_eps": frac_str(cert.tail_eps),
                                    "full_lo": frac_str(cert.full_lo),
                                    "full_hi": frac_str(cert.full_hi),
                                    "sq_lb": frac_str(cert.sq_lb),
                                    "sq_ub": frac_str(cert.sq_ub),
                                    "denominator": cert.denominator,
                                    "contribution": frac_str(cert.contribution),
                                    "contribution_upper": frac_str(cert.contribution_upper),
                                }
                            )

                    if lower > upper:
                        raise RuntimeError(f"computed lower bound exceeds upper bound for k={k}, m={m}")
                    c_lower[k].append(lower)
                    c_upper[k].append(upper)
                    if progress:
                        print(
                            f"  c_lower[{k},{m}]={frac_str(lower)} ~ {float(lower):.17g}",
                            file=sys.stderr,
                        )
                        print(
                            f"  c_upper[{k},{m}]={frac_str(upper)} ~ {float(upper):.17g}",
                            file=sys.stderr,
                        )
        finally:
            if detail_file is not None:
                detail_file.close()

        return c_lower, c_upper
    finally:
        ctx.prec = old_prec


def write_coefficient_json(
    path: str | Path,
    *,
    M: int,
    ks: Sequence[int],
    R: Fraction,
    prec: int,
    digits: int,
    rel_tol: str | None,
    abs_tol: str | None,
    deg_limit: int | None,
    eval_limit: int | None,
    depth_limit: int | None,
    c_lower: dict[int, list[Fraction]],
    c_upper: dict[int, list[Fraction]],
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema": "fj-hermite-coefficient-bounds-v1",
        "method": "acb.integral",
        "M": int(M),
        "ks": [int(k) for k in sorted(ks)],
        "R": frac_str(R),
        "prec": int(prec),
        "digits": int(digits),
        "rel_tol": rel_tol,
        "abs_tol": abs_tol,
        "deg_limit": deg_limit,
        "eval_limit": eval_limit,
        "depth_limit": depth_limit,
        "c_lower": {
            str(k): [frac_str(x) for x in c_lower[k]] for k in sorted(c_lower)
        },
        "c_upper": {
            str(k): [frac_str(x) for x in c_upper[k]] for k in sorted(c_upper)
        },
    }
    with path.open("w") as f:
        import json
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
