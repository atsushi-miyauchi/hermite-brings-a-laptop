#!/usr/bin/env python3
"""Certify the quadratic Gaussian noise-stability bound in Section 5.4.3.

For rho = 3/4 and every p in [0, 1], the manuscript proves
    Lambda_rho(p) <= (3/4)*p**2 + (453/1000)*p.
This script checks the two numerical inequalities in that proof, at
r = 56/65, using exact rational intervals and alternating Taylor series.
It uses only the Python standard library and needs no coefficient files.

Write H(p) = 4*Lambda_rho(p) - 3*p**2 - (453/250)*p. Its second
derivative has a convex/concave/convex sign pattern. The tangent T at
p0 = Phi(-r) lies in the concave region. Certifying T(0) < 0 and T(1) < 0,
together with H(0) = 0 and H(1) < 0, proves the bound on the whole interval.
The Jensen lower bound on D = p0 - Lambda_rho(p0) suffices for both checks;
no bivariate Gaussian integration is performed.
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction as F
from math import factorial, isqrt

DIGITS = 30
TERMS = 40


def require(condition: bool, message: str) -> None:
    # Unlike assert, certificate checks must remain active under python -O.
    if not condition:
        raise ArithmeticError(message)


class Interval:
    """Closed rational interval; every operation encloses its exact value."""

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        require(self.lo <= self.hi, "Reversed interval")

    @staticmethod
    def cast(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.cast(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        other = self.cast(other)
        products = [x * y for x in (self.lo, self.hi)
                    for y in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        require(not self.lo <= 0 <= self.hi, "Division by an interval containing zero")
        return Interval(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * self.cast(other).reciprocal()

    def __rtruediv__(self, other):
        return self.cast(other) * self.reciprocal()

    def scaled_endpoints(self, digits):
        scale = 10**digits
        lo = (self.lo.numerator * scale) // self.lo.denominator
        hi = -((-self.hi.numerator * scale) // self.hi.denominator)
        return lo, hi, scale

    def rounded(self):
        # Outward rational rounding keeps later series denominators small.
        lo, hi, scale = self.scaled_endpoints(DIGITS)
        return Interval(F(lo, scale), F(hi, scale))

    def sqrt(self):
        require(self.lo >= 0, "Negative square root")
        scale = 10**DIGITS

        def lower_sqrt(q):
            return F(isqrt((q.numerator * scale * scale) // q.denominator), scale)

        return Interval(lower_sqrt(self.lo), lower_sqrt(self.hi) + F(1, scale))

    def decimal_enclosure(self, digits=14):
        lo, hi, scale = self.scaled_endpoints(digits)

        def fmt(n):
            sign = "-" if n < 0 else ""
            n = abs(n)
            return f"{sign}{n // scale}.{n % scale:0{digits}d}"

        return f"[{fmt(lo)}, {fmt(hi)}]"


def alternating_bounds(terms):
    """Bound an alternating series by its last odd/even partial sums.

    Each caller documents why its unsigned terms decrease to zero for all
    indices, not just for the finite prefix checked here.
    """
    total = F(0)
    last_even = last_odd = previous = None
    for n, term in enumerate(terms):
        require(term >= 0 and (previous is None or term <= previous),
                "Alternating-series terms must be nonnegative and decreasing")
        previous = term
        total += term if n % 2 == 0 else -term
        if n % 2 == 0:
            last_even = total
        else:
            last_odd = total
    require(last_even is not None and last_odd is not None, "Too few series terms")
    return Interval(last_odd, last_even)


def atan_small(x):
    require(0 <= x < 1, "arctan argument outside [0, 1)")
    # Successive unsigned terms have ratio x**2*(2*n+1)/(2*n+3) < 1.
    return alternating_bounds(x**(2*n + 1) / F(2*n + 1) for n in range(TERMS))


def gaussian_integral_series(u):
    """Enclose integral_0^1 exp(-u*t**2) dt for rational 0 <= u < 1."""
    require(0 <= u < 1, "Gaussian series argument outside [0, 1)")
    # Term ratio u*(2*n+1)/((n+1)*(2*n+3)) <= u < 1.
    return alternating_bounds(u**n / F(factorial(n) * (2*n + 1))
                              for n in range(TERMS))


def exp_minus(x):
    """Enclose exp(-x) for rational 0 <= x < 1."""
    require(0 <= x < 1, "Exponential series argument outside [0, 1)")
    # Term ratio x/(n+1) <= x < 1.
    return alternating_bounds(x**n / F(factorial(n)) for n in range(TERMS))


def certify():
    # Machin's identity and tan(arccos(3/4)/2) = 1/sqrt(7).
    pi = (16*atan_small(F(1, 5)) - 4*atan_small(F(1, 239))).rounded()
    sqrt7 = Interval(7).sqrt()
    sqrt2pi = (2*pi).sqrt()
    invsqrt7 = 1/sqrt7
    theta = (2*Interval(atan_small(invsqrt7.lo).lo,
                        atan_small(invsqrt7.hi).hi)).rounded()

    r = F(56, 65)
    c = F(453, 250)
    u = r*r/2
    p0 = (F(1, 2) - r*gaussian_integral_series(u)/sqrt2pi).rounded()
    q0 = (F(1, 2) - (r/sqrt7)*gaussian_integral_series(u/7)/sqrt2pi).rounded()

    # D >= L := theta/(2*pi) * exp(-r**2/(theta*sqrt(7))) by Jensen.
    exponent = (r*r/(theta*sqrt7)).rounded()
    exponential = Interval(exp_minus(exponent.hi).lo, exp_minus(exponent.lo).hi)
    lower_bound = (theta/(2*pi)*exponential).rounded()

    # Replacing D by L gives upper bounds on the actual tangent values.
    t0_upper = 4*p0 - 4*lower_bound + 3*p0*p0 - 8*p0*q0
    t1_upper = t0_upper + 8*q0 - 6*p0 - c

    # The manuscript's analytic estimate ensures H''(p0) < 0:
    # exp(3*r**2/7) < exp(3/7) < 7/4, so H''(p0) < 2*sqrt(7)-6 < 0.
    require(0 < r < 1 and 2*sqrt7.hi < 6, "Tangent point not certified concave")
    require(c > 1, "H(1) must be negative")
    require(t0_upper.hi < -F(4, 100000), "Could not certify T(0) < -4e-5")
    require(t1_upper.hi < -F(5, 100000), "Could not certify T(1) < -5e-5")

    return {
        "Phi(-56/65)": p0,
        "Phi(-56/(65*sqrt(7)))": q0,
        "Jensen lower-bound expression L": lower_bound,
        "Upper-bound expression for T(0)": t0_upper,
        "Upper-bound expression for T(1)": t1_upper,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    try:
        values = certify()
    except ArithmeticError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    if not args.quiet:
        for name, value in values.items():
            print(f"{name}: {value.decimal_enclosure()}")
        print("PASS: Lambda_3/4(p) <= (3/4)*p^2 + (453/1000)*p for all p in [0,1].")
        print("Resulting limiting integrality ratio bound: 401/500 = 0.802.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
