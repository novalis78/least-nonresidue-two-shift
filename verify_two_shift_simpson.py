#!/usr/bin/env python3
"""Independently verify a two-shift Fourier witness for the CMQR constant C(1).

Install: python -m pip install python-flint==0.9.0 sympy==1.14.0
Run:     python verify_two_shift_simpson.py path/to/witness.json

For the input scale a, shifts [c, 0], coefficients b[j,n], and lambda=1+a,
put P_j(u)=sum_n b[j,n] u^(2n-1)/(2n-1)!, with n starting at 1.
The numerator of the witness quotient equals

    N = 2a sum_{j,n} exp(c_j) b[j,n]/lambda^(2n)
        - 4a exp(c) integral_0^(c/a) exp(-lambda*u) max(P_1(u),0) du.

The verifier isolates every interior sign change using exact rational
polynomials and evaluates the exponential-polynomial integrals in closed
form. Write H(z)=sum_{j,n} b[j,n] exp(i*(c_j-c)*z/a)/(1+i*z)^(2n).
The denominator is D=(1/pi) integral_0^infinity |H(z)| dz.

On a real panel of width h, Taylor interval arithmetic bounds the fourth
derivative of |H| wherever |H|^2 stays positive. The Simpson error is at
most h^5*sup|(|H|)''''|/2880. The degree-four Taylor coefficient includes
1/4!, so its absolute bound is multiplied by h^5*4!/2880 = h^5/120.
When the square root cannot be differentiated on a panel, direct interval
evaluation bounds the integral by h times the range of |H|. If a bound is
too wide, the panel is subdivided. This also handles possible real zeros.

Beyond R=96, integral_R^infinity |H(z)| dz <= sum_{j,n} |b[j,n]| /
((2n-1)*R^(2n-1)), by |1+i*z| >= z. A positive lower endpoint for N
and an upper endpoint for D give the claimed upper bound on (D/N)^2.
All polynomial root isolation is exact. Numerical enclosures use FLINT
balls at 150-bit precision. This script does not import the archived
verifiers, call FLINT's complex quadrature, or use the GRH assumption.
"""

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path

import sympy as sp
from flint import acb, acb_series, arb, arb_series, ctx, fmpq


TAIL_START = 96
SLAB_EDGES = (0, 1, 2, 4, 8, 16, 32, 64, TAIL_START)
ERROR_PER_UNIT_LENGTH = Fraction(5, 10**12)
TARGET = Fraction(760790289, 10**9)


def ball(value):
    """Enclose a rational without a binary floating-point conversion."""
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def parse_rational(value):
    if type(value) not in (str, int):
        raise ValueError("Rational inputs must be decimal strings or integers")
    return Fraction(value)


def read_witness(path):
    source = json.loads(path.read_text(encoding="utf-8"))
    scale = parse_rational(source["scale"])
    shifts = [parse_rational(x) for x in source["shifts"]]
    if scale <= 0 or len(shifts) != 2 or shifts[0] <= 0 or shifts[1] != 0:
        raise ValueError("Expected positive scale a and shifts [c, 0] with c > 0")
    rows = source["coefficients"]
    if not isinstance(rows, list) or len(rows) != 2 or not all(
        isinstance(row, list) and len(row) == 12 for row in rows
    ):
        raise ValueError("Expected two rows of 12 coefficients for this witness")
    return scale, shifts[0], [[parse_rational(x) for x in row] for row in rows]


def positive_segments(row, end):
    """Return enclosures of the endpoints where P_1 is positive in (0,end)."""
    u = sp.Symbol("u")
    # Removing the root u=0 leaves the same sign in the open interval.
    polynomial = sp.Poly(sum(
        sp.Rational(b.numerator, b.denominator * factorial(2*n - 1)) * u**(2*n - 2)
        for n, b in enumerate(row, 1)
    ), u)
    last = sp.Rational(end.numerator, end.denominator)
    if polynomial.is_zero:
        return [], 0
    if polynomial.eval(0) == 0 or polynomial.eval(last) == 0:
        raise ArithmeticError("Unexpected endpoint root of P_1(u)/u")
    zeros = polynomial.intervals(inf=0, sup=last, eps=sp.Rational(1, 10**32))
    if any(mult != 1 for _, mult in zeros):
        raise ArithmeticError("This certificate expects simple interior roots")
    brackets = [(sp.Integer(0), sp.Integer(0)), *(pair for pair, _ in zeros), (last, last)]
    segments = []
    for (lo1, hi1), (lo2, hi2) in zip(brackets, brackets[1:]):
        if hi1 >= lo2:
            raise ArithmeticError("Interior roots must have disjoint isolating intervals")
        sign = sp.sign(polynomial.eval((hi1 + lo2) / 2))
        if sign not in (-1, 1):
            raise ArithmeticError("Failed to determine the sign between two roots")
        if sign == 1:
            segments.append((ball(lo1).union(ball(hi1)), ball(lo2).union(ball(hi2))))
    return segments, len(zeros)


def exponential_tail(x, row, rate):
    """Enclose integral_x^infinity exp(-rate*u) P_j(u) du."""
    y = ball(rate) * x
    term = partial_sum = arb(1)
    total = arb(0)
    for n, coefficient in enumerate(row, 1):
        degree = 2*n - 1
        term = term * y / degree
        partial_sum += term
        total += ball(coefficient / rate**(degree + 1)) * partial_sum
        term = term * y / (degree + 1)
        partial_sum += term
    return (-y).exp() * total


def enclose_numerator(scale, shift, rows):
    rate = 1 + scale
    a = ball(scale)
    exp_shift = ball(shift).exp()
    numerator = 2*a * (
        exp_shift*exponential_tail(0, rows[0], rate) + exponential_tail(0, rows[1], rate)
    )
    segments, root_count = positive_segments(rows[0], shift / scale)
    for lo, hi in segments:
        positive_mass = exponential_tail(lo, rows[0], rate) - exponential_tail(hi, rows[0], rate)
        numerator -= 4*a*exp_shift*positive_mass
    if not (numerator.is_finite() and numerator.lower() > 0):
        raise ArithmeticError("Failed to certify that N is positive")
    return numerator, root_count


def amplitude_function(scale, shift, rows):
    """Return H for either a complex number or a Taylor series argument."""
    left, right = [[ball(x) for x in row] for row in rows]
    imaginary_unit = acb(0, 1)
    frequency = -ball(shift / scale)

    def evaluate(z):
        v = (1 + imaginary_unit*z)**-2
        first = second = 0
        for b1, b2 in zip(reversed(left), reversed(right)):
            first = b1 + first*v
            second = b2 + second*v
        return v * (first + (imaginary_unit*frequency*z).exp()*second)

    return evaluate


def integrate_finite_norm(evaluate):
    """Enclose integral_0^R |H| with validated adaptive real Simpson panels."""
    @lru_cache(maxsize=None)
    def point_value(x):
        return abs(evaluate(acb(ball(x))))

    simpson_count = 0
    direct_count = 0
    tolerance = ball(ERROR_PER_UNIT_LENGTH)

    def panel(lo, hi, depth_left):
        nonlocal simpson_count, direct_count
        width = ball(hi - lo)
        mid = (lo + hi) / 2
        real_interval = ball(lo).union(ball(hi))
        # The expansion is about an arbitrary point of this real interval.
        series = evaluate(acb_series([real_interval, 1], 5))
        real = arb_series([series[k].real for k in range(5)], 5)
        imag = arb_series([series[k].imag for k in range(5)], 5)
        squared_norm = real*real + imag*imag
        if squared_norm[0] > 0:
            fourth_coefficient = squared_norm.sqrt()[4]
            error = abs(fourth_coefficient) * width**5 / 120
            if error.is_finite() and error < tolerance*width:
                simpson_count += 1
                estimate = width * (point_value(lo) + 4*point_value(mid) + point_value(hi)) / 6
                return estimate + arb(0, error.upper())

        value_range = abs(series[0])
        if value_range.is_finite() and value_range.upper() - value_range.lower() < tolerance:
            direct_count += 1
            return width * value_range
        if depth_left == 0:
            raise ArithmeticError(f"Failed to enclose the panel [{lo}, {hi}]")
        return panel(lo, mid, depth_left - 1) + panel(mid, hi, depth_left - 1)

    integral = sum(
        panel(Fraction(lo), Fraction(hi), 55) for lo, hi in zip(SLAB_EDGES, SLAB_EDGES[1:])
    )
    return integral, simpson_count, direct_count


def tail_upper_bound(rows):
    """Calculate the rational triangle-inequality bound for the unscaled tail."""
    return sum(
        abs(b) / ((2*n - 1) * TAIL_START**(2*n - 1))
        for row in rows for n, b in enumerate(row, 1)
    )


def main():
    parser = argparse.ArgumentParser(description="Verify the two-shift Fourier witness by real Simpson quadrature.")
    parser.add_argument("witness", type=Path, help="path to the archive's witness.json")
    ctx.prec = 150
    ctx.cap = 5
    scale, shift, rows = read_witness(parser.parse_args().witness)
    numerator, root_count = enclose_numerator(scale, shift, rows)
    integral, simpson_count, direct_count = integrate_finite_norm(amplitude_function(scale, shift, rows))
    norm_upper = ((integral + ball(tail_upper_bound(rows))) / arb.pi()).upper()
    if not (norm_upper.is_finite() and norm_upper > 0):
        raise ArithmeticError("Failed to certify a positive, finite upper bound for D")
    bound = ((norm_upper / numerator.lower())**2).upper()
    if not (bound.is_finite() and bound < ball(TARGET)):
        raise ArithmeticError(f"The certified ratio bound does not reach {TARGET}")
    print("Interior polynomial roots:", root_count)
    print("Accepted Simpson panels:", simpson_count)
    print("Accepted direct range panels:", direct_count)
    print("Numerator enclosure:", numerator)
    print("Norm upper endpoint:", norm_upper)
    print("Squared-ratio upper endpoint:", bound)
    print("Certified (D/N)^2 < 0.760790289")


if __name__ == "__main__":
    main()
