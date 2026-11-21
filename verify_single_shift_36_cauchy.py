#!/usr/bin/env python3
"""Certify C(1)^(-2) < 0.760777857 for an explicit single-shift witness.

Run: python verify_single_shift_36_cauchy.py
Dependencies tested: Python 3, SymPy 1.14.0 and mpmath 1.3.0.

For the Fourier extremal problem (1.1) of Carneiro, Milinovich,
Quesada-Herrera and Ramos (arXiv:2404.08380), put R(z)=sum_j r_j z^j
and let d_k be the coefficient of w^(k-1) in R(2w-1). With
Q(u)=sum_{k>=1} d_k u^k/k!, take Fhat(t)=exp(-u) Q(u) for
u=(c-pi*t)/a >= 0 and zero elsewhere. Its extremal quotient is

    2*pi*a*exp(c)*N/D,
    N = integral_0^infinity exp(-(1+a)*u) Q(u) du
        - 2*integral_0^(c/a) exp(-(1+a)*u) max(Q(u),0) du,
    D = (1/2)*integral_0^pi |R(exp(i*theta))| dtheta.

Indeed, inverse Fourier transformation gives F(x)=(a/pi) exp(2icx)
times the sum of d_k/(1+2iax)^(k+1). Substitution 2ax=tan(theta)
gives ||F||_1=D/pi; F is integrable since k starts at one.

The numerator is integrated in closed form after exact rational root
isolation. A panelwise Cauchy-Schwarz inequality bounds D from above.
All real-valued enclosures use mpmath's interval context. The root
isolations and polynomial transformations use exact arithmetic.
"""

from fractions import Fraction
from math import comb, factorial

import sympy as sp
from mpmath import iv


# Integer numerators of the six-decimal witness. Its common denominator
# cancels from the quotient, so all calculations use the integer vector.
COEFFICIENTS = (
    1000000, 4815715, 8710552, 7617483, 4487466, 3434816,
    2345152, 1189056, 938349, 510409, 258500, 235574,
    72629, 71024, 46027, 3286, 28100, -2855, 6332,
    5339, -4625, 6832, -4287, 2822, -479, -732, 1601,
    -1653, 1541, -1050, 1019, -464, 1320, 215, 2070, 5064,
)
DILATION = Fraction(3, 10)
SHIFT = Fraction(4, 5)
PANEL_COUNT = 2048
TARGET = Fraction(760777857, 10**9)


def enclose(value):
    """Convert an exact rational value to a containing real interval."""
    value = Fraction(value)
    return iv.mpf(value.numerator) / value.denominator


def shifted_coefficients(coefficients):
    """Return the coefficients of R(2w-1) in ascending powers of w."""
    return [
        2**k * sum(
            coefficients[j] * comb(j, k) * (-1)**(j - k)
            for j in range(k, len(coefficients))
        )
        for k in range(len(coefficients))
    ]


def positive_intervals(coefficients, cutoff):
    """Isolate every root and check the sign on each complementary interval."""
    u = sp.Symbol("u")
    # The factor u is positive in the interior, so isolate zeros of Q(u)/u.
    poly = sp.Poly(sum(
        sp.Rational(d, factorial(k)) * u**(k - 1)
        for k, d in enumerate(coefficients, 1)
    ), u)
    end = sp.Rational(cutoff.numerator, cutoff.denominator)
    if poly.eval(0) == 0 or poly.eval(end) == 0:
        raise ArithmeticError("Unexpected endpoint root")
    roots = poly.intervals(inf=0, sup=end, eps=sp.Rational(1, 10**32))
    if len(roots) != 2 or any(mult != 1 for _, mult in roots):
        raise ArithmeticError("Expected exactly two simple roots")

    brackets = [(sp.Integer(0), sp.Integer(0))]
    brackets.extend(pair for pair, _ in roots)
    brackets.append((end, end))
    intervals = [
        iv.mpf([enclose(lo).a, enclose(hi).b]) for lo, hi in brackets
    ]
    positive = []
    for j, (left, right) in enumerate(zip(brackets, brackets[1:])):
        if left[1] >= right[0]:
            raise ArithmeticError("Root isolates must be disjoint and interior")
        sample = (left[1] + right[0]) / 2
        sign = sp.sign(poly.eval(sample))
        if sign != (1, -1, 1)[j]:
            raise ArithmeticError("Unexpected polynomial sign in a root-free gap")
        if sign > 0:
            positive.append((intervals[j], intervals[j + 1]))
    return positive


def exponential_tail(x, coefficients, rate):
    """Enclose the integral of exp(-rate*u) Q(u) over u >= x."""
    y = enclose(rate) * x
    term = partial = iv.mpf(1)
    total = iv.mpf(0)
    for k, d in enumerate(coefficients, 1):
        term = term * y / k
        partial += term
        total += enclose(Fraction(d) / rate**(k + 1)) * partial
    return iv.exp(-y) * total


def enclose_numerator():
    coefficients = shifted_coefficients(COEFFICIENTS)
    rate = 1 + DILATION
    positive = positive_intervals(coefficients, SHIFT / DILATION)
    numerator = exponential_tail(iv.mpf(0), coefficients, rate)
    for lo, hi in positive:
        positive_mass = (
            exponential_tail(lo, coefficients, rate)
            - exponential_tail(hi, coefficients, rate)
        )
        numerator -= 2 * positive_mass
    if not (numerator.a > 0):
        raise ArithmeticError("Failed to establish a positive numerator")
    return numerator


def denominator_upper_bound():
    """Enclose a Cauchy-Schwarz upper bound for D using exact panel means."""
    # |R(exp(it))|^2 = corr[0] + 2*sum_{j>=1} corr[j]*cos(j*t).
    corr = [
        sum(x * y for x, y in zip(COEFFICIENTS, COEFFICIENTS[j:]))
        for j in range(len(COEFFICIENTS))
    ]
    width = iv.pi / PANEL_COUNT
    weights = []
    for j in range(1, len(corr)):
        half_angle = j * width / 2
        weights.append(2 * corr[j] * iv.sin(half_angle) / half_angle)

    upper = iv.mpf(0)
    for i in range(PANEL_COUNT):
        midpoint = (2*i + 1) * width / 2
        mean_square = corr[0] + sum(
            weight * iv.cos(j * midpoint) for j, weight in enumerate(weights, 1)
        )
        if not (mean_square.a > 0):
            raise ArithmeticError(f"Could not certify positivity on panel {i}")
        # Integral_I |R| <= width * sqrt(mean_I |R|^2).
        upper += width * iv.sqrt(mean_square) / 2
    return upper


def main():
    iv.dps = 65
    numerator = enclose_numerator()
    upper = denominator_upper_bound()
    score_bound = 2 * iv.pi * enclose(DILATION) * iv.exp(enclose(SHIFT)) * numerator / upper
    if not (score_bound.a > 0):
        raise ArithmeticError("Failed to bound the quotient away from zero")
    coefficient_bound = 1 / score_bound**2
    if not (coefficient_bound.b < enclose(TARGET).a):
        raise ArithmeticError(f"Certificate did not reach {TARGET}")
    print("Numerator enclosure:", numerator)
    print("Enclosure of a denominator upper bound:", upper)
    print("Enclosure of a quotient lower bound:", score_bound)
    print("Enclosure of a coefficient upper bound:", coefficient_bound)
    print("Certified: C(1)^(-2) < 0.760777857")


if __name__ == "__main__":
    main()
