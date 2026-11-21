"""Enclose (1.1) of arXiv:2404.08380 using outward-rounded interval arithmetic.

Run: python3 certify.py witness.json [--panels 32768]
The JSON decimal strings denote exact rational numbers. This verifier handles
the two shifts [c, 0] and a common positive scale. Only mpmath is required.
"""
import argparse
import json
import math
from pathlib import Path
from mpmath import iv

iv.dps = 35
I = iv.mpf


def positive_upper(x):
    return max(I(0), x.b)


def polynomial(x, coefficients):
    result = coefficients[-1]
    for coefficient in reversed(coefficients[:-1]):
        result = result*x + coefficient
    return result


def numerator(scale, shift, rows, panels):
    lam = 1 + scale
    signed = I(0)
    for c, row in zip((shift, I(0)), rows):
        signed += 2*scale*iv.exp(c)*sum(b/lam**(2*n) for n, b in enumerate(row, 1))
    # On w=pi*t>0 only the first copy contributes. Bound the positive part of
    # exp(-lam*u)*P(u) on [0,shift/scale]. The positive part is convex and
    # 1-Lipschitz, so the same interpolation bound as for the norm applies.
    coeff = [b/math.factorial(2*n-1) for n, b in enumerate(rows[0], 1)]
    end = shift/scale

    def value(u):
        return positive_upper(iv.exp(-lam*u)*u*polynomial(u**2, coeff))

    second = I(0)
    for n, b in enumerate(coeff, 1):
        m = 2*n-1
        bound = lam**2*end**m + 2*lam*m*end**(m-1)
        if m >= 2:
            bound += m*(m-1)*end**(m-2)
        second += abs(b)*bound
    h = end/panels
    trapezoid = (value(I(0)) + value(end))/2
    for j in range(1, panels):
        trapezoid += value(j*h)
    penalty = h*trapezoid + end*h**2*second/12
    return signed - 4*scale*iv.exp(shift)*penalty


def denominator(scale, shift, rows, panels):
    omega = -shift/scale
    # z=2*a*x. A common phase of modulus 1 has been removed.
    def amplitude(z):
        v = 1/(1 + iv.j*z)**2
        phase = iv.mpc(iv.cos(omega*z), iv.sin(omega*z))
        return abs(v*(polynomial(v, rows[0]) + phase*polynomial(v, rows[1]))).b

    def second_derivative(left):
        r = iv.sqrt(1 + I(left)**2)
        total = I(0)
        for freq, row in zip((I(0), abs(omega)), rows):
            for n, b in enumerate(row, 1):
                k = 2*n
                total += abs(b)/r**k * (freq**2 + 2*freq*k/r + k*(k+1)/r**2)
        return total.b

    integral = I(0)
    # On each panel linearly interpolate the complex amplitude. Convexity
    # bounds the norm of the interpolant by the endpoint trapezoid. The
    # integrated interpolation remainder is at most sup|H''|*h^3/12.
    slabs = [0, 1, 2, 4, 8, 16, 32, 64]
    previous = amplitude(I(0))
    for left, right in zip(slabs, slabs[1:]):
        # Derivatives decay rapidly on the outer slabs.
        count = max(128, panels // max(1, left//2))
        h = I(right-left)/count
        trapezoid = previous/2
        for j in range(1, count):
            trapezoid += amplitude(left+j*h)
        previous = amplitude(I(right))
        integral += h*(trapezoid + previous/2)
        integral += (right-left)*h**2*second_derivative(left)/12

    # |1+i*z| >= z, so each term's tail has an elementary integral.
    end = I(slabs[-1])
    for row in rows:
        for n, b in enumerate(row, 1):
            integral += abs(b)/((2*n-1)*end**(2*n-1))
    return integral/iv.pi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('witness', type=Path)
    parser.add_argument('--panels', type=int, default=32768)
    args = parser.parse_args()
    witness = json.loads(args.witness.read_text())
    scale = I(witness['scale'])
    shift, second = map(I, witness['shifts'])
    if not (scale > 0 and shift > 0 and second == 0 and args.panels > 0):
        parser.error('expected a positive scale and shifts [c, 0] with c > 0')
    rows = [[I(b) for b in row] for row in witness['coefficients']]
    if len(rows) != 2 or not all(rows):
        parser.error('expected two nonempty coefficient lists')
    num = numerator(scale, shift, rows, args.panels).a
    print('Numerator lower endpoint:', num, flush=True)
    if not num > 0:
        raise SystemExit('numerator enclosure too wide')
    norm = denominator(scale, shift, rows, args.panels).b
    bound = ((norm/num)**2).b
    print('Norm upper endpoint:', norm)
    print('Reciprocal-square upper endpoint:', bound)
    if not bound < I('0.7608'):
        raise SystemExit('did not certify the target')
    print('Certified below 0.7608')


if __name__ == '__main__':
    main()
