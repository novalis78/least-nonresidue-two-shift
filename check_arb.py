"""Verify the two-shift witness with Arb ball quadrature (python-flint).

Usage: python check_arb.py witness.json
"""
import json
import sys
from math import factorial
from pathlib import Path

from flint import acb, arb, ctx, fmpq

ctx.prec = 128
tol = arb("1e-14")


def polynomial(z, coefficients):
    value = coefficients[-1]
    for b in reversed(coefficients[:-1]):
        value = value*z + b
    return value


def run(path):
    source = json.loads(path.read_text())
    scale = arb(source["scale"])
    c, c0 = map(arb, source["shifts"])
    if not (scale > 0 and c > 0 and c0 == 0):
        raise ValueError("Expected positive scale and shifts [c, 0] with c > 0")
    left, right = [[arb(b) for b in row] for row in source["coefficients"]]
    # N: integrate a smooth upper bound for twice the positive part.
    # For real y, y+sqrt(y^2+epsilon^2) >= y+abs(y).
    eps2 = arb(fmpq(1, 10**24))
    coeff = [b/factorial(2*n-1) for n, b in enumerate(left, 1)]

    def correction(u, analytic):
        y = u * polynomial(u*u, coeff)
        return (-(1+scale)*u).exp() * (y + (y*y+eps2).sqrt(analytic=analytic))

    full = arb(0)
    for shift, row in [(c, left), (c0, right)]:
        full += 2*scale*shift.exp()*sum(b/(1+scale)**(2*n) for n, b in enumerate(row, 1))
    integral = acb.integral(correction, 0, c/scale, abs_tol=tol, rel_tol=tol)
    n_lower = (full - 2*scale*c.exp()*integral.real).lower()
    if not n_lower > 0:
        raise ArithmeticError("Positive numerator required")

    # D: H(z)H(-z) = |H(z)|^2 for real z. Respect the branch cut on
    # complex integration boxes when enclosing the analytic continuation.
    i = acb(0, 1)

    def H(z):
        v = 1/(1+i*z)**2
        return v*(polynomial(v, left) + (-i*c*z/scale).exp()*polynomial(v, right))

    def magnitude(z, analytic):
        return (H(z)*H(-z)+eps2).sqrt(analytic=analytic)

    stops = [0, 1, 2, 4, 8, 16, 32, 64]
    finite = arb(0)
    for lo, hi in zip(stops, stops[1:]):
        panels = max(32, 4*(hi-lo))
        step = arb(hi-lo)/panels
        for j in range(panels):
            part = acb.integral(magnitude, lo+j*step, lo+(j+1)*step,
                                abs_tol=tol, rel_tol=tol)
            finite += part.real
    end = arb(stops[-1])
    tail = sum(abs(b)/((2*n-1)*end**(2*n-1))
               for row in (left, right) for n, b in enumerate(row, 1))
    d_upper = ((finite+tail)/arb.pi()).upper()
    if not d_upper > 0:
        raise ArithmeticError("Positive norm upper bound required")
    bound = ((d_upper/n_lower)**2).upper()
    print("Lower bound on N:", n_lower)
    print("Upper bound on D:", d_upper)
    print("Upper bound on (D/N)^2:", bound)
    if not bound < arb("0.76079030"):
        raise SystemExit("Failed to verify 0.76079030")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
