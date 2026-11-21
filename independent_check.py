"""Independent (non-rigorous) check of the two-shift witness.

Evaluates the CMQR (1.1) score directly from the definitions:
  N = 2*pi*( int_{-inf}^0 Fhat(t) e^{pi t} dt - int_0^inf |Fhat(t)| e^{pi t} dt )
  D = int_R |F(x)| dx
using Fhat as defined in the note (t-variable, no change of variables), and
F obtained (a) from the closed form and (b) by direct numerical Fourier
inversion of Fhat at sample points, to test the closed form.
Does not reuse the H(z)/omega reduction or the lambda^{-2n} identity.
"""
import json
import sys
from mpmath import mp, mpf, mpc, exp, pi, quad, factorial, findroot, fabs, sqrt, inf, linspace, cos, sin

mp.dps = 30
w = json.load(open(sys.argv[1]))
a = mpf(w["scale"])
cs = [mpf(c) for c in w["shifts"]]
B = [[mpf(b) for b in row] for row in w["coefficients"]]


def P(j, u):
    return sum(b * u**(2*n-1) / factorial(2*n-1) for n, b in enumerate(B[j], 1))


def Fhat(t):
    s = mpf(0)
    for j, c in enumerate(cs):
        u = (c - pi*t)/a
        if u >= 0:
            s += exp(-u) * P(j, u)
    return s


def F_closed(x):
    s = mpc(0)
    for j, c in enumerate(cs):
        for n, b in enumerate(B[j], 1):
            s += b * exp(2j*c*x) * (1 + 2j*a*x)**(-2*n)
    return a/pi * s


def F_numeric(x):
    # F(x) = int Fhat(t) e^{2 pi i x t} dt over t <= max(c)/pi
    tmax = max(cs)/pi
    brk = [-inf, -40, -20, -10, -5, -2, -1, 0] + ([tmax] if tmax > 0 else [])
    re = quad(lambda t: Fhat(t)*cos(2*pi*x*t), brk, maxdegree=10)
    im = quad(lambda t: Fhat(t)*sin(2*pi*x*t), brk, maxdegree=10)
    return mpc(re, im)


# (a) closed form vs numerical inversion
print("closed-form vs numeric Fourier inversion:")
for x in [mpf(0), mpf("0.3"), mpf("-1.7"), mpf("2.5")]:
    fc, fn = F_closed(x), F_numeric(x)
    print(f"  x={float(x):+.2f}  |diff|={float(abs(fc-fn)):.2e}  |F|={float(abs(fc)):.6f}")

# (b) numerator directly in t
neg = quad(lambda t: Fhat(t)*exp(pi*t), [-inf, -40, -20, -10, -5, -2, -1, -0.5, 0], maxdegree=10)
c1 = max(cs)
tmax = c1/pi
grid = linspace(mpf(0), tmax, 4001)
vals = [Fhat(t) for t in grid]
roots = []
for t0, t1, v0, v1 in zip(grid, grid[1:], vals, vals[1:]):
    if v0 == 0 or v0*v1 < 0:
        roots.append(findroot(Fhat, (t0, t1), solver="anderson"))
pts = [mpf(0)] + roots + [tmax]
pos = sum(quad(lambda t: fabs(Fhat(t))*exp(pi*t), [l, r]) for l, r in zip(pts, pts[1:]))
N = 2*pi*(neg - pos)
print("sign changes of Fhat on (0, c1/pi):", len(roots))
print("N =", N)

# (c) L1 norm of F on the real line, closed form, |F| integrated directly
edges = [mpf(0)] + [mpf(k)/8 for k in range(1, 8*16+1)] + [mpf(32), mpf(64), mpf(128), mpf(512), mpf(4096), inf]
half = lambda s: sum(quad(lambda x: abs(F_closed(s*x)), [l, r]) for l, r in zip(edges, edges[1:]))
D = half(1) + half(-1)
print("D =", D)
print("score N/D =", N/D, "  (D/N)^2 =", (D/N)**2)
