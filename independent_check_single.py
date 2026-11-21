"""Independent (non-rigorous) check of the 36-coefficient single-shift witness.

Works in the original variables t and x: builds Fhat(t)=1_{u>=0} e^{-u} P(u),
u=(c-pi t)/a, evaluates the CMQR (1.1) numerator by quadrature in t, checks
the closed form of F against numerical Fourier inversion, and integrates |F|
over the real line in x. Does not use the theta substitution or the R(z) form
of the norm.
"""
from fractions import Fraction
from math import comb
from mpmath import mp, mpf, mpc, exp, pi, quad, factorial, findroot, fabs, inf, linspace, cos, sin

mp.dps = 30
q = (1000000, 4815715, 8710552, 7617483, 4487466, 3434816, 2345152, 1189056,
     938349, 510409, 258500, 235574, 72629, 71024, 46027, 3286, 28100, -2855,
     6332, 5339, -4625, 6832, -4287, 2822, -479, -732, 1601, -1653, 1541,
     -1050, 1019, -464, 1320, 215, 2070, 5064)
r = [Fraction(x, 10**6) for x in q]
a, c = mpf(3)/10, mpf(4)/5
n = len(r)
# d[k] for k=1..36 = coefficient of w^(k-1) in R(2w-1), computed by direct expansion
d = {k + 1: sum(r[j] * comb(j, k) * 2**k * (-1)**(j - k) for j in range(k, n)) for k in range(n)}
d = {k: mpf(v.numerator) / v.denominator for k, v in d.items()}


def P(u):
    return sum(dk * u**k / factorial(k) for k, dk in d.items())


def Fhat(t):
    u = (c - pi*t)/a
    return exp(-u)*P(u) if u >= 0 else mpf(0)


def F_closed(x):
    return a/pi * exp(2j*c*x) * sum(dk * (1 + 2j*a*x)**(-(k + 1)) for k, dk in d.items())


def F_numeric(x):
    brk = [-inf, -40, -20, -10, -5, -2, -1, 0, c/pi]
    re = quad(lambda t: Fhat(t)*cos(2*pi*x*t), brk, maxdegree=10)
    im = quad(lambda t: Fhat(t)*sin(2*pi*x*t), brk, maxdegree=10)
    return mpc(re, im)


print("closed-form vs numeric Fourier inversion:")
for x in [mpf(0), mpf("0.3"), mpf("-1.7"), mpf("2.5")]:
    fc, fn = F_closed(x), F_numeric(x)
    print(f"  x={float(x):+.2f}  |diff|={float(abs(fc - fn)):.2e}  |F|={float(abs(fc)):.6f}")

neg = quad(lambda t: Fhat(t)*exp(pi*t), [-inf, -40, -20, -10, -5, -2, -1, -0.5, 0], maxdegree=10)
tmax = c/pi
grid = linspace(mpf(0), tmax, 4001)
vals = [Fhat(t) for t in grid]
roots = [findroot(Fhat, (t0, t1), solver="anderson")
         for t0, t1, v0, v1 in zip(grid, grid[1:], vals, vals[1:]) if v0*v1 < 0]
pts = [mpf(0)] + roots + [tmax]
pos = sum(quad(lambda t: fabs(Fhat(t))*exp(pi*t), [lo, hi]) for lo, hi in zip(pts, pts[1:]))
N = 2*pi*(neg - pos)
print("sign changes of Fhat on (0, c/pi):", len(roots))
edges = [mpf(0)] + [mpf(k)/8 for k in range(1, 8*16 + 1)] + [mpf(32), mpf(64), mpf(128), mpf(512), mpf(4096), inf]
D = sum(quad(lambda x: abs(F_closed(s*x)), [lo, hi]) for s in (1, -1) for lo, hi in zip(edges, edges[1:]))
print("N =", N, " D =", D)
print("score N/D =", N/D, "  (D/N)^2 =", (D/N)**2)
