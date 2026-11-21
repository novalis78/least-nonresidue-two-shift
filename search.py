"""Two translated Fourier test functions for arXiv:2404.08380, (1.1) and (6.4).

Work in w = pi*t, write u = (shift-w)/scale, and expand the transform in
u**(2*n-1)*exp(-u)/(2*n-1)! for u >= 0. Coefficients use this sign convention.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

SCALE = 0.246
SHIFT = 0.626
TABLE = np.array([-0.0027383, 0, 4.1716, -0.6464, 3.4098, 1.0923, -1.4628,
                  2.5377, 0.94904, -2.5121, 2.1423, -0.3164])


def quadrature(breaks, order):
    x, w = roots_legendre(order)
    nodes, weights = [], []
    for left, right in zip(breaks, breaks[1:]):
        nodes.append(left + (right-left)*(1+x)/2)
        weights.append(w*(right-left)/2)
    return np.concatenate(nodes), np.concatenate(weights)


class Functional:
    def __init__(self, shifts, degree, order):
        powers = 2 * np.arange(1, degree + 1)
        # Integrate the signed transform on the entire line analytically.
        self.signed = np.concatenate([2*SCALE*np.exp(c)/(1+SCALE)**powers for c in shifts])
        # Subtract twice its positive part on the positive half-line.
        breaks = sorted(set([0.0] + [c for c in shifts if c > 0]))
        w, weights = quadrature(breaks, order)
        self.penalty_weight = 4 * weights * np.exp(w)
        blocks = []
        for c in shifts:
            u = np.maximum((c-w)/SCALE, 0)
            terms = [np.exp(-u)*u**(k-1)/math.factorial(int(k-1)) for k in powers]
            blocks.append(np.column_stack(terms))
        self.frequency = np.hstack(blocks)
        # z = 2*a*x = tan(theta). Absorb the Jacobian into the basis.
        theta, self.norm_weight = quadrature([0, .5, 1, 1.35, 1.5, np.pi/2], order)
        self.norm_weight /= np.pi
        cosine = np.cos(theta)[:, None]
        rational = cosine**(powers-2) * np.exp(-1j*theta[:, None]*powers)
        self.space = np.hstack([rational*np.exp(1j*(c-shifts[0])/SCALE*np.tan(theta))[:, None]
                               for c in shifts])

    def evaluate(self, coefficients):
        transform = self.frequency @ coefficients
        numerator = self.signed @ coefficients - self.penalty_weight @ np.maximum(transform, 0)
        norm = self.norm_weight @ np.abs(self.space @ coefficients)
        return float(numerator/norm)


def optimize(shifts, degree, order, iterations):
    f = Functional(shifts, degree, order)
    # The slowest decaying oscillatory term is unnecessary in the second copy.
    keep = np.ones(degree*len(shifts), dtype=bool)
    if len(shifts) == 2:
        keep[degree] = False
    space = f.space[:, keep]
    frequency = f.frequency[:, keep]
    weight = f.penalty_weight / 2
    # 2*max(v,0) = v + abs(v). Normalize the remaining linear term to 1.
    linear = f.signed[keep] - weight @ frequency
    start = np.zeros(len(linear))
    start[:len(TABLE)] = TABLE
    start /= linear @ start

    # Remove nearly dependent directions before reweighted least squares.
    stacked = np.vstack((space.real * np.sqrt(f.norm_weight[:, None]),
                         space.imag * np.sqrt(f.norm_weight[:, None]),
                         frequency * np.sqrt(weight[:, None])))
    _, singular, vectors = np.linalg.svd(stacked, full_matrices=False)
    stable = singular > 2e-6 * singular[0]
    transform = vectors[stable].T / singular[stable]
    y = (vectors[stable] @ start) * singular[stable]
    constraint = linear @ transform
    y /= constraint @ y
    space = space @ transform
    frequency = frequency @ transform

    for i in range(iterations):
        v = space @ y
        h = frequency @ y
        norm = f.norm_weight @ np.abs(v)
        penalty = weight @ np.abs(h)
        ratio = (1 - penalty)/norm
        if ratio <= 0:
            raise ArithmeticError('positive starting score required')
        smoothing = max(2e-8, 1e-4 * .94**i)
        w_space = f.norm_weight / np.hypot(np.abs(v), smoothing)
        w_freq = weight / np.hypot(h, smoothing)
        gram = ((space.real.T * w_space) @ space.real
                + (space.imag.T * w_space) @ space.imag
                + (frequency.T * (w_freq/ratio)) @ frequency)
        update = np.linalg.solve(gram, constraint)
        update /= constraint @ update
        if np.linalg.norm(update-y) < 1e-11:
            break
        y = update
    coeff = np.zeros(len(keep))
    coeff[keep] = transform @ y
    return coeff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shifts', type=float, nargs='+', default=[0.85, 0.0])
    parser.add_argument('--degree', type=int, default=12)
    parser.add_argument('--order', type=int, default=192)
    parser.add_argument('--iterations', type=int, default=300)
    parser.add_argument('--out', type=Path, help='save rationalized coefficients as decimal strings')
    args = parser.parse_args()
    if args.degree < len(TABLE) or len(args.shifts) not in (1, 2) or max(args.shifts) <= 0:
        parser.error('use degree >= 12 and one or two shifts, at least one positive')
    if min(args.order, args.iterations) < 1:
        parser.error('order and iterations must be positive')
    baseline = Functional([SHIFT], len(TABLE), 512).evaluate(TABLE)
    print(f'Table 2 score {baseline:.10f}, reciprocal square {baseline**-2:.10f}')
    coeff = optimize(args.shifts, args.degree, args.order, args.iterations)
    coeff = np.round(coeff, 7)
    score = Functional(args.shifts, args.degree, 800).evaluate(coeff)
    print(f'Candidate score {score:.10f}, reciprocal square {score**-2:.10f}')
    if args.out:
        witness = dict(scale=str(SCALE), shifts=[format(c, '.15g') for c in args.shifts],
                       coefficients=[[format(b, '.7f') for b in row]
                                     for row in coeff.reshape(len(args.shifts), -1)])
        args.out.write_text(json.dumps(witness, indent=2)+'\n')


if __name__ == '__main__':
    main()
