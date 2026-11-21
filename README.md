# A two-shift witness for the least quadratic nonresidue

Lennart Lopin (Euler's Identity LLC, ORCID [0009-0001-2751-2060](https://orcid.org/0009-0001-2751-2060))

Archived version: Zenodo DOI to follow.

For an odd prime $p$, let $n_p$ be its least positive quadratic nonresidue.
The explicit function below gives $\mathcal C(1)^{-2}<0.76079030$. Theorem 1 of
Carneiro, Milinovich, Quesada-Herrera and Ramos then gives
$\limsup_{p\to\infty} n_p/(\log p)^2<0.76079030$ under GRH for Dirichlet
$L$-functions. Their published bound is $0.7615$. The interval computation
uses no hypothesis about zeros. GRH enters in the cited theorem connecting
the Fourier optimization problem to primes. [Paper, Theorem 1 and §6](https://arxiv.org/html/2404.08380v2)

Context: L. Alpöge's [announcement](https://alpo.ge/remark.pdf) of a Lean
formalization covering Dirichlet GRH, with A. V. Sutherland's
[verification report](https://math.mit.edu/~drew/grh/), prompted this search.
The bound on $\mathcal C(1)^{-2}$ does not depend on that work. If the
Dirichlet GRH result is accepted, it discharges the hypothesis of the two
corollaries about primes. The note states both corollaries conditionally.

## Reproduce the certificate

From this directory, run the following command with Python 3 and mpmath
(tested with 1.3.0).

```sh
python3 certify.py witness.json
```

The program uses 35-decimal-place outward-rounded interval arithmetic. With
the default 32,768 panels on the first slab it returned the following bounds,
rounded outward here.

| Quantity | Bound |
| --- | ---: |
| Weighted numerator $N$ | $>0.91779528593$ |
| Norm $D$ | $<0.80053609157$ |
| $(D/N)^2$ | $<0.760799446$ |

The [coefficients](witness.json) are exact decimal rationals. The
[verifier](certify.py) reads them directly and runs independently of the
optimizer. A separate 45-digit quadrature calculation gave approximately
$0.76079029$ for the witness. The proposition in the paper uses the tighter Arb upper bound below.
The Lean file proves the tail lemma and the final algebraic comparison.
The Python verifiers establish the finite-domain enclosures.
A third, separately written verifier is described below.

## Second certificate with Arb

The same exact witness gives $(D/N)^2<0.76079030$ with a second
implementation using python-flint 0.9.0. Run the following in this directory.

```sh
python3 check_arb.py witness.json
```

The [Arb verifier](check_arb.py) uses complex ball quadrature at 128-bit
precision. In the numerator it bounds $2\max(y,0)$ from above by
$y+\sqrt{y^2+10^{-24}}$ on the real axis. For the norm integral over
$[0,64]$, it integrates $\sqrt{H(z)H(-z)+10^{-24}}$, which bounds $|H(z)|$
from above for real $z$. Both square roots pass the analyticity flag to
Arb's integral routine so it checks the branch cut on complex boxes. The
infinite tail uses the analytic bound given below. The run returned
$N>0.9177970191817$, $D<0.8005327903879$, and $(D/N)^2<0.760790297865$,
rounded outward. Both verifiers use the same Fourier reduction.

## Third certificate with exact roots and Simpson panels

`verify_two_shift_simpson.py` shares no code with the first two verifiers.
It isolates the sign changes of $P_1$ by exact rational root isolation
(SymPy), integrates the positive part in closed form, and bounds the norm
with adaptive Simpson panels on $[0,96]$, where FLINT Taylor balls at
150 bits bound the fourth derivative of $|H|$. Run it with python-flint 0.9.0
and SymPy 1.14.0.

```sh
python3 verify_two_shift_simpson.py witness.json
```

It reports $(D/N)^2<0.760790289$.

## Non-rigorous cross-check

`independent_check.py` (mpmath only) evaluates (1.1) from the definitions,
without the closed forms in Section 2. It compares the closed form for $F$
with numerical Fourier inversion of $\widehat F$, integrates
$|\widehat F|e^{\pi t}$ in $t$ between its sign changes, and integrates
$|F|$ over the real line. It is a floating-point sanity check, not a proof.
It gives $(D/N)^2\approx0.7607902846$.

Theorem 2 of the cited paper also bounds the limsup for the least *prime*
quadratic residue by $\mathcal C(1)^{-2}$ under Dirichlet GRH. It therefore
inherits the same improved upper bound.

## Witness and normalization

Use the Fourier convention $\widehat F(t)=\int F(x)e^{-2\pi ixt}\,dx$.
Set $a=0.246$, $c_1=0.85$, $c_2=0$ and write $b_{jn}$ for row $j$,
entry $n$ of the coefficient file, starting at $n=1$. Put

$$
P_j(u)=\sum_{n=1}^{12} b_{jn}\frac{u^{2n-1}}{(2n-1)!},\qquad
\widehat F(t)=\sum_{j:u_j\ge0}e^{-u_j}P_j(u_j),\qquad
u_j=(c_j-\pi t)/a .
$$

This definition fixes the signs of the odd powers by using a nonnegative
variable on each component's support. Inverting the transform termwise by the
gamma integral gives

$$
F(x)=\frac a\pi\sum_{j,n} e^{2ic_j x} b_{jn}(1+2iax)^{-2n}.
$$

In particular $F$ is integrable, its Fourier transform is real, and the
positive numerator certified above makes it nonzero. These are the
admissibility conditions for (1.1) of the paper. In that formula the score
for this witness is $N/D$ and $\mathcal C(1)$ is the supremum of the scores.

## Numerator

Let $w=\pi t$. The numerator in (1.1) is the signed integral over all $w$
with weight $2e^w$, less twice the contribution of the positive part on
$w>0$ with that same weight. Only the component with $c_1>0$ contributes
to the correction. Integrating the signed terms gives

$$
N=2a\sum_{j,n}\frac{e^{c_j}b_{jn}}{(1+a)^{2n}}
 -4ae^{c_1}\int_0^{c_1/a}
       \max\{e^{-(1+a)u}P_1(u),0\}\,du .
$$

The verifier bounds the correction from above. For each monomial
$d u^m$ in $P_1$, put $U=c_1/a$ and $\lambda=1+a$. An upper bound on the
absolute second derivative of $d e^{-\lambda u}u^m$ over $[0,U]$ is

$$
|d|\big(\lambda^2U^m+2\lambda m U^{m-1}+m(m-1)U^{m-2}\big).
$$

Omit the last summand for $m=1$. Add these bounds and apply the panel
estimate below with the positive-part function.

## Norm and quadrature error

Set $z=2ax$, $\omega_j=(c_j-c_1)/a$, and define

$$
H(z)=\sum_{j,n} b_{jn}e^{i\omega_j z}(1+iz)^{-2n},\qquad
D=\|F\|_1=\frac1\pi\int_0^\infty|H(z)|\,dz .
$$

The norm identity uses conjugate symmetry and removal of the common phase.
On a slab starting at $l\ge0$, set $r=\sqrt{1+l^2}$ and $k=2n$. The
product rule gives a uniform second-derivative bound for each summand,

$$
|b_{jn}|r^{-k}\big(\omega_j^2+2|\omega_j|k/r+k(k+1)/r^2\big).
$$

For a panel of length $h$, the linear interpolant to a real or complex
function $g$ differs from $g$ at position $x$ within the panel by at most
$\sup|g''|\,x(h-x)/2$. This follows by integrating the interpolation
remainder kernel and holds for the complex modulus as well. Both modulus
and positive part are convex and 1-Lipschitz. Thus the integral of either
applied to $g$ is bounded above by its endpoint trapezoid plus
$h^3\sup|g''|/12$. This works across zeros of $g$.

The verifier uses this estimate on $[0,64]$, accumulating interval upper
bounds for the endpoint values. Beyond $R=64$ it bounds each remaining
integral by $|b_{jn}|R^{1-k}/(k-1)$ using $|1+iz|\ge z$. Keeping
the numerator lower endpoint and norm upper endpoint gives the stated
upper bound for $(D/N)^2$.

## Reproduce the search

The search needs NumPy and SciPy (tested with 2.3.5 and 1.18.1).
The first line of output reproduces the $A=1$ test in Table 2 with
score about $1.14600669$ and reciprocal square about $0.76142248$.
The second uses two independently weighted translates of the same basis.
The numerical optimizer normalizes the linear part of the objective and
solves successive weighted least-squares problems.

```sh
OPENBLAS_NUM_THREADS=1 python3 search.py --out /tmp/nonresidue-witness.json
python3 certify.py /tmp/nonresidue-witness.json
```

The saved witness is sufficient to repeat the interval proof across
different optimizer or BLAS versions.

## Check the Lean theorems

Use Lean and Mathlib v4.24.0. The following commands assume that the
matching Lean toolchain and Lake are on PATH. Run from the directory with
the replication files. The cache downloads prebuilt Mathlib dependencies.

```sh
note_dir="$PWD"
git clone --depth 1 --branch v4.24.0 https://github.com/leanprover-community/mathlib4.git mathlib4
cd mathlib4
lake exe cache get Mathlib/Analysis/SpecialFunctions/ImproperIntegrals.lean Mathlib/Tactic/NormNum.lean Mathlib/Tactic/Linarith.lean Mathlib/Tactic/Positivity.lean
lake env lean -DwarningAsError=true "$note_dir/WitnessBounds.lean"
```

`TwoShift.tail_bound` proves integrability and the infinite-tail inequality
for an arbitrary finite sum with the paper's phases and denominator powers.
`TwoShift.witness_tail_arithmetic` checks that the unnormalized tail
majorant computed from Table 1 at 64 is less than 0.000002035.
`TwoShift.constant_of_enclosures` transfers the rounded Arb enclosures and
the score upper bound to the stated constant. The file prints the axioms
of those declarations. The Fourier inversion, compact quadrature and
prime-number corollaries are justified in the paper and accompanying code.

## Paper source

Compile `two-shift-note.tex` with pdfLaTeX or Tectonic. It needs amsmath,
amssymb, booktabs, geometry and hyperref. `two-shift-note.pdf` is the
compiled version deposited on Zenodo.

## Verification record (2026-11-20)

| Check | Result |
| --- | --- |
| `check_arb.py` (python-flint 0.9.0) | $(D/N)^2<0.760790297865$ |
| `certify.py` (mpmath 1.3.0) | $(D/N)^2<0.760799446$ |
| `verify_two_shift_simpson.py` (python-flint 0.9.0, SymPy 1.14.0) | $(D/N)^2<0.760790289$ |
| `independent_check.py` (floating point) | $(D/N)^2\approx0.7607902846$ |
| `WitnessBounds.lean` (Lean/Mathlib v4.24.0, `-DwarningAsError=true`) | compiles; axioms `propext`, `Classical.choice`, `Quot.sound` |

## Provenance and licence

OpenAI Codex assisted with the numerical search, verification code and the
note. Anthropic's Claude independently reran all verifiers and the Lean check
before deposit. Code is released under the MIT License (see `LICENSE`). The
note (`two-shift-note.tex`, `two-shift-note.pdf`) is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
