# A₂ integrability screen for a general third-derivative perturbation

On a connected domain with `U>0), consider the exact truncated system

```text
U_t = V_x,
V_t = 3 U U_x + ε² (A(U) U_xxx + B(U) U_xx U_x + C(U) U_x³).
```

In particular the ε⁴ input coefficient is **zero**, not unknown. The
answer does not classify systems with freely adjustable corrections at
ε⁴ or later. The dependent-variable and ε conventions are those in
[the A₂ benchmark](README.md). The test follows the polynomiality of
transported symmetries in [Liu–Zhang §5](https://arxiv.org/pdf/nlin/0510019#page=16).

## Complete result in this family

The differential equations forced by the screen, in exact order, are

```text
ε²:  A′ = 0, B′ = 0, C = 0.
ε⁴:  B² = 0.
```

Hence the survivors are exactly **`(A,B,C)=(c,0,0)` for constant `c`**.
The equations already force this conclusion for arbitrary smooth
coefficient functions, before making a power-law ansatz. In particular
for `(A,B,C)=(c₁ U^m₁,c₂ U^m₂,c₃ U^m₃)` with real independent exponents,
the first screen is `c₁ m₁=c₂ m₂=c₃=0`, and the second is `c₂²=0`.
Thus `c₂=c₃=0`; if `c₁≠0` then `m₁=0`, and if `c₁=0` its exponent
is irrelevant. The exponents of zero terms are always irrelevant.
There are no omitted branches where one divides by `A`, `B`, or `C`.
The original normalization corresponds to `c=1/4`. The degenerate
`c=0` branch is the dispersionless system, not a new dispersive family.

## How the obstructions are certified

Transform first in `r,s=V±6 b³`, `U=3 b²` with `b>0` and nonzero first
jets. Both canonical components receive the same dispersive correction.
The homological solver projects onto the forced triangular hodograph
primitive and checks the remainder in physical jets. Coefficient jets
of `A,B,C` are independent symbolic variables, with the usual total
derivative rule. Through order two the remainder generates the ideal
`(A′,A″,B′,C,C′)`; as differential conditions on an interval these
reduce to `A′=B′=C=0`. Under these conditions the second-order
quasi-Miura generator is verified exactly. At order four the source is
`[G₂,X₂]/2` for the homological convention `[X₀,Gₙ]=Rₙ`; the remainder
ideal is `(B²)`. For example, the coefficient of
`r_x^(-2) s_x (r_xx)^3` in the first physical bracket residual is
`B²/(540 b²)` in the convention `[X₀,G₄]−R₄`.

An independent polynomial-symmetry certificate prevents treating a
failed primitive representation as a proof of nonintegrability.
Choose the admissible hydrodynamic symmetry
`Y₀=(V U_x+U V_x, V V_x+3 U² U_x)`, whose canonical speeds are
`V±9 b³`. The equations for its homogeneous corrections are
`[X₀,Y₂]=[Y₀,X₂]` and `[X₀,Y₄]=[Y₂,X₂]`. Their right sides are
translation-invariant in the flat coordinate `V`. So are any
polynomial solutions: differentiating with respect to `V` gives a
homogeneous centralizer element, which vanishes by the nonzero
eigenspeed gaps and self-slopes (see `SymmetryVanishing.lean` in the
project's Lean reference). The same argument forces the odd corrections
to vanish through the orders used here.

At order two, comparing jet coefficients and imposing compatibility
between each coefficient and its U derivative gives `A′=B′=C=0`.
For constants `A,B` the actual coefficient is

```text
Y₂^U = (2A/3) V_xxx + (B/3) (U_x V_xx + U_xx V_x),
Y₂^V = 3UA U_xxx + (6A+2UB) U_x U_xx + B U_x³ + (B/3) V_x V_xx.
```

At order four the `A²`, `AB`, and `B²` source terms have distinct
weights for the Euler vector `S=(2U,3V)`. Since `[S,X₀]=X₀` and
the homogeneous centralizer vanishes, each source sector must have a
separate solution of its forced weight. Thus each coefficient in the
polynomial ansatz has a forced U power, and the calculation involves
only rational linear algebra. On the normalized `B²` sector a fixed
15-row combination of jet coefficient equations cancels all unknown
coefficients and equals `−2/3`. Restoring the amplitude, it gives
`−2B²/3=0`. This rules out a polynomial lift and is not just a special
sample or an unsuccessful search.

## Positive check

Every nonzero surviving system is the classical Boussinesq flow
with a scaled dispersion constant in the [cubic Gelfand–Dickey hierarchy](https://arxiv.org/abs/nlin/0407016).
In the Lax normalization of the A₂ benchmark, put `σ²=−3c` instead of
`−3/4`. Since all series are even, the
verified A₂ quasi-Miura generator at order 2g scales by `(4c)^g`;
no real square root of a negative `c` is needed for this formal test.
Transport with the full Tsarev family, treating the pure speed
derivatives as independent and eliminating the mixed derivatives,
gives no negative jet powers or jet logarithms through ε⁶ and verifies
commutation exactly through that order. The exact generic jet calculation
uses pure speed derivatives up to order nine; it does not replace the
Tsarev family by selected test functions. By scaling, the same finite
universal check holds for every `c`; at `c=0` the transform is the identity.
This is an exhaustive classification within this restricted truncated
family, not a classification of unrestricted two-component deformations.

## Reproduce

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m examples.a2.perturbation_screen
PYTHONDONTWRITEBYTECODE=1 python -m examples.a2.perturbation --symmetry
PYTHONDONTWRITEBYTECODE=1 python -m examples.a2.generic_symmetries
PYTHONDONTWRITEBYTECODE=1 python -m examples.a2.generic_symmetries --dispersion-coefficient=-3/2
PYTHONDONTWRITEBYTECODE=1 python -m unittest tests.test_a2_screen -v
```

The first command solves the power-law branches after replaying the
independent certificate and checking both symbolic residual ideals. The
second exposes all monomial coefficients and the selected direct lift.
The third and fourth check the full generic symmetry family at the
original and at a negative constant dispersion coefficient.
