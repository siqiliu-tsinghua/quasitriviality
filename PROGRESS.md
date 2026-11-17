# Progress

## 2026-11-17: Phase 1

- Created the repository and documented the mathematics, interfaces and
  staged acceptance checks before implementing.
- Read the Lean source at `6b53168`, Liu–Zhang §4–5, and the Joel Ekstrand
  reference. The Mathematica notebook was not provided; the source was
  read but not executed. The Lean checkout was not rebuilt.
- Shipped the transform-first Python/SymPy API. It computes the
  quasi-Miura transformation and its inverse, verifies the substitution,
  and transports arbitrary admissible leading symmetries. All outputs
  use the same explicit cutoff. It reports polynomiality and obstructions
  separately from the verified transformation.
- Scalar checks cover KdV, Burgers, non-Hopf speeds, logarithmic primitives,
  and a nonpolynomial symmetry from the Liu–Zhang §5 class. Coupled checks
  cover the inverse jet chart, normalized brackets, symbolic arbitrary
  Tsarev speeds, and a manufactured genuinely coupled deformation with
  two nontrivial generator orders.
- Multicomponent homological solvability is implemented in a configurable
  finite span. Unsupported inputs are reported as unresolved. The full
  smooth Darboux construction and DZ cross-validation remain follow-up work.

Run `python -m unittest discover -v` from the repository root to reproduce
verification. Changes are committed locally in documented stages.

Verification on Python 3.12.14 / SymPy 1.14.0: all **51 tests passed**
in the full discovery run. The README example was executed separately
and returned the documented cutoff and polynomiality result.

## 2026-11-17: constructive coupled revision

- Re-read the Lean normalized action, order and transport lemmas, the
  free-jet cocycle integration, and the inhomogeneous Darboux construction.
  The earlier finite-span implementation above is superseded and has been
  removed. The plan was committed first as `18be940`.
- The default coupled solver now integrates the uniquely forced highest
  transverse jet partials, then recovers own-jet powers by lower-order
  scalar transport. At degree one it checks the state-only Darboux
  remainder and returns a verified primitive. General state coefficients
  defined by the local Volterra equation are exact integral representations,
  not truncated approximations or freely chosen unknown functions.
- End-to-end scope remains transform first, including the inverse, then
  arbitrary Tsarev symmetry transport and polynomiality at the same cutoff.
  No input symmetry lift is required. The separate cocycle test recovers
  the same positive-weight Laurent–log primitive from a known generic lift.
- Full checkout run: `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -v`,
  Python 3.12.14 / SymPy 1.14.0, **69 tests passed** in 99.751 seconds.
  In addition to scalar regressions, tests cover two- and three-component
  coupling, exact Laurent–log recovery with arbitrary state coefficients,
  compatible and incompatible Darboux data, and singular basepoint rejection.
- Separately executed the README Burgers example and an end-to-end
  non-inverse-difference coupled example,
  `v=(2r+3s,r+2s)` with basepoint `(1,0)`, at `O(ε³)`. It returned
  one Volterra definition for the transformation and a polynomial
  transported generic Tsarev symmetry at the same cutoff.

## 2026-11-17: Toda / Dubrovin–Zhang benchmark

- Deposited all three requested PDFs in `/workspace/reference/toda-dubrovin-zhang`.
  Read the Toda difference and tau maps, Example 4, and the CP¹ genus-two
  formula. The PDFs are reference data, not vendored in the repository.
- Computed the complete homogeneous reducing and dressing maps from the
  equation through `O(ε⁷)`. The exact specialized Laurent backend uses the
  same hodograph quadratures as the coupled solver and no free-energy input.
  Its bracket and direct chain-rule checks included the original lattice
  variables, retaining the raw odd powers. A second SymPy physical-jet
  calculation independently verified Lie reconstruction at ε², ε⁴, ε⁶.
- In the tau Miura gauge, both components agree exactly with F₁ and F₂
  at ε² and ε⁴. Reconstructing the genus-two density gives zero difference
  from (3.10.116) itself. Reconstructing the genus-three density gives
  290 terms of jet order at most seven and verifies both tau-coordinate
  identities and the order-six tau equation. No external F₃ comparison
  is claimed.
- Complete coefficients, reproduction commands, assumptions and equations:
  [examples/toda/README.md](examples/toda/README.md). Exact byte-for-byte
  regeneration of the JSON also passed.
  Full regression suite: **72 tests passed**.

## 2026-11-17: transport other Toda symmetries

- Reused the archived `G₂,G₄,G₆`; there was no new normalization solve.
  Transported the independent dispersionless leaders of `t^(2,1)`,
  `t^(2,2)`, and the first nontranslation logarithmic flow `t^(1,1)`.
- After undoing the half shift, independently expanded the CDZ difference
  Lax equations (2.10) and the explicit logarithmic equation (2.32).
  All 42 raw component/order comparisons through epsilon^6 vanished
  symbolically, including the raw odd orders. All three Toda commutators
  vanished at the same cutoff, and no space-jet denominators survived.
- [Equations and reproduction](examples/toda/SYMMETRIES.md);
  the complete coefficient file regenerated byte for byte.
  **74 tests passed**, including rejection of a changed archived generator.

## 2026-11-17: Toda DZ loop equation through genus three

- Specialized (3.10.95) directly via the CP¹ periods and Gram matrix.
  In (3.10.102), the order-zero source should be `−W/𝒟²`, and the
  prefactor of the final block should be `+ε²`, not `−ε²/2`.
  The printed `κ₀²` should be `κ₀`; both vanish for CP¹.
- Substituted the previously verified genus-one density and the unchanged
  archived genus-two and genus-three densities. Clearing spectral
  denominators gives polynomial numerators identically zero at
  `ε⁰, ε², ε⁴`. Each original substantive misprint separately produces a
  nonzero numerator. A changed genus-three coefficient is also detected.
- [Derivation and reproduction](examples/toda/LOOP_EQUATION.md).

## 2026-11-17: A₂ / Boussinesq through genus three

- Normalized the exact primary Boussinesq flow to
  `U_t=V_x`, `V_t=3 U U_x+ε² U_xxx/4`. The two canonical speeds
  `±3 b` with `r−s=12 b³` require nonzero higher pure state derivatives
  in the normalized hodograph transport.
- Constructed the map in both directions through `O(ε⁷)` by exact
  quadratures, checking physical brackets and both direct PDE chain
  rules. Recovered `F₁,F₂,F₃` by two total-derivative primitives of the
  dressed `U` coefficient and checked the `V` and tau evolution
  identities. The archive regenerated byte for byte.
- The independent period formulation of the A₂ DZ loop equation gave
  identically zero numerators at `ε⁰, ε², ε⁴`; a changed genus-three
  coefficient fails. The independently published A₂ genus-one and
  genus-two formulas matched termwise. The generic SymPy coupled solver
  also agrees with the sparse result at the first dispersive order.
- [Formulas, assumptions, and reproduction](examples/a2/README.md).

## 2026-11-17: A₂ third-derivative perturbation screen

- For exact truncated inputs `A(U) U_xxx+B(U) U_xx U_x+C(U) U_x³`
  at ε², the triangular quasi-Miura residuals require
  `A′=B′=C=0`, then `B²=0` at ε⁴.
- An independent flat polynomial lift of the nontrivial principal
  symmetry gives the same conditions, ruling out incompatibility
  artifacts from the chosen primitive representation.
- The complete power-law answer, including zero amplitudes, is
  `(A,B,C)=(constant,0,0)`. The original `A=1/4` is included.
  The generic Tsarev lift was verified in physical jets to `O(ε⁷)`
  without poles or logs; homogeneity transfers the check to any
  constant dispersion, including negative values and zero.
- [Proof and reproduction](examples/a2/SCREEN.md).

## Toda high-order performance benchmark

The centered Toda generator and both maps were computed through
`O(epsilon^11)` using the existing hodograph solver with an isolated
packed rational chart. Total construction took 268.753 s with a peak RSS
of 247.75 MiB. The direct physical reducing identity vanishes through
epsilon^10, and all coefficients through epsilon^8 reproduced with the
repository runner. See [per-order timings and reproduction instructions](examples/toda/HIGH_ORDER.md).

The same ten-order dressing map now also passes the corrected CP¹ DZ loop
equation through genus five (ε⁸ in the loop equation). F₄ and F₅ have 1,824
and 9,199 terms. All three tau identities and two independent exact loop
residual computations passed; each misprint and a genus-five mutation
were detected. See [the extended loop check](examples/toda/LOOP_EQUATION.md).
