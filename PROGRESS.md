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
