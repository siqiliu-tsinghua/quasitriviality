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
