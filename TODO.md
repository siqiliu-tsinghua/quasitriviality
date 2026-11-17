# TODO

## Phase 1: implemented and tested

- [x] Document the scope, formal cutoff, and direction of the transform.
- [x] Implement exact physical and hodograph jet algebra, formal logs,
  total/evolutionary derivations, brackets, and truncated Lie actions.
- [x] Compute the scalar transformation directly in the original state,
  without Hopf point normalization; return and verify both directions.
- [x] Reuse the transform for arbitrary leading symmetries and report
  jet-polynomiality conditions in physical jets.
- [x] Build the coupled Tsarev family reduction and deterministic
  homological quadratures, including the affine Darboux residual.
- [x] Cross-check the equation-only inverse against the Lean-style
  positive free-jet cocycle construction for a known generic lift.
- [x] Test independent scalar benchmarks, a negative polynomiality case,
  and a genuinely coupled manufactured example through `O(ε³)`.

## Next extensions

- [ ] Add more explicit state-space quadrature reductions and support
  for inspecting or evaluating the exact local Volterra representation.
- [ ] Increase supported domains and implement better accounting of
  coefficient singularities and parameter specializations.
- [x] Cross-check the CP¹ Toda quasi-Miura map through genus two against
  Dubrovin–Zhang, and compute and verify its genus-three extension.
- [x] Verify the same three densities in the corrected Toda DZ loop equation.
- [ ] Add a second Dubrovin–Zhang benchmark and generalize the fast
  exact Laurent arithmetic beyond the Toda state-gap coefficients.
- [ ] Consider nondiagonal input and automatic Riemann coordinates only after
  the diagonal pipeline is stable.
