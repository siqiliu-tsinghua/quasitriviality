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
- [ ] Collect Dubrovin–Zhang free-energy examples, align gauges and compare
  quasi-Miura transformations order by order.
- [ ] Consider nondiagonal input and automatic Riemann coordinates only after
  the diagonal pipeline is stable.
