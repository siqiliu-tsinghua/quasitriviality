# TODO

## Development setup

- [ ] Choose supported Python and SymPy versions and a test runner.
- [ ] Define package layout, input/output notation and cutoff metadata.
- [ ] Add small fixtures derived independently of the implementation.

## Phase 1A: scalar

- [ ] Implement and test weighted physical jets, `Dₓ`, evolutionary
      derivations, formal logs and the bracket.
- [ ] Implement exact series truncation and Lie action with the physical
      `O(ε^N)` convention.
- [ ] Implement inverse jets, Joel invariants and exact
      Laurent–log primitives without a Hopf point normalization.
- [ ] Construct the graded generator and reducing substitution; record
      normalization and verify both the Lie and substitution identities.
- [ ] Transport a symbolic arbitrary `h` and verify the commutator.
- [ ] Separate polynomiality from commutation; extract order/component
      conditions after cancellation in physical jets.
- [ ] Test boundary cutoffs, non-Hopf speed, log cancellation, KdV,
      Burgers and a nonintegrable scalar fixture.

## Phase 1B: coupled

- [ ] Translate the Lean recurrences for Tsarev coefficients, hodograph
      jets and the normalized bracket; test scalar degeneration and a
      genuinely coupled chart.
- [ ] Encode pure free jets and reduce mixed derivatives with Tsarev
      equations; test against the explicit two-component family.
- [ ] Determine a supported symbolic strategy for a common primitive,
      including the degree-two Darboux residual. Document conditions
      for exact success, certified incompatibility and unresolved output.
- [ ] Verify and document the manufactured coupled fixture with arbitrary
      `A,B`, plus an unsupported or obstructed case.

## Deferred

- [ ] Collect DZ/free-energy examples, pin source formulas, align gauges
      and compare the transforms order by order.
- [ ] Extend symbolic coverage and consider automatic Riemann coordinates
      for nondiagonal inputs only after the diagonal pipeline is stable.
