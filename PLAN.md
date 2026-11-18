# Original development plan

This is the historical implementation plan. The planned pipeline is complete;
see [PROGRESS.md](PROGRESS.md) for results and [TODO.md](TODO.md) for limits.

## Phase 1: a working transform-first pipeline

**Input:** an evolutionary PDE specified modulo `ε^N` and an arbitrary
admissible leading symmetry. **Output:** (1) a verified quasi-Miura
transformation, (2) the transported symmetry, both modulo the same `ε^N`,
and (3) the finite-order polynomiality result or obstruction conditions.
A nonpolynomial symmetry does not invalidate the transformation. Keep the
transformation separately inspectable and reusable for another leading
symmetry without recomputing it.

### 1A. Scalar reference implementation

1. Implement exact finite jet calculus, Laurent–log arithmetic and truncation.
2. Implement Joel’s generating-vector route with the original speed `f(u)`;
   implement the Lie action and expose the actual substitution. Establish a
   documented primitive/gauge choice.
3. Transport the symmetry with symbolic `h(u)`; extract and verify conditions
   for polynomiality in physical jets. Confirm the result scales with the
   requested cutoff and does not return unknown terms.
4. Check the known KdV and Burgers examples in Liu–Zhang §4–5. Include a
   non-Hopf speed, for example take a local polynomial point change of a known
   scalar example and compare after changing variables back. Include an
   example needing a log in the transformation and a negative example from
   §5 whose polynomiality conditions fail. A log in the transformation is
   not itself a negative test for integrability.

### 1B. Multicomponent extension and first supported slice

1. Port the Lean triangular hodograph jet maps and normalized bracket.
   Verify both substitution directions and bracket conjugacy, including
   genuinely coupled speeds. Confirm the one-component restriction recovers
   the scalar homological operator and, after matching gauge, its primitive.
2. Build a representation of generic Tsarev data with constrained
   transverse derivatives and free pure derivatives. Demonstrate closure
   under the differentiation needed at the requested cutoff.
3. Implement the constructive recurrence in
   [CONSTRUCTIVE_PLAN.md](CONSTRUCTIVE_PLAN.md). Descend the forced
   transverse derivatives and own-jet coefficients in normalized jets.
   Treat the degree-two source's Darboux step explicitly, including the
   exact local radial integral representation when a finite quadrature
   is not available. Do not use an ansatz search, or require the user to
   provide an already known symmetry lift.
4. Verify the full pipeline on a genuinely coupled manufactured example.
   Candidate: `v=(2r+s,r+2s)` on `r≠s`,
   `w=grad((A(r)+B(s))/(r−s))` with independent arbitrary functions,
   and `G=ε(rₓ,0)`. Manufacture the input via the truncated
   `exp(ad_G)X₀`; the first correction is
   `(−rₓsₓ,+rₓsₓ)`. Pass criteria are the reconstructed input,
   the reduction identity and commutation for the symbolic family, not
   recovering the identical generator in a different gauge. Also verify
   that an unresolved coupled problem is not reported as nonintegrable.

### Acceptance checks across both tracks

- Test `N=1`: identity substitution, leading symmetry and no known
  dispersive correction. Test a nonconsecutive series where known internal
  zeros do not move the cutoff. Test the highest known coefficient and
  require that missing input terms are never fabricated.
- Check grading, signs, inverse composition where exposed,
  `exp(ad_G)X₀=X`, the reduction identity in physical jets, and
  `[X,Y]=0` through the stated cutoff.
- Confirm merged or exactly cancelled denominators/logs do not raise a
  false obstruction; a certified surviving one must be reported with its
  component and order. Coefficient denominators depending only on the
  states are allowed on an appropriate domain.
- Reuse one transformation with different admissible leading terms.
  For a universal claim preserve arbitrary-function dependence rather
  than substituting a handful of test speeds.
- Document the input domain, supported cases, exact verification and any
  unresolved calculation. A finite test is not a proof at every order.

## Dubrovin–Zhang cross-validation

Completed with the [Toda / CP¹](examples/toda/README.md) and
[A₂ Boussinesq](examples/a2/README.md) benchmarks. These comparisons
validate matched orders and gauges; DZ structure is not a solver dependency.
