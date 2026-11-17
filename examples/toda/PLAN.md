# Toda / CP1 cross-validation through O(epsilon^7)

1. Start from the shift identities of Carlet–Dubrovin–Zhang (2.14), not
   the inconsistent printed backward Taylor sign. Keep all raw powers
   through epsilon^6. Center the shift using V=exp(−epsilon D/2)v, U=u,
   then use r=V+2 exp(U/2), s=V−2 exp(U/2). In these coordinates the leading
   speeds are (r−s)/4 and −(r−s)/4 and the flow is even.
2. Apply the same constructive hodograph and transport inversions as the
   general solver, with no free-energy data in the computation. For this
   benchmark use the exact, sparse Laurent algebra over the single state
   gap d=r−s to avoid expression swell. Establish the Lie generator,
   reduction and inverse at all three positive even orders. Verify the
   brackets and directly verify the coordinate identities through cutoff.
3. Convert the inverse map (dispersionless to dispersive) to centered
   flat variables. The raw map is obtained by v=exp(epsilon D/2)V, u=U.
   Compare in the Dubrovin–Zhang gauge t¹=B V, t²=B² U, where
   B=(epsilon D)/(2 sinh(epsilon D/2)). Expand this map and reproduce
   Example 4 (7.41)–(7.42).
4. Independently transcribe F1 and F2. Compare both components of the
   transformed inverse with (D E0 Fg, D² Fg) at epsilon^(2g), for g=1,2,
   using exact coefficientwise symbolic cancellation rather than point
   samples. The epsilon^6 map is checked against Toda itself, not
   against a published F3 formula. Also try reconstructing F3 itself in
   the tau gauge by exact total-derivative primitives, and check both
   components independently if that reconstruction succeeds.
5. Keep a reproducible runner, complete machine-readable coefficients of
   both maps and the generator, and a readable derivation and record.
   Run the original regression suite and the new benchmark checks
   before reporting. No PDF is vendored in Git.

## Reuse the Toda transformation on other flows

- Read the archived generator at epsilon orders 2, 4, and 6. Do not
  solve a new homological equation for any secondary flow.
- Specify the dispersionless canonical characteristics of the ordinary
  Toda flows with p=1 and p=2 (t^{2,1}, t^{2,2}) independently. Apply the
  truncated Lie action of the saved generator.
- Obtain the comparison expressions from the raw difference Lax product
  (2.10), not from the quasi-Miura calculation. Keep the commutator
  through epsilon^7 before dividing by epsilon; compare both raw
  components at all orders 0 through 6. Also try the first nontranslation
  extended flow (2.32), with an independently expanded reference side.
- Record symbolic residuals and changes of variables; ensure cancellation
  of jet poles in the transported flows. Add a reproducible check and
  regression test.

## High-order performance experiment

1. Extend the same centered Toda input to `O(epsilon^11)`. Generate the
   homological sources by a truncated adjoint recurrence, solve each with
   the existing hodograph quadrature, and exponentiate the same generator
   in both directions. Use a larger private sparse chart without changing
   the variable indices of the archived six-order results.
2. In one fresh single-process run, record elapsed wall time for each
   positive even coefficient, cumulative time, and the process high-water
   RSS at the end of each coefficient. Finish both map directions before
   stopping the timer for each coefficient. Keep checks intrinsic to the
   homological inversion timed; perform archive export and the additional
   direct PDE check separately.
3. Compare the first three map and generator coefficients to the exact
   archive, verify the reducing chain rule in physical jets through the
   new cutoff, and record the command, numbers, and scope of the test.
   Free energies, secondary flows, and loop identities at the new orders
   are not part of this experiment.

## CP¹ loop cross-check through genus five

- Start with the unchanged O(ε¹¹) centered dressing map. Apply the established
  shifts B and B² to its flat coordinates, extract F₂–F₅ by two exact total
  primitives, and verify the two coordinate tau identities and the t,t
  identity at each genus. F₁ is the known density checked against this map.
- Extend the previous CP¹ check of the corrected DZ equation to the
  coefficients of ε⁶ and ε⁸, using only these recovered densities. Reuse the
  verified period kernels and both corrections to (3.10.102); compare over
  the exact jet and spectral function field, not at sample values.
- Recheck the earlier three genera and include independent checks that each
  restored misprint and a changed high-genus density are detected. Publish
  the exact computed densities and record the scope, inputs and reproduction
  steps. Do not infer a genus-six identity from the cutoff.
