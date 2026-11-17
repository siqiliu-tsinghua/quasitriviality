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
