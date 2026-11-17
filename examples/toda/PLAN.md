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
