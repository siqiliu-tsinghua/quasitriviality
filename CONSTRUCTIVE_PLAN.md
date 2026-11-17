# Replace the coupled search with a constructive homological inverse

The previous finite-span solver is not the multicomponent algorithm. Remove
its candidate, state-basis and degree knobs. Preserve the public transform-first,
same-cutoff contract and the scalar Joel implementation.

## Reference and extraction route

Read in this order: Lean `NormalizedBracket`, `HodographConjugacy`,
`HomologicalOrder`, `TransportVanishing`, `CocycleClosedForm`,
`CocycleValueReduction`, `CocycleDarbouxExactness`, and
`Darboux/InhomogeneousOpenBall`. Its homological exactness theorem assumes a full
Tsarev cochain. The closed form is integrated in the positive auxiliary jets;
in weight one the value-only remainder is an inhomogeneous Darboux system.
Do not assume existence for an arbitrary source or treat unsuccessful symbolic
calculation as an integrability obstruction.

The input to this project is the **equation**, not a supplied lift of its
arbitrary symmetry. Recover the positive-jet part directly from the
single-speed normalized action, exploiting the same triangularity and
transport injectivity as Lean's order lemma:

1. Write `S = normalize(R)`, `H = normalize(F)`; solve `A_v H = S`,
   equivalent to `[v p, F] = R`. Put `m = weight(F) − 1 >= 0`.
   Each auxiliary z has weight −1; auxiliary logs have weight zero.
2. Begin at the source jet order (at least one). No finite-jet primitive
   can have higher order: the highest transverse symbol in the bracket
   determines transverse partials, and the highest own derivative obeys
   the shifted, injective homogeneous transport equation.
3. At order N, for j ≠ i, read `∂H_i/∂z_j,N` unambiguously from the
   coefficient of `z_j,N+1` in the equation, by multiplying by
   `z_j,1/(v_j−v_i)`. Use exact termwise primitives, accounting for log
   derivatives when N=1. Subtract the action of the recovered vector.
4. For N≥2, the unknown same-component order-N jet occurs polynomially.
   Its velocity is `−α_i z_i,N + β`, with `α_i = Σ_{j≠i} ∂_j v_i/z_j,1`
   and β independent of this jet. No other unknown component has this jet.
   Descend powers k>0 in this variable. Solve the coefficient by the lower
   jet order scalar transport inverse with weight m+k and multiplier
   (1−k)α_i, then subtract. The same operation recurses for scalar transport
   with coefficient cα_i (shift c−k). Positive weight makes it injective:
   at order one its pure own-power projection is the invertible polynomial
   operator `−∂_i v_i (∂log−weight)`. Thus the top power, when present, is
   detected in the source and no polynomial degree limit is guessed.
5. At N=1, after integrating transverse gradients, project onto the pure
   own first-jet power. Invert the Euler–log operator; at weight zero take
   the log primitive with zero constant. Positive weight leaves no
   unknown; at weight zero an unknown state vector h remains.
6. At weight zero, the leftover must have the form
   `Σ_{j≠i}(v_j−v_i) f_ij/z_j,1`. Check this and the affine Darboux
   compatibility identities, and solve
   `∂_j h_i = Γ_ij(h_j−h_i)+f_ij`. Use exact quadratures where available
   and otherwise the local radial Volterra construction of the reference,
   zero axis data on a specified admissible basepoint. Do not truncate
   an infinite state-space integral as if it were exact. Keep coefficient
   differentiation and locality explicit. If a coefficient backend
   cannot represent the solution, report that distinctly.
7. Denormalize and verify the bracket in physical jets. Verify
   `exp(ad_G)X_0=X`, the actual coordinate substitutions and inverse
   composition, then transport and check the generic Tsarev symmetry at
   the same cutoff. Never solve for that symmetry as a prerequisite of
   returning an independently verified transformation.

## Commit and verification sequence

- Commit this plan first.
- Commit the normalized transport primitives and order-by-order solver,
  with direct recovery of manufactured Laurent/log generators without
  any candidate input.
- Commit Darboux integration and end-to-end testing: genuinely coupled
  examples beyond the previous span, nontrivial Darboux forcing, arbitrary
  symbolic leaders, exact input cutoff, and a source outside the image.
  Compare against the Lean cocycle construction when a generic lift is known.
- Re-run the full suite, document what is computed as a finite expression
  versus as a local integral representation, record verification and commit.
