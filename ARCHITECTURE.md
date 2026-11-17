# Architecture

## 1. Input, conventions and accuracy

In physical PDE notation, write the input characteristic as
`X = X₀ + ε X₁ + … + ε^(N−1) X_(N−1) + O(ε^N)`, with integer `N ≥ 1`.
Use `X₀ⁱ = vⁱ(u) uₓⁱ` in Riemann coordinates; the scalar case is
`X₀ = f(u) uₓ`. Coefficient `X_k` has differential weight `k+1`.
A jet `∂ₓ^s uⁱ` has weight `s`, a function of the states has weight
zero, and a formal log of a first jet has weight zero. Lean and the
Liu–Zhang paper put a common ε in front of this characteristic: their
coefficient of `ε^(k+1)` is our `X_k`. Thus our `O(ε^N)` corresponds
to their `O(ε^(N+1))` for the rescaled characteristic.

A leading symmetry is `Y₀ⁱ = wⁱ(u) uₓⁱ`. For one component `w=h(u)`
is an arbitrary smooth function. In several components it is any
Tsarev speed: `∂ⱼwⁱ = Γᵢⱼ(wʲ−wⁱ)` for `i≠j), with
`Γᵢⱼ = (∂ⱼvⁱ)/(vʲ−vⁱ)`. For universal tests represent this constrained
family symbolically, for example by values and free pure derivatives
with transverse derivatives reduced using Tsarev's equations. A particular
admissible `w` is useful for testing, but certifies only that choice.
For the coupled case assume a semihamiltonian diagonal leading system,
pairwise distinct velocities and nonvanishing `∂ᵢvⁱ` on the domain.
Scalar assumptions reduce to `f′(u) ≠ 0`. General matrix input and finding
Riemann coordinates are not part of the initial interface.

An order-aware series stores the **first unknown** order. All requested
results are modulo `ε^N`; no series result claims a coefficient outside
the known range. Carry exact assumptions and formal-domain information
alongside the result.

## 2. Jet calculus and characteristic transport

The initial representation is SymPy expressions with distinct symbols
for each component and jet order. The physical jet ring is polynomial in
all positive-order jets, over local smooth functions of the state. The
intermediate ring also permits inverses of the first jets `pᵢ=uₓⁱ` and
polynomial dependence on formal `ℓᵢ=log|pᵢ|`, with
`∂ℓᵢ/∂pᵢ=1/pᵢ`. Do not use branch-sensitive log rewrites without an
explicit assumption or a formal rule.

Implement total differentiation `Dₓ` and evolutionary derivations
`E_Q = Σᵢ,s Dₓ^s(Qⁱ) ∂/∂uⁱ_(s)`, where `s=0` acts on the states.
Use the same bracket convention as the reference code:
`[Q,P] = E_Q(P) − E_P(Q)`. Truncated convolution and Lie exponential
should avoid computing terms at or above the cutoff. Record the actual
grading, not just the largest derivative appearing.

## 3. Transformation first

Use the direction `X = exp(ad_G) X₀ mod ε^N`, with
`G = ε G₁ + … + ε^(N−1) G_(N−1) + O(ε^N)` and `G_k` of weight `k`.
In this convention the reducing differential substitution is
`v = T(u) = exp(E_G)u`; its inverse is represented by `exp(−E_G)`.
Check the direction by verifying `E_X(T(u)) = X₀[T(u)] mod ε^N`, where
the notation on the right includes differential prolongation of `T`.
Expose the reducing substitution as a first-class output, together with
its generator. Preserve the selected normalization or gauge so examples
can be compared reproducibly; different valid gauges need not yield
identical formulas.

At order `k`, subtract the nested-bracket contributions already fixed
by `G₁,…,G_(k−1)` from `X_k` to obtain a residual `R_k`.
With the stated direction solve `[G_k,X₀]=R_k`
(equivalently `[X₀,G_k]=−R_k`). Verify each residual identity and
reconstruct the input by the full truncated action.

### Scalar solver

Follow Joel’s `reducingVector` and use `reducingTransformation` as a
second reference for the substitution. Do not require a point change
`u→f(u)`. In inverse jets `z_s=∂ᵤ^s x` the direction
`E_(f p)−f Dₓ` fixes `u` and acts by `z_s↦−f^(s)`.
Use `z_s−(f^(s)/f′)z₁` for `s≥2` as invariant coordinates.
Holding them and `u` fixed gives
`[f p,F]=f′ p³ ∂ₚ(F/p)`. For `[f p,F]=S` one may use
`F=(p/f′) ∫ S/p³ dp`; pass `S=−R_k` in the generator step.
Use a reproducible primitive convention, initially zero homogeneous
integration constant in the chosen invariant coordinates. Integrate
Laurent–log monomials exactly, including the exponent `−1` case.

### Coupled solver

Use the local generalized hodograph jet chart of the Lean reference:
triangular inverse jets with Tsarev transverse corrections, invertible
even with formal logs. Normalize each characteristic by dividing its
component by `pᵢ`; use the transported total and evolutionary
derivations to form the normalized bracket. This keeps the physical
speeds and does not replace a coupled homological system by independent
copies of the scalar integral.

Relevant Lean sources under `SemiHamiltonian/Jet/`:
`HodographCoefficients.lean`, `HodographJets.lean`,
`HodographSubstitution.lean`, `HodographConjugacy.lean`,
`NormalizedCharacteristic.lean`, `NormalizedBracket.lean`,
`NormalizationStep.lean` and `HomologicalExactness.lean`.
Lean uses an all-Tsarev cocycle in the homological existence proof:
integrate the closed vertical form in the positive inverse jets, then
solve the remaining affine Darboux equation in weight one. The optional
`primitive_from_lift` independently reproduces this when a generic lift
is known. It is not needed by the equation-only API.

In the equation-only construction, write `A_v H = S` in the same chart.
The coefficient of the next transverse jet is the positive-jet
partial derivative of `H_i` times the speed gap over `z_j,1`. Taking
its zero-constant primitive strips transverse dependence. In the same
order, the unknown positive powers of the own jet reduce to scalar
shifted transport at lower order. Its weight stays positive; hence an
own leading power cannot disappear from the projected source. At the
bottom, solve the projected Euler–log equation, permitting just a
state-only remainder in weight zero. The full descent and multiplier
shifts are specified in [CONSTRUCTIVE_PLAN.md](CONSTRUCTIVE_PLAN.md).

The remaining state system, when needed, is
`∂_j h_i = Γ_ij (h_j − h_i) + f_ij`. The reference's radial Volterra
operator supplies an exact local definition `h = (I−T)^−1 A`, with the
basepoint, forcing, radial integrals and derivative identities retained.
On a small enough regular cube the operator is contractive; the
compatibility check allows using the affine Darboux equations as
exact differentiation rules. Decoupled cases use finite quadratures.
A special two-component finite quadrature is available with a different gauge.
This is not a right inverse on arbitrary sources: verify residuals and
the physical bracket, and keep failed construction separate from the
polynomiality test.

## 4. Symmetry and integrability conditions

After verifying the transformation, compute
`Y = exp(ad_G) Y₀ mod ε^N`. Verify the leading term and check
`[X,Y]=0 mod ε^N` in physical jets. If testing polynomiality, normalize
each coefficient in the physical variables (not the hodograph variables).
Terms with surviving negative first-jet powers, nonzero formal-log
powers or other nonpolynomial jet dependence do not belong to the
differential polynomial ring. With arbitrary admissible `w`, reduce
its transverse derivatives before declaring any symbolic coefficient an
obstruction. Independence of the free data justifies coefficient
splitting; do not assume all mixed derivatives are independent.

Return the transformation and the transported symmetry whether or not
the latter is polynomial. Optionally return order-, component- and
monomial-tagged conditions for their polynomiality, with statuses such
as verified through cutoff, obstructed, parameter conditions required,
or undetermined. An unevaluated integral or unresolved simplification
is not proof of a nonpolynomial term. Agreement to a finite order is
not a conclusion about all orders or about convergence.
