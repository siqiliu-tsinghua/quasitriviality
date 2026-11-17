# Toda / CP¹ loop-equation check

## Verification plan

1. Derive the CP¹ specialization from the period form of the Dubrovin–Zhang
   loop equation (3.10.95), using the periods and Gram matrix of Example 3.10.27.
   Record any differences from the printed (3.10.102).
2. Use the genus-one density already checked against the Toda quasi-Miura map
   and load the computed genus-two and genus-three densities from
   `toda-epsilon6.json`. Do not fit or re-solve their coefficients.
   Differentiate in flat jet coordinates by the exact chain rule, and hold the
   spectral parameter fixed in every total x derivative.
3. Compare the coefficients of `epsilon^0`, `epsilon^2`, and `epsilon^4` in
   the loop equation, with exact arithmetic and independent symbolic jets.
   The loop residual must be zero as a rational function of the spectral
   parameter, not just numerically at one or more sample points. Keep a
   reproducible command and a negative control.

## Specialization and errata

Use $V=v^1=(r+s)/2$, $U=v^2=2\log((r-s)/4)$,
$W=e^U$, $a=V-\lambda$, and
$\mathcal D=a^2-4W=(r-\lambda)(s-\lambda)$.
Here $D_x$ is total differentiation at fixed $\lambda$.
The relevant source is the [period form (3.10.95), PDF p. 150](https://arxiv.org/pdf/math/0108160#page=150);
the specialization is [Example 3.10.27, PDF pp. 152–153](https://arxiv.org/pdf/math/0108160#page=152).

The periods satisfy
`p1_lambda = 2/sqrt(𝒟)`, `p2_lambda = 0`,
`grad p1 = (a, −2)/sqrt(𝒟)`, and `grad p2 = (1, 0)`,
where the flat metric interchanges components. The $p_2$ term
drops out of the jet sums because its gradient is constant, and
the remaining Gram entry is $G^{11}=-1/2$.
The flat basis has multiplication
`e1*e1=e1, e1*e2=e2, e2*e2=W e1`.

The flat partial derivatives of the contracted period pairing, in order
$V,U$, are obtained by contracting the two $\lambda$-differentiated
gradients via the flat product and the Gram matrix. They are
`8Wa/𝒟³` and `−2W(a²+4W)/𝒟³`.
Integrating in the flat coordinates and using the normalization at
spectral infinity gives the pairing `−2W/𝒟²`. The expression printed
just above (3.10.102) is its **$U$-derivative**, not the pairing. In
(3.10.102) the order-zero source is therefore `−W/𝒟²`. Also the
printed subtraction `κ₀²/λ²` should be `κ₀/λ²` as in (3.10.95), but
`κ₀=0` here by (3.10.106) with `μ̂=diag(−1/2,1/2)`.

The other correction is in the last block of (3.10.102). Let
`b = grad(p1_lambda) = (4W, −2a)/𝒟^(3/2)`. Its square in the flat product
is `(4W(a²+4W), −16Wa)/𝒟³`. Thus the last term of (3.10.95)
becomes `+ε² Σ_(γ,k) F_(γ,k) D_x^(k+1) C^γ`, with

```text
C¹ = W [−(a²+4W) V_x + 4Wa U_x] / 𝒟³
C² = W [4a V_x − (a²+4W) U_x] / 𝒟³.
```

The original brackets and placement of $W$ inside $D_x^{k+1}$
stay the same; the printed overall `−ε²/2` is replaced by `+ε²`.

For clarity, the full corrected formula in the flat jet coordinates is
below. Write $F=\Delta F$, $F_{\gamma,k}$ for its flat jet partial,
$P=(a,-2)/\sqrt{\mathcal D}$, and $J=(a,-2)/\mathcal D$.

```text
Σ_(γ,r) F_(γ,r) [
    D_x^r J^γ + Σ_(j=1..r) binom(r,j) D_x^(j−1)(1/√𝒟) D_x^(r−j+1) P^γ
] = −W/𝒟²
    − ε²/4 Σ_(γ,k;ρ,l) [F_(γ,k;ρ,l) + F_(γ,k) F_(ρ,l)]
                  D_x^(k+1) P^γ D_x^(l+1) P^ρ
    + ε² Σ_(γ,k) F_(γ,k) D_x^(k+1) C^γ.
```

The verification uses
`F₁=(log|r_x|+log|s_x|−U)/24`, the same genus-one density already
checked against both components of the tau map. The higher densities
are read directly from the `free_energy` entries of the archived
map output; no published density is substituted for the computed
$F_2$ or $F_3$.

## Reproduce

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m examples.toda.loop_equation
PYTHONDONTWRITEBYTECODE=1 python -m unittest tests.test_toda_loop_equation -v
```

The implementation evaluates all flat partials on the canonical jets by
the chain rule. It keeps the spectral symbol $a=V-\lambda$ independent
of the gap and the positive-order jets, but sets $D_x a=V_x$.
Square roots cancel between paired period derivatives. Each residual is
cleared by a power of `𝒟` before comparing all monomials in that symbol
and the jets:

| Genus | Coefficient in the loop equation | Clearing factor | Numerator |
| --- | --- | --- | --- |
| 1 | ε⁰ | 𝒟² | 0 |
| 2 | ε² | 𝒟⁵ | 0 |
| 3 | ε⁴ | 𝒟⁸ | 0 |

In particular, the only orders verified are those determined by
`ΔF=F₁+ε² F₂+ε⁴ F₃+O(ε⁶)`. The benchmark gives no all-genus claim.

As negative controls, restoring just the printed order-zero source gives
numerator `−2 W a²` for the genus-one residual with denominator `𝒟³`.
Restoring just the printed prefactor of the last block gives nonzero
residuals at genera two and three. Changing a highest-jet term of
the archived genus-three density leaves genera one and two unchanged,
but fails the genus-three check.
