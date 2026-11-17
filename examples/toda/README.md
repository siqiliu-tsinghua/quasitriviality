# Toda and the Dubrovin–Zhang CP¹ map

This is an exact local differential-jet benchmark, modulo
`O(epsilon**7)`. We work where `d=r−s=4 exp(U₀/2)>0` and both first
jets are nonzero. It makes no convergence assertion.

## Input and canonical coordinates

Use `D=D_x`, `z=εD` and the **shift form** of [Carlet–Dubrovin–Zhang
(2.14), PDF p. 4](https://arxiv.org/pdf/nlin/0306060#page=4). The raw
coefficients are
`u_t=Σ_{k=0}^6 (−1)^k ε^k D^(k+1)v/(k+1)!`,
`v_t=Σ_{k=0}^6 ε^k D^(k+1)exp(u)/(k+1)!`,
up to the common remainder. The second line of the published (2.14)
prints `(−1)^(k+1)`, inconsistent with its own backward difference and
with (2.4); the difference gives `(−1)^k`, starting with `+v_x`.

Center the mesh: `U=u`, `V=exp(−z/2)v`. Put
`A(z)=2 sinh(z/2)/z`. Then
`U_t=D A V`, `V_t=D A exp(U)`, with
`A=1+z²/24+z⁴/1920+z⁶/322560` at the cutoff. In Riemann coordinates
`r=V+2 exp(U/2)`, `s=V−2 exp(U/2)`, let `h=(r−s)/4`. The order
`k=0,2,4,6` coefficient of the centered flow is

```text
X_k^± = (D^(k+1)(h²) ± h D^(k+1)((r+s)/2)) / (2^k (k+1)!)
```

In particular the leading speeds are `+h,−h`. Centered odd coefficients
are zero; raw odd coefficients are not discarded.

## Calculation and verification

`compute.py` implements the same normalized hodograph integrations as
`CoupledSolver` in the special sparse Laurent ring of this benchmark.
The state gap and jet monomials are exact, and no genus data or candidate
span occurs in the solve. Write `G=ε²G₂+ε⁴G₄+ε⁶G₆` with
`exp(ad_G)X₀=X`. The reducing and dressing maps in `(r,s)` are
`exp(E_G)(r,s)` and `exp(−E_G)(r,s)`, respectively. The latter takes
dispersionless to centered lattice variables. Their composition is the
identity by the inverse exponentials. The runner checks all three
homological brackets, the reducing chain rule, and the dressing chain
rule independently in centered and raw lattice coordinates (through
every raw power `0..6`). Time and spatial operators use exact rationals.

The full coefficients are in [toda-epsilon6.json](toda-epsilon6.json).
Each row is a component of a series coefficient, encoded as a list of
`[exponents, rational]` terms; each exponent is `[variable_index,
power]`. Variable 0 is `d`; 1,2 are `r_x,s_x`; 3,4 are
`r_xx,s_xx`, and so on. `reduce`, `dress`, and `generator` use
`(r,s)`; the flat tables use `(u,v)` order. The leading coordinate
terms are implicit. `free_energy["4"]` and `free_energy["6"]` are
single-component rows for the reconstructed `F₂` and `F₃`.

| Per-component terms in centered `(r,s)` | ε² | ε⁴ | ε⁶ |
| --- | ---: | ---: | ---: |
| Reducing map | 9 | 107 | 842 |
| Dressing map | 9 | 110 | 850 |

To use the reducing map on raw data, set `U=u`, `V=exp(−z/2)v`,
form its canonical coordinates, and apply `reduce`. Convert the outputs
`(r₀,s₀)` to `U₀=2 log((r₀−s₀)/4)`, `V₀=(r₀+s₀)/2`.
Conversely, `original_flat_dress` already includes the undoing of the
half shift. All compositions use the same explicit cutoff.

## Exact comparison in the tau gauge

Set `B=1/A`. Use `t¹=B V`, `t²=B² U` for the computed dressing
map. This is the linear Miura map dictated by the difference operators
on the logarithm of the tau function; see [CDZ (4.9), PDF p. 14](https://arxiv.org/pdf/nlin/0306060#page=14)
and [(6.30), PDF p. 19](https://arxiv.org/pdf/nlin/0306060#page=19).
After translating both normal fields by `exp(z/2)` the formulas are
`t¹=B v`, `t²=exp(z/2)B² u`, whose coefficients begin exactly as
[(7.41)–(7.42), CMP Example 4](https://people.sissa.it/~dubrovin/Papers/1998/dubrovin_zhang_1998_cmp.pdf#page=46):

```text
t¹ = v − z² v/24 + O(ε⁴)
t² = u + z u/2 + z² u/24 − z³ u/48 + O(ε⁴)
```

Let `E₀` be the dispersionless time derivative. At order `ε^(2g)`
the predicted tau-coordinate dressing correction is
`(D E₀ F_g, D² F_g)` in `(t¹,t²)) order. From the genus-one
determinant formula and the CP¹ `G=−U₀/24` of Example 4,
`F₁=(log|r_x s_x|−U₀)/24` up to a constant. For `F₂` use the
independent transcription of [DZ (3.10.116), PDF p. 161](https://arxiv.org/pdf/math/0108160#page=161).
SymPy finds that both coefficient differences are identically zero at
`ε²` and `ε⁴`. Moreover, integrating the computed `[ε⁴]t²` twice
with respect to the total derivative gives the same `F₂` itself.

At `ε⁶` the same constructive integration gives a tau-gauge `F₃`
with 290 Laurent monomials and jet order 7. Its order-7 part is
`d² (D⁷r/r_x³ + D⁷s/s_x³)/331776`. Both
`D² F₃=[ε⁶]t²` and `D E₀ F₃=[ε⁶]t¹` are checked, as is the
coefficient of the Toda tau equation at `ε⁶`. This is a reconstruction
from the computed map, **not** an independent comparison with a
published genus-three formula.

## Reproduce from the project root

```sh
python -m examples.toda.compute --output /tmp/toda-epsilon6.json
python -m examples.toda.compare /tmp/toda-epsilon6.json
python -m examples.toda.verify_generator_sympy /tmp/toda-epsilon6.json
```

The third command independently reconstructs the flow via the package's
physical jet bracket in SymPy. Pass `--direct` to the comparison for an
additional SymPy check in the flat fields `V` and `exp(U)`. Run
`python -m unittest tests.test_toda -v` to check regeneration of the
checked-in artifact, both literature comparisons, and agreement with the
standard SymPy solver on the first nonzero correction.
