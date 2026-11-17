# Transport other Toda flows with the same transformation

## Inputs

The runner reads `G₂,G₄,G₆` from [toda-epsilon6.json](toda-epsilon6.json).
It never solves a new homological equation. Set `h=(r−s)/4=e^(U/2)`,
`W=h²`, `V=(r+s)/2`, and `c=V²/2+W). Use the following *leading
speeds* in canonical coordinates, in the time normalization of CDZ:

| Flow | Speed in r | Speed in s |
| --- | --- | --- |
| Ordinary `t^(2,1)` | `h(V+h)` | `−h(V−h)` |
| Ordinary `t^(2,2)` | `VW+hc` | `VW−hc` |
| Logarithmic `t^(1,1)` | `V+hU` | `V−hU` |

The leading characteristic uses these speeds times `r_x,s_x).
The first two pairs are the dispersionless limits of the ordinary
Lax flows with powers 2 and 3; the third is the dispersionless limit
of the explicit nontranslation extended flow.

## Independent reference and comparison

On the transport side we compute `exp(ad_G)Y₀` modulo `O(ε⁷)`.
On the reference side we multiply the original difference operators,
using `(aΛ^i)(bΛ^j)=a(Λ^i b)Λ^(i+j)`, and form
`[(L^m)_+,L]/(ε m!)` for `m=2,3`, where
`L=Λ+v+e^uΛ⁻¹`. The coefficients of `Λ^0` and `Λ⁻¹` give
`v_t` and `(e^u)_t`, respectively. This is
[CDZ (2.10), PDF p. 3](https://arxiv.org/pdf/nlin/0306060#page=3).
For the logarithmic flow we instead independently expand
[CDZ (2.30)–(2.32), PDF p. 8](https://arxiv.org/pdf/nlin/0306060#page=8),
including the Bernoulli operators. We keep the numerator through
`ε⁷` before every division by `ε`.

For the two ordinary flows one can state the matched equations compactly
without listing hundreds of Taylor monomials. With `W=e^u`,
`v_±=Λ^(±1)v`, `W₊=ΛW`, and `W₊₊=Λ²W`:

```text
u_t = (1−Λ⁻¹) P_m /(ε m!)
v_t = (Λ−1)(W Λ⁻¹Q_m)/(ε m!)

P_2 = v² + W + W₊
Q_2 = v + v₊
P_3 = v³ + W₊(2v+v₊) + W(2v+v₋)
Q_3 = v² + v v₊ + v₊² + W + W₊ + W₊₊
```

The transported equations agree with these difference expressions
expanded through `O(ε⁷)`; no such formulas are used on the transport side.

The transported canonical output is converted to
`U_t=2(r_t−s_t)/(r−s)`, `V_t=(r_t+s_t)/2`, then
`u_t=U_t`, `v_t=exp(εD/2)V_t`. The independent raw reference is
evaluated on the same centered jet inputs. Both raw components of all
three comparisons have zero symbolic coefficients at every order
`0,1,…,6`, including the nonzero odd-order raw terms in `v_t`.
There are no jet poles or formal jet logs in the transported results.
We also check the commutator with the original Toda flow directly
through the same cutoff for all three examples.

The complete raw and centered coefficients are in
[toda-symmetries-epsilon6.json](toda-symmetries-epsilon6.json).
Encoding follows the [primary Toda benchmark](README.md), with two
extra state symbols for `V,U`. In the raw table the output is
`(u_t,v_t)`, but the input jets are still those of the centered
Riemann coordinates, so no additional expansion of a raw lattice
coordinate is needed to reproduce the comparison.

```sh
python -m examples.toda.symmetries --output /tmp/toda-symmetries-epsilon6.json
python -m unittest tests.test_toda_symmetries -v
```
