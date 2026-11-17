# Toda quasi-Miura benchmark through O(ε¹¹)

Run in a fresh Linux Python process from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m examples.toda.high_order --verify --output /tmp/toda-epsilon10.json
```

The benchmark uses the same centered input as the sixth-order example,
$U=u$, $V=e^{-\epsilon D_x/2}v$, and $r,s=V\pm2e^{U/2}$.
The centered shift multiplier is
$2\sinh(z/2)/z=\sum_{m\geq0}z^{2m}/(2^{2m}(2m+1)!)$.
The equations and both maps are therefore even in ε. No free energies are
used as inputs, and no free energies, loop equations, or secondary flows
are computed for this benchmark.

For the truncated generator $G$, the source at order $n$ is the input
flow coefficient minus the order-$n$ part of
$\sum_{k\geq2}\mathrm{ad}_G^k X_0/k!$.
The source uses only generators of smaller order. The original hodograph
quadratures solve it with an exact physical-jet bracket check. The two maps
are the order-$n$ coefficients of $\exp(\pm E_G)(r,s)$.
For performance the benchmark allocates a private 20-jet chart, packs
exponents, groups jet substitutions, and performs multiplication via
integer numerators over common denominators. These optimizations do not
change the quadratures or the original example's chart.

## Measured run

Linux, Python 3.12.14, AMD EPYC 9V74, one Python process. These are
wall-clock stages, each including construction of the generator and
both maps and the bracket check built into the solver. Peak RSS is the
process high-water mark through the end of the stage. Startup,
serialization and the additional physical PDE chain-rule check are untimed.

| Order | Stage (s) | Cumulative (s) | Peak RSS (MiB) | Reduce / dress terms per component |
| --- | ---: | ---: | ---: | ---: |
| ε² | 0.010 | 0.010 | 12.75 | 9 / 9 |
| ε⁴ | 0.098 | 0.108 | 13.25 | 107 / 110 |
| ε⁶ | 1.293 | 1.402 | 16.88 | 842 / 850 |
| ε⁸ | 15.218 | 16.620 | 44.72 | 4,791 / 4,811 |
| ε¹⁰ | 252.133 | 268.753 | 247.75 | 22,289 / 22,336 |

These timings are from a single run of the packed-chart recurrence, not
averages and not timings of the original tuple-based sixth-order runner.

The generator and both maps agree exactly with the existing archive at
orders 2, 4 and 6. An additional, direct substitution of the reducing map
into the centered flow gives zero symbolic residual in both components
at every order through 10; the order-10 check took about 21 seconds.
The dressing map is the inverse by construction (exponentiation of the
same derivation with opposite sign). The exported JSON stores all
rational coefficients and the jet-index basis.

## Independent DZ check of the highest derivatives

```sh
python -m examples.toda.check_dz_top_jets /tmp/toda-epsilon10.json
```

[Shadrin–Wang, Theorem 1.1](https://arxiv.org/pdf/2405.00276#page=3)
expresses the dependence of $F_g$ on the jets of order $3g-2$ in terms of
the point-theory energy weighted by $f_i^{1-g}$ for each canonical metric
entry $f_i$. The coefficient of the highest jet in the point theory is
$1/(24^g g! u_x^g)$, by the one-point intersection number. At genus 1
the same formula means the derivative of $\log |u_x|/24$ with respect
to $u_x$. For Toda the flat metric $2\,dU\,dV$ gives
$f_r=2/\Delta$, $f_s=-2/\Delta$, where $\Delta=r-s$. Thus the entire
jet-order-$3g$ slice of $D_x^2 F_g$ should be
[
  \frac{\Delta^{g-1}}{2^{g-1}24^g g!}
  \left(\frac{D_x^{3g}r}{r_x^g}
        +(-1)^{g-1}\frac{D_x^{3g}s}{s_x^g}\right).
]
Recover that slice from the two canonical dressing components by
linearizing $U=2\log((r-s)/4)$ at the dispersionless gap. Products of
lower-order corrections have jet order strictly less than $3g$. The
gauge change $B^2=\left(\epsilon D_x/(2\sinh(\epsilon D_x/2))\right)^2$
does not affect this slice: its order-$2(g-h)$ term applied to the
order-$2h$ correction has jet order at most $3h+2(g-h)<3g$ for $h<g$
(including the dispersionless term).

All five slices match exactly, with no additional terms at these jet
orders. In particular, the scalar factors multiplying
$\Delta^{g-1}$ at genera 4 and 5 are $1/63700992$ and $1/15288238080$.
Deleting one highest-jet term in the genus-5 archive triggers a
mismatch. This check does not reconstruct $F_4,F_5$ or verify the full
loop equation at those genera.
