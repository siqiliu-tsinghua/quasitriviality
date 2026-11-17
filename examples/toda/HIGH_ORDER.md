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
