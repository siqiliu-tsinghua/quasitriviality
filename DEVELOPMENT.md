# Development and verification notes

The transform-first, same-cutoff contract in [ARCHITECTURE.md](ARCHITECTURE.md)
is authoritative. Use supported Python 3.11+, SymPy and the standard-library
unittest runner. Run `python -m unittest discover -v` from the checkout.
See [RELEASING.md](RELEASING.md) for packaging and release checks.

## Completed integration sequence

1. Jet algebra: normalized exact expressions, separate state and jet symbols,
   formal logarithms, grading, total/evolutionary derivatives and brackets.
   Check scalar and vector chain rules, inverses and cancellation.
2. Series and scalar normalization: order-aware vector series, Lie action,
   Joel inverse-jet quadrature, generator and reducing/inverse substitutions.
   Verify reduction independently of generator construction.
3. Symmetry and obstruction API: transport arbitrary symbolic leaders by
   the same generator, verify commutation to the input cutoff, and classify
   physical-jet polynomiality separately. Add examples and regression tests.
4. Coupled hodograph and constructive primitives: port the normalized
   action, invert its forced positive-jet transport by descending exact
   integrals, and solve the affine state-only Darboux remainder when present.
   Verify in physical jets and on the symbolic Tsarev family. See
   [CONSTRUCTIVE_PLAN.md](CONSTRUCTIVE_PLAN.md) for the separate commits and
   acceptance checks for this replacement.
5. Run the complete suite, check reproducible usage examples, record
   implemented support and limitations, and commit final documentation.

The implementation was committed locally in stages. The longer benchmark
checks have separate reproduction commands linked from [README.md](README.md).

## Internal interfaces

- `jets.py`: `JetSpace(names)` with `states`, `dimension`, `jet(i,s)`,
  `log(i)`, `partial(expr,i,s)`, `dx(expr,times=1)`,
  `evolutionary(Q,expr)`, `bracket(Q,P)`, `is_homogeneous(expr,weight)`,
  `active_jets(expr)` (pairs including order zero when present).
  Module helpers `canonical(expr)` and `is_zero(expr)`.
- `series.py`: `Series(coefficients, cutoff)` with vector coefficients from
  order zero, `coefficient(k)`, `dimension`, `cutoff`; reject explicit
  coefficients beyond the cutoff, pad internal known zeros. Truncated Lie
  action requires a generator of strictly positive epsilon order. It must
  not infer an extra known order for the result.
- `scalar.py`: `scalar_primitive(space,speed,source)` solving
  `[f(u)u_x,F]=source` via Joel-type inverse-jet invariants, with bounded
  termwise Laurent–log integration and explicit unsupported errors.
- `hodograph.py`: a distinct inverse-jet coordinate space, triangular
  forward/inverse jet substitutions, normalized characteristics, and
  the transported bracket. Keep logs consistent in both directions.
- Public pipeline: validate the leading equation and grading, solve in
  increasing epsilon order, return a reusable verified transformation,
  transport admissible leaders, and distinguish polynomial, obstruction,
  conditional and undetermined results. Verify each solver witness before
  using it. Exact cancellation precedes analysis of jet dependence.

The original acceptance fixtures are in [PLAN.md](PLAN.md). Prefer small exact tests that pinpoint a
sign, chain-rule, degree or cutoff error. Save no external copies or API calls.
