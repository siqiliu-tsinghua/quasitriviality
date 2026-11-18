# quasitriviality

Exact, finite-order quasi-Miura transformations and symmetry tests for
evolutionary PDEs, written in Python and SymPy.

Start with a flow $X=X_0+\sum_{k=1}^{N-1}\epsilon^k X_k+O(\epsilon^N)$.
On supported inputs, the package constructs a verified transformation and
its inverse, then uses the **same** transformation to transport an
admissible dispersionless symmetry. It verifies commutation and reports
jet polynomiality separately. Every result has the input's cutoff:
$O(\epsilon^N)$ means the coefficient of $\epsilon^N$ is the first unknown.

## Install and run

Python 3.11 or newer is required. From a checkout, preferably in a virtual
environment:

```sh
python -m pip install .
```

```python
from sympy import O, Symbol
from quasitriviality import JetSpace, Series, check_integrability, transport_symmetry

eps = Symbol("eps")
jets = JetSpace("u")
u, = jets.states
ux, uxx = jets.jet(0, 1), jets.jet(0, 2)
equation = Series.from_expr(u * ux + eps * uxx + O(eps**3), eps)

result = check_integrability(jets, equation)
T = result.transformation
print(T.reducing.as_expr(eps)[0])                # v = T(u)
print(T.inverse.as_expr(eps)[0])                 # inverse at the same cutoff
print(result.symmetry.symmetry.as_expr(eps)[0])  # arbitrary leading h(u) * ux
print(result.symmetry.polynomiality.status.value)  # polynomial

specific = transport_symmetry(T, u**2)          # leading u**2 * ux
```

For several fields, use e.g. `JetSpace(("r", "s"))` and give
`Series.from_expr` a tuple of expressions with the same explicit cutoff.
Build inputs as SymPy expressions using jets from the factory, not as
strings. Without an explicit leader the coupled solver uses arbitrary local
Tsarev speed values and their free pure derivatives. A supplied speed must
commute with the dispersionless flow, and its result covers only that choice.

Run `python -m unittest discover -v` from the checkout. The example runners
and coefficient archives are repository material, not importable modules
installed by the wheel. The longest benchmark runs are separate from the
default regression suite; follow the example instructions below to reproduce them.

## What the results mean

- Scalar input has leading term $f(u)u_x$, with nonzero $f'(u)$ on the
  chosen domain. The multicomponent path starts with a diagonal
  semi-Hamiltonian system in Riemann coordinates, with distinct speeds
  and nonzero self-derivatives. It does not find those coordinates for you.
- The scalar and coupled solvers use exact local quadratures in a finite
  Laurent–log jet class; negative first-jet powers and formal first-jet logs
  are allowed in the transformation. Choose a regular local domain where
  the required denominators and first jets do not vanish.
- A coupled first-order correction may involve a state-dependent affine
  Darboux equation. The default solver can return an exact local integral
  representation, not a numerical approximation or an elementary formula.
  See [Architecture](ARCHITECTURE.md#3-transformation-first) for the
  choice of basepoint and access to its defining equations.
- `UnresolvedTransformation` means construction or verification did not
  complete within the supported class. It is not, by itself, a proof of
  nonexistence. A verified transformation is still a valid result when a
  transported symmetry is certified nonpolynomial. The other polynomiality
  statuses are `polynomial`, `conditional`, and `undetermined`;
  `obstructed` denotes a certified nonzero condition.
- These are local formal calculations. Finite-order agreement does not
  establish all-order integrability, convergence, or global existence.

## Reproducible examples

| Example | Verified result |
| --- | --- |
| [Toda / CP¹](examples/toda/README.md) | Original and centered lattice maps through $O(\epsilon^7)$, including comparison with published genus-one and genus-two formulas. |
| [Toda symmetry transport](examples/toda/SYMMETRIES.md) | Three additional Toda flows through $O(\epsilon^7)$, checked against their difference equations. |
| [High-order Toda](examples/toda/HIGH_ORDER.md) | Centered maps through $O(\epsilon^{11})$, with timings and peak memory. This uses a benchmark-specific sparse implementation, not the generic solver. |
| [Corrected CP¹ loop equation](examples/toda/LOOP_EQUATION.md#from-the-ten-order-quasi-miura-map-through-genus-five) | Computed densities $F_1,\ldots,F_5$ satisfy the corrected equation through genus five. |
| [A₂ Boussinesq](examples/a2/README.md) and [parameter screen](examples/a2/SCREEN.md) | Maps and loop check through genus three; within the specified exact truncated family the screen leaves constant dispersion. |

## Documentation and references

- [Architecture](ARCHITECTURE.md) for conventions, mathematical assumptions,
  and verification.
- [Development notes](DEVELOPMENT.md), [original plan](PLAN.md),
  [constructive coupled design](CONSTRUCTIVE_PLAN.md), and
  [verification history](PROGRESS.md).
- [Known limitations](TODO.md), [release checklist](RELEASING.md), and
  [security policy](SECURITY.md).

The scalar method follows Liu–Zhang,
[On Quasitriviality and Integrability of a Class of Scalar Evolutionary PDEs](https://arxiv.org/abs/nlin/0510019),
§§4–5, with Joel Ekstrand's Mathematica implementation as an algorithmic
reference. The coupled method uses the generalized hodograph and exactness
results in [SemiHamiltonian in Lean](https://github.com/siqiliu-tsinghua/semi-hamiltonian/tree/6b53168).
Reference source files and papers are not bundled or needed at runtime.
The benchmark pages cite their mathematical formulas separately.

## License

[MIT](LICENSE).
