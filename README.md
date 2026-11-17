# quasitriviality

A SymPy tool for constructing quasi-Miura transformations, transporting
arbitrary dispersionless symmetries, and extracting finite-order integrability
conditions for evolutionary PDEs.

Given an equation through `O(ε^N)`, the package returns a verified
quasi-Miura transformation and a symmetry for an arbitrary admissible leader
through the **same** `O(ε^N)`. The transformation is itself a research
output, not just an intermediate calculation. Its relation to free energies in
Dubrovin–Zhang theory motivates later cross-validation; a DZ or bihamiltonian
structure is not required.

## Quick start

From the repository root with SymPy installed, or after `pip install -e .`:

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
print(T.reducing.as_expr(eps)[0])          # v=T(u), through O(eps**3)
print(T.inverse.as_expr(eps)[0])           # u=T^{-1}(v) with same cutoff
print(result.symmetry.symmetry.as_expr(eps)[0])  # arbitrary h(u)
print(result.symmetry.polynomiality.status.value)  # polynomial

# Reuse exactly the same transformation with a selected leader.
specific = transport_symmetry(T, u**2)
```

For several components use e.g. `JetSpace(("r", "s"))` and pass one
expression per component to `Series.from_expr`; each needs the same explicit
big-O. Inputs must be constructed as SymPy expressions, with jet symbols
from the `JetSpace` factory. Default symmetry transport keeps arbitrary
local Tsarev speed values and their free pure derivatives symbolically.
An explicitly supplied speed must pass the leading commutator check and its
result applies only to that chosen leader.

Run `python -m unittest discover -v` from the checkout. The test suite
covers scalar KdV and Burgers, a scalar obstruction with a valid
transformation, a non-Hopf speed with a logarithmic transformation,
and a genuinely coupled two-component example through `O(ε³)`.

## Supported scope and limits

- Scalar solving uses exact Laurent–log quadrature for the nondegenerate
  leading speed `f(u)u_x`; coefficients depend on state and parameters.
- The coupled path takes diagonal semi-Hamiltonian systems in Riemann
  coordinates. The generalized hodograph action determines each unknown
  top transverse derivative and top own-jet coefficient. Integrate these
  exactly, descending jet orders, with no candidate generator basis or
  guessed exponents. The input is the equation, not an already known
  symmetry lift.
- A weight-one generator can need a state-dependent correction. The
  program checks the residual shape and affine Darboux compatibility.
  It uses exact quadrature in decoupled cases, and otherwise defines the
  correction by the exact local radial Volterra integral equation with
  zero axis data. For the two-component inverse-difference case,
  `DarbouxProblem.difference_solution()` also gives a finite quadrature.
  For the latter no finite approximation is substituted for the solution.
  Supply a regular coordinate basepoint with
  `quasitrivialize(jets, equation, solver=CoupledSolver(basepoint=(...)))`;
  by default the basepoint is symbolic. Choose a sufficiently small
  nonsingular cube centered there. `volterra_definitions(*coefficients)`
  retrieves the defining problems from returned SymPy expressions.
- Every answer is checked in physical jets. `UnresolvedTransformation`
  reports a failed step and its order and residual, not a conclusion
  of nonintegrability. Finite Laurent–log jet support is required;
  no global solution or elementary closed form for arbitrary
  smooth state functions is promised.
- Pick a local domain where speed differences, self-derivatives and any
  denominators in the answer are nonzero. Logs are formal first-jet
  logarithms, represented by `jets.log(i)`. The package does not select
  analytic branches or infer global domain assumptions.
- Polynomiality is assessed only through the cutoff. Conditional identities
  and undetermined expressions are kept distinct from certified obstructions.
  A result through finite order asserts neither all-order integrability nor
  convergence.

## Development documents

- [AGENTS.md](AGENTS.md): mathematical and engineering rules.
- [ARCHITECTURE.md](ARCHITECTURE.md): conventions, operations and result contract.
- [DEVELOPMENT.md](DEVELOPMENT.md): implementation and verification protocol.
- [PLAN.md](PLAN.md): milestone sequence and acceptance tests.
- [CONSTRUCTIVE_PLAN.md](CONSTRUCTIVE_PLAN.md): the revised coupled algorithm.
- [TODO.md](TODO.md): actionable backlog.
- [PROGRESS.md](PROGRESS.md): current status and verification history.

## References

1. [SemiHamiltonian in Lean](https://github.com/siqiliu-tsinghua/semi-hamiltonian/tree/6b53168):
   local copy at `../reference/semi-hamiltonian`, pinned for development at
   `6b53168`. The files on hodograph jets, normalized brackets and homological
   exactness specify the multicomponent operations and their hypotheses.
2. Liu–Zhang, [On Quasitriviality and Integrability of a Class of Scalar
   Evolutionary PDEs](https://arxiv.org/abs/nlin/0510019), §4–5:
   scalar construction, arbitrary leading symmetries, polynomiality criterion
   and examples for validation.
3. Joel Ekstrand’s Mathematica package, local reference
   `../reference/quasitriviality-joel/quasitriviality.m`.
   Use `reducingVector`, `reducingTransformation`, `transformVector` and
   `extractEquationsFromExpression` to compare algorithmic steps. The notebook
   contains examples and is not required or present.
