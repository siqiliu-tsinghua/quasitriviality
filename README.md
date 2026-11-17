# quasitriviality

A planned SymPy tool for constructing quasi-Miura transformations, transporting
arbitrary dispersionless symmetries, and extracting finite-order integrability
conditions for evolutionary PDEs.

Given an equation through `O(ε^N)`, the first milestone should return the
quasi-Miura transformation and the symmetry for an arbitrary admissible leading
term through the **same** `O(ε^N)`. The transformation is itself a research
output, not just an intermediate calculation. Its relation to free energies in
Dubrovin–Zhang theory motivates later comparison with examples from that theory;
the first milestone does not assume a DZ or bihamiltonian structure.

## Development documents

- [AGENTS.md](AGENTS.md): mathematical and engineering rules.
- [ARCHITECTURE.md](ARCHITECTURE.md): conventions, operations and result contract.
- [PLAN.md](PLAN.md): milestone sequence and acceptance tests.
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
