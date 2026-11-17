# Progress

## 2026-11-17: references and design

- Initialized the local project repository on `main`.
- Reviewed the pinned `semi-hamiltonian` Lean reference, Liu–Zhang
  §4–5, and Joel Ekstrand’s Mathematica source stored as an adjacent
  reference. The Mathematica notebook was not supplied. The package
  was read, not executed; the Lean reference was not rebuilt.
- Established that the main workflow is **transform → transport arbitrary
  leading symmetry → test jet polynomiality**. The transformation is an
  independent output and a future point of comparison with DZ free energies.
- Compared the two scalar constructions: the original paper applies its
  integral in Hopf coordinates; Joel keeps the general nondegenerate speed.
  Lean’s generalized hodograph jets and normalized action extend the
  latter approach structurally, not by componentwise scalar integration.
- Documented the common cutoff, sign convention, staged implementation
  and acceptance tests. The general coupled symbolic solve and its
  coefficient domain remain to be established.

**Status:** design only. No Python package or tests have been written or run.
Next implementation step: Phase 1A jet calculus and cutoff tests.
