# Quasitriviality: repository guidance

This project implements a symbolic, finite-order quasi-Miura construction and
formal integrability test in Python/SymPy. Read [ARCHITECTURE.md](ARCHITECTURE.md)
for conventions, [PLAN.md](PLAN.md) for acceptance criteria, and [TODO.md](TODO.md)
and [PROGRESS.md](PROGRESS.md) before choosing the next task.

## Mathematical requirements

- **Transform first.** Compute and return a quasi-Miura transformation
  independently of the test. Then use that same transformation to transport an
  arbitrary admissible dispersionless symmetry. Testing jet polynomiality is a
  third, separate step. A transformation is a valuable output even when a
  transported symmetry is not polynomial.
- Preserve the leading speeds in the input coordinates. Use Joel’s adapted
  hodograph integration for one component; use the Lean reference’s generalized
  hodograph jets and normalized action for multiple components. Construct by
  integration; do not enumerate candidates or impose a guessed jet basis. Do not assume
  coupled speeds can be made into independent Hopf coordinates.
- Work formally to an explicit cutoff. Input `O(ε^N)` means the coefficient of
  `ε^N` and everything above it are unknown. Do not silently set them to zero
  and report a result at a higher order.
- Use exact symbolic arithmetic. Keep state-dependent coefficients distinct
  from positive-order physical jets and auxiliary hodograph jets. Preserve
  arbitrary functions and their derivative relations. Do not test universality
  by sampling a few functions or points.
- Verify results in physical jet coordinates using the evolutionary bracket.
  Keep transformation existence, symmetry commutation, and jet polynomiality
  distinct. Unsupported symbolic cases are inconclusive, not obstructions.
- Scope the calculus to local domains with the required nonvanishing speed
  differences and self-derivatives. First-jet inverses and formal logs are
  allowed in the intermediate algebra. Jet polynomiality is tested after
  returning to physical jets and simplifying exactly.
- Do not claim convergence, full-order integrability from finite checks, or
  unconditional multicomponent quasitriviality.

## Development

- Keep each stage independently testable: jet calculus, truncated series,
  hodograph transport, normalization, symmetry transport, obstruction reporting.
- Use the three references in [README.md](README.md) for cross-checks, not as
  runtime dependencies. Do not edit the reference copies or rely on Mathematica
  being installed. Prefer a fresh implementation to copying reference code.
- Keep examples exact and record their assumptions and order. Update
  `PROGRESS.md` with what actually ran and `TODO.md` with unresolved work.
