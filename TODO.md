# Known limitations

The planned finite-order implementation and benchmark validation are
complete. This file records scope limits, not commitments to further
development.

- The generic solver requires diagonal semi-Hamiltonian leading systems in
  Riemann coordinates. It does not compute those coordinates automatically.
- Solving locally with distinct speeds, nonzero self-derivatives and
  nonvanishing first jets does not give a global existence or convergence
  result. A failed symbolic construction is not a nonexistence theorem.
- General affine Darboux corrections may be returned as exact local
  Volterra definitions, without an elementary closed form.
- The high-order Toda and A₂ benchmarks use specialized sparse charts.
  Their achieved orders and timing do not predict general solver performance.
- The A₂ screen classifies the stated exact truncated perturbation; it does
  not classify families with independently adjustable higher coefficients.

For what was implemented and verified, see the [README](README.md) and
[verification history](PROGRESS.md).
