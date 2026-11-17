"""Reproduce the complete power-law screen for the exact truncated A2 family."""
from dataclasses import dataclass

import sympy as sp

from .perturbation import FUNCTIONS, U, transform_residuals
from .perturbation_certificate import verify_certificate


@dataclass(frozen=True)
class PowerScreen:
    amplitudes: tuple
    powers: tuple
    first_conditions: tuple
    fourth_condition: object
    nonzero_first_amplitude_branch: dict
    zero_first_amplitude_branch: dict


def solve_power_laws():
    """Solve for independent real exponents, including identically zero terms."""
    conditions, obstruction = verify_certificate()
    amplitudes = sp.symbols("c1:4", real=True)
    powers = sp.symbols("m1:4", real=True)
    a, b, c = amplitudes
    m, n, ell = powers
    base = next(iter(conditions[2].free_symbols))
    functions = (conditions[0].expr, conditions[1].expr, conditions[2])
    substitution = dict(zip(functions, (a*base**m, b*base**n, c*base**ell)))
    # On U>0 each of the extracted powers is a unit; no amplitude
    # is used as a divisor. Evaluating at 1 merely reads the multiplier
    # of the single forced power of U after the symbolic substitution.
    first = tuple(expr.subs(substitution, simultaneous=True).doit().subs(base, 1)
                  for expr in conditions)
    if first != (a*m, b*n, c) or obstruction != -sp.Rational(2, 3):
        raise ArithmeticError("the differential certificate and power conditions disagree")
    return PowerScreen(amplitudes, powers, first, b**2,
                       {b: 0, c: 0, m: 0}, {a: 0, b: 0, c: 0})


def verify_transform_screen():
    """Check that the coupled triangular computation agrees with the certificate."""
    residuals = transform_residuals()
    A, B, C = FUNCTIONS
    if set(residuals.epsilon2_ideal) != {A.diff(U), A.diff(U, 2), B.diff(U), C, C.diff(U)}:
        raise ArithmeticError("unexpected transform residual ideal at order two")
    if tuple(residuals.epsilon4_ideal) != (B**2,):
        raise ArithmeticError("unexpected transform residual ideal at order four")
    return residuals


def main():
    result = solve_power_laws()
    verify_transform_screen()
    print("epsilon^2 necessary conditions:", result.first_conditions, "= (0, 0, 0)")
    print("epsilon^4 necessary condition:", result.fourth_condition, "= 0")
    print("c1 != 0:", result.nonzero_first_amplitude_branch)
    print("c1 = 0:", result.zero_first_amplitude_branch)
    print("All surviving flows have (a1,a2,a3) = (constant,0,0).")


if __name__ == "__main__":
    main()
