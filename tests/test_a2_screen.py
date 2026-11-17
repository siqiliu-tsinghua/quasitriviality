"""Tests for the exact truncated A₂ power-law classification."""
from fractions import Fraction
import unittest

import sympy as sp

from examples.a2 import generic_symmetries, perturbation
from examples.a2.perturbation_screen import solve_power_laws, verify_transform_screen


class A2ScreenTests(unittest.TestCase):
    def test_exact_homological_residual_ideals(self):
        exact = verify_transform_screen()
        A, B, C = perturbation.FUNCTIONS
        b = perturbation.BROOT
        # First physical component, r_x^(-2) s_x (r_xx)^3.
        witness = exact.epsilon4_coefficients[0, ((1, -2), (2, 1), (3, 3))]
        self.assertEqual(witness, B**2/(540*b**2))
        for residual in exact.epsilon2_coefficients.values():
            self.assertEqual(sp.simplify(residual.subs({A: sp.Rational(1, 4), B: 0, C: 0}).doit()), 0)
        direct = perturbation.symmetry_residuals()
        self.assertEqual(set(direct[2]["ideal"]), set(exact.epsilon2_ideal))
        self.assertEqual(direct[4]["ideal"], (B**2,))

    def test_power_laws_include_zero_amplitudes(self):
        result = solve_power_laws()
        a, b, c = result.amplitudes
        m, n, _ = result.powers
        self.assertEqual(result.first_conditions, (a*m, b*n, c))
        self.assertEqual(result.fourth_condition, b**2)
        self.assertEqual(result.nonzero_first_amplitude_branch, {b: 0, c: 0, m: 0})
        self.assertEqual(result.zero_first_amplitude_branch, {a: 0, b: 0, c: 0})

    def test_negative_dispersion_and_zero_through_order_six(self):
        negative = generic_symmetries.verify(dispersion=Fraction(-3, 2))
        self.assertEqual(tuple(negative), (0, 2, 4, 6))
        for order in negative:
            self.assertEqual(negative[order]["obstructions"], [0, 0])
            self.assertEqual(negative[order]["commutator_terms"], [0, 0])
        self.assertEqual(negative[6]["monomial_counts"], [1069, 1069])
        zero = generic_symmetries.verify(dispersion=Fraction(0))
        self.assertTrue(all(zero[order]["monomial_counts"] == [0, 0]
                            for order in (2, 4, 6)))


if __name__ == "__main__":
    unittest.main()
