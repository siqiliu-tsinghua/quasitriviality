import unittest

import sympy as sp

from quasitriviality.darboux import DarbouxProblem, DarbouxUnresolved, solve_darboux, volterra_definitions
from quasitriviality.hodograph import HodographChart
from quasitriviality.jets import JetSpace, is_zero


class DarbouxTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        self.r, self.s = self.space.states

    def test_public_expression_helpers_reject_expression_text(self):
        chart = HodographChart(self.space, self.space.states)
        problem = DarbouxProblem(chart, {}, basepoint=(1, 0))
        for text in ("r + s", b"r + s"):
            with self.subTest(text=text):
                with self.assertRaisesRegex(TypeError, "not expression text"):
                    volterra_definitions(text)
                with self.assertRaisesRegex(TypeError, "not expression text"):
                    problem.contraction_condition(text, 2)
                with self.assertRaisesRegex(TypeError, "not expression text"):
                    problem.contraction_condition(sp.Rational(1, 64), text)
        self.assertEqual(volterra_definitions(self.r + self.s, 2), ())
        self.assertTrue(problem.contraction_condition(sp.Rational(1, 64), 2))

    def test_explicit_two_component_quadrature(self):
        r, s = self.space.states
        chart = HodographChart(self.space, (2*r+s, r+2*s))
        forcing = {(0, 1): r*s+1, (1, 0): r**2}
        answer = DarbouxProblem(chart, forcing, basepoint=(1, 0)).difference_solution()
        self.assertFalse(any(t.has(sp.Integral) for t in answer))
        self.assertFalse(volterra_definitions(*answer))
        DarbouxProblem(chart, forcing, basepoint=(1, 0)).verify(answer)

    def test_non_elementary_finite_quadrature(self):
        r, s = self.space.states
        chart = HodographChart(self.space, (2*r+s, r+2*s))
        forcing = {(0, 1): sp.Function("a")(r,s), (1, 0): sp.exp(r*s)}
        answer = DarbouxProblem(chart, forcing, basepoint=(1, 0)).difference_solution()
        self.assertTrue(any(t.has(sp.Integral) for t in answer))
        DarbouxProblem(chart, forcing, basepoint=(1, 0)).verify(answer)

    def test_certified_exact_volterra_backend(self):
        r, s = self.space.states
        chart = HodographChart(self.space, (2*r+3*s, r+2*s))
        forcing = {(0, 1): r*s+1, (1, 0): r**2}
        h = solve_darboux(chart, forcing, basepoint=(1, 0))
        self.assertEqual(len(volterra_definitions(*h)), 1)
        problem, = volterra_definitions(*h)
        self.assertEqual(h[0].subs(s, 0), 0)
        self.assertEqual(h[1].subs(r, 1), 0)
        self.assertTrue(is_zero(h[0].diff(r,s)-h[0].diff(s,r)))
        self.assertTrue(is_zero(h[1].diff(r,s)-h[1].diff(s,r)))
        self.assertTrue(any(t.has(sp.Integral) for t in problem.fixed_point_rhs()))
        self.assertTrue(problem.contraction_condition(sp.Rational(1,64), 2))

    def test_decoupled_three_component_coordinate_path(self):
        space = JetSpace(("r", "s", "t"))
        r,s,t = space.states
        forcing = {(0,1): 3*s**2+t, (0,2): s}
        chart = HodographChart(space, space.states)
        self.assertEqual(solve_darboux(chart, forcing, basepoint=(1,2,3)),
                         (s**3+s*t-14, 0, 0))

    def test_rejects_incompatible_forcing(self):
        space = JetSpace(("r", "s", "t"))
        r,s,t = space.states
        chart = HodographChart(space, space.states)
        with self.assertRaisesRegex(DarbouxUnresolved, "compatibility"):
            solve_darboux(chart, {(0,1): t}, basepoint=(1,2,3))

    def test_known_singular_basepoints_are_rejected(self):
        r, s = self.space.states
        straight = HodographChart(self.space, (r, s))
        with self.assertRaisesRegex(DarbouxUnresolved, "speed gap"):
            solve_darboux(straight, {(0, 1): 1}, basepoint=(1, 1))
        with self.assertRaisesRegex(DarbouxUnresolved, "forcing"):
            solve_darboux(straight, {(0, 1): 1/(r-1)}, basepoint=(1, 0))
        zero_self = HodographChart(self.space, (r**2, s))
        with self.assertRaisesRegex(DarbouxUnresolved, "self derivative"):
            solve_darboux(zero_self, {(0, 1):1}, basepoint=(0, 1))
        cross_singular = HodographChart(self.space, (2*r+sp.sqrt(s), r+2*s))
        with self.assertRaisesRegex(DarbouxUnresolved, "Gamma"):
            solve_darboux(cross_singular, {(0, 1):1}, basepoint=(1, 0))

    def test_high_mixed_derivatives_of_volterra_coefficients(self):
        from itertools import permutations
        r, s = self.space.states
        chart = HodographChart(self.space, (2*r+3*s, r+2*s))
        h = solve_darboux(chart, {(0,1): r*s+1, (1,0):r**2}, basepoint=(1,0))
        for value, own, transverse in zip(h, (r,s), (s,r)):
            expected = value.diff(transverse, own, own, own)
            for order in set(permutations((own,own,own,transverse))):
                actual = value
                for variable in order:
                    actual = actual.diff(variable)
                self.assertTrue(is_zero(actual - expected))

    def test_general_volterra_coefficient_in_a_full_transform(self):
        from quasitriviality import CoupledSolver, Series, quasitrivialize
        from quasitriviality.series import lie_exponential

        r, s = self.space.states
        p, q = self.space.jet(0, 1), self.space.jet(1, 1)
        speeds = (2*r + 3*s, r + 2*s)
        leading = Series(((speeds[0]*p, speeds[1]*q),), 3)
        g = Series(((0, 0), (r**2*p, r*s*q)), 3)
        equation = lie_exponential(self.space, g, leading)
        result = quasitrivialize(self.space, equation,
                                 solver=CoupledSolver(basepoint=(1,0)))
        self.assertEqual(result.cutoff, 3)
        expressions = [e for vector in result.generator.coefficients for e in vector]
        self.assertEqual(len(volterra_definitions(*expressions)), 1)
        self.assertEqual(lie_exponential(self.space, result.generator, leading), equation)

    def test_custom_solver_verified(self):
        r, s = self.space.states
        chart = HodographChart(self.space, (r,s))
        self.assertEqual(solve_darboux(chart, {(0,1):1}, solver=lambda chart,f: (s,0)), (s,0))
        with self.assertRaisesRegex(DarbouxUnresolved, "did not verify"):
            solve_darboux(chart, {(0,1):1}, solver=lambda chart,f: (r,0))


if __name__ == "__main__":
    unittest.main()
