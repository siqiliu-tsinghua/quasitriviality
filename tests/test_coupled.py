import unittest

import sympy as sp

from quasitriviality.coupled import CoupledSolveUnresolved, CoupledSolver
from quasitriviality.darboux import DarbouxProblem, volterra_definitions
from quasitriviality.hodograph import HodographChart
from quasitriviality.jets import JetSpace, is_zero


class CoupledSolverTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        self.r, self.s = self.space.states
        self.p, self.q = (self.space.jet(i, 1) for i in range(2))
        self.speeds = (2 * self.r + self.s, self.r + 2 * self.s)
        self.x0 = (self.speeds[0] * self.p, self.speeds[1] * self.q)

    def solve(self, source, **options):
        result = CoupledSolver(**options)(self.space, self.speeds, source)
        actual = self.space.bracket(self.x0, result)
        self.assertTrue(all(is_zero(a - b) for a, b in zip(actual, source)))
        return result

    def manufactured(self, prescribed, **options):
        return self.solve(self.space.bracket(self.x0, prescribed), **options)

    def test_coupled_weight_one_with_forced_darboux_residual(self):
        f = (self.p * self.q, -self.p * self.q)
        self.assertNotEqual(self.solve(f, basepoint=(2, 0)), (0, 0))

    def test_coupled_weight_two_recovery(self):
        expected = (self.space.jet(0, 2), self.p * self.q)
        self.assertEqual(self.manufactured(expected), expected)

    def test_arbitrary_state_coefficient_and_formal_logs(self):
        a = sp.Function("a")(self.r, self.s)
        expected = (a * self.p**2 * self.space.log(0)**3,
                    self.q**3 / self.p * self.space.log(1))
        self.assertEqual(self.manufactured(expected), expected)

    def test_highest_own_and_transverse_jet_powers_with_laurent_dependence(self):
        chart = HodographChart(self.space, self.speeds)
        z, l = chart.jet, chart.log
        h = (z(0, 2)**2 / z(0, 1)**3 + z(1, 2) / z(0, 1)**2,
             z(1, 2) * l(0) / z(1, 1)**2)
        expected = chart.denormalize(h)
        self.assertTrue(all(is_zero(a-b) for a, b in zip(self.manufactured(expected), expected)))

    def test_weight_one_log_generators_with_no_candidate_input(self):
        a = sp.Function("a")(self.r)
        expected = (0, a * self.q * self.space.log(1))
        # The own-log normalization has zero log-constant, so here no
        # state-space Darboux correction is required.
        self.assertEqual(self.manufactured(expected), expected)

    def test_incompatible_source_does_not_return_a_witness(self):
        with self.assertRaises(CoupledSolveUnresolved):
            CoupledSolver()(self.space, self.speeds, (self.q**2, 0))

    def test_weight_and_exact_input_validation(self):
        with self.assertRaisesRegex(ValueError, "homogeneous"):
            CoupledSolver()(self.space, self.speeds, (self.p + self.p**2, 0))
        with self.assertRaisesRegex(ValueError, "floating-point"):
            CoupledSolver()(self.space, self.speeds, (0.5 * self.p**2, 0))

    def test_three_genuinely_coupled_components(self):
        space = JetSpace(("a", "b", "c"))
        a, _, _ = space.states
        p = tuple(space.jet(i, 1) for i in range(3))
        speeds = tuple(sum(space.states)+i for i in range(3))
        leading = tuple(v*q for v, q in zip(speeds, p))
        source = space.bracket(leading, (a*p[0], 0, 0))
        result = CoupledSolver(basepoint=(0, 0, 0))(space, speeds, source)
        self.assertTrue(all(is_zero(a-b) for a, b in zip(space.bracket(leading, result), source)))
        self.assertEqual(len(volterra_definitions(*result)), 1)

    def test_general_volterra_declares_exact_state_integral(self):
        chart = HodographChart(self.space, (3*self.r+self.s, self.r+4*self.s))
        problem = DarbouxProblem(chart, {(0, 1): 1, (1, 0): 2}, basepoint=(2, 0))
        h = problem.zero_axis_solution()
        self.assertEqual(h[0].subs(self.s, 0), 0)
        self.assertEqual(h[1].subs(self.r, 2), 0)
        self.assertTrue(all(q.has(sp.Integral) for q in problem.inhomogeneity()))
        self.assertEqual(problem.verify(h), h)


if __name__ == "__main__":
    unittest.main()
