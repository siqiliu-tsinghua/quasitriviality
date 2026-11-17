import unittest

import sympy as sp

from quasitriviality.coupled import CoupledSolveUnresolved, CoupledSolver
from quasitriviality.jets import JetSpace, is_zero


class CoupledSolverTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        self.r, self.s = self.space.states
        self.p, self.q = (self.space.jet(i, 1) for i in range(2))
        self.speeds = (2 * self.r + self.s, self.r + 2 * self.s)
        self.x0 = (self.speeds[0] * self.p, self.speeds[1] * self.q)

    def assertSolves(self, solver, source):
        result = solver(self.space, self.speeds, source)
        actual = self.space.bracket(self.x0, result)
        self.assertTrue(all(is_zero(a - b) for a, b in zip(actual, source)))
        return result

    def test_coupled_weight_one_recovery_in_hodograph_chart(self):
        source = (self.p * self.q, -self.p * self.q)
        result = self.assertSolves(CoupledSolver(), source)
        self.assertNotEqual(result, (0, 0))

    def test_coupled_weight_two_recovery_in_hodograph_chart(self):
        prescribed = (self.space.jet(0, 2), self.p * self.q)
        source = self.space.bracket(self.x0, prescribed)
        self.assertSolves(CoupledSolver(), source)

    def test_state_coefficient_basis(self):
        prescribed = (self.r * self.p, sp.S.Zero)
        source = self.space.bracket(self.x0, prescribed)
        self.assertSolves(CoupledSolver(state_basis=(1, self.r)), source)

    def test_log_candidate(self):
        prescribed = (self.p * self.space.log(0), sp.S.Zero)
        source = self.space.bracket(self.x0, prescribed)
        self.assertSolves(CoupledSolver(max_log_degree=1), source)

    def test_extra_vector_candidate_and_parameters(self):
        a = sp.Symbol("a")
        prescribed = (self.r * self.p, self.s * self.q)
        source = tuple(a * part for part in self.space.bracket(self.x0, prescribed))
        self.assertSolves(CoupledSolver(state_basis=(), extra_candidates=(prescribed,)), source)

    def test_inverse_exponents_in_finite_span(self):
        prescribed = (self.space.jet(0, 2) / self.p, sp.S.Zero)
        source = self.space.bracket(self.x0, prescribed)
        self.assertSolves(
            CoupledSolver(max_inverse_power=1, max_jet_order=2, use_hodograph=False),
            source,
        )

    def test_unresolved_is_not_an_obstruction(self):
        with self.assertRaisesRegex(CoupledSolveUnresolved, "selected"):
            CoupledSolver()(self.space, self.speeds, (self.q**2, 0))
        with self.assertRaisesRegex(CoupledSolveUnresolved, "bound"):
            CoupledSolver(max_candidates=1)(self.space, self.speeds, (self.p * self.q, 0))

    def test_weight_and_exact_input_validation(self):
        with self.assertRaisesRegex(ValueError, "homogeneous"):
            CoupledSolver()(self.space, self.speeds, (self.p + self.p**2, 0))
        with self.assertRaisesRegex(ValueError, "floating-point"):
            CoupledSolver()(self.space, self.speeds, (0.5 * self.p**2, 0))


if __name__ == "__main__":
    unittest.main()
