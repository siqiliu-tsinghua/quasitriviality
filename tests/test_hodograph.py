import unittest

import sympy as sp

from quasitriviality.hodograph import HodographChart
from quasitriviality.jets import JetSpace, is_zero


class HodographChartTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        r, s = self.space.states
        self.chart = HodographChart(self.space, (2 * r + s, r + 2 * s))

    def assert_expression_equal(self, left, right):
        self.assertTrue(is_zero(left - right), f"{left} != {right}")

    def test_coupled_triangular_second_order_and_round_trip(self):
        a, b = self.space.jet(0, 1), self.space.jet(1, 1)
        aa = self.space.jet(0, 2)
        gamma = self.chart.tsarev(0, 1)
        self.assert_expression_equal(self.chart.inverse_jet(0, 2), -aa / a**3 + gamma * b / a**2)
        self.assert_expression_equal(self.chart.jet_coefficient(0, 1, 1, 0), -gamma)
        self.assert_expression_equal(self.chart.jet_coefficient(0, 1, 2, 1), -gamma)
        f = aa / a**2 + self.space.jet(1, 3) * self.space.log(0) / b
        self.assert_expression_equal(self.chart.to_physical(self.chart.to_hodograph(f)), f)
        z = self.chart.jet(1, 3) / self.chart.jet(0, 1) ** 2 + self.chart.log(1)
        self.assert_expression_equal(self.chart.to_hodograph(self.chart.to_physical(z)), z)

    def test_conjugacy_including_log(self):
        a, b = self.space.jet(0, 1), self.space.jet(1, 1)
        r, s = self.space.states
        w = (3 * r**2 + 2 * r * s + s**2, r**2 + 2 * r * s + 3 * s**2)
        for res in self.chart.commuting_residual(w).values():
            self.assert_expression_equal(res, 0)
        flow = (w[0] * a, w[1] * b)
        f = self.space.jet(0, 2) / a + b * self.space.log(1)
        z = self.chart.to_hodograph(f)
        self.assert_expression_equal(self.chart.to_physical(self.chart.evolutionary(w, z)),
                                     self.space.evolutionary(flow, f))
        self.assert_expression_equal(self.chart.to_physical(self.chart.total(z)),
                                     self.space.dx(f))

    def test_arbitrary_tsarev_family(self):
        r, s = self.space.states
        potential = (sp.Function("A")(r) + sp.Function("B")(s)) / (r - s)
        for res in self.chart.commuting_residual((potential.diff(r), potential.diff(s))).values():
            self.assert_expression_equal(res, 0)

    def test_normalized_bracket(self):
        p, q = self.space.jet(0, 1), self.space.jet(1, 1)
        F = (q * self.space.log(0), self.space.jet(1, 2) / q)
        flow = (p * self.chart.speeds[0], q * self.chart.speeds[1])
        h = self.chart.normalize(F)
        for actual, expected in zip(self.chart.denormalize(h), F):
            self.assert_expression_equal(actual, expected)
        transported = self.chart.denormalize(self.chart.lie_action(self.chart.speeds, h))
        for actual, expected in zip(transported, self.space.bracket(flow, F)):
            self.assert_expression_equal(actual, expected)

    def test_scalar_joel_invariant(self):
        space = JetSpace(("u",))
        u, = space.states
        chart = HodographChart(space, (u**2,))
        p = space.jet(0, 1)
        self.assert_expression_equal(chart.inverse_jet(0, 2), -space.jet(0, 2) / p**3)
        z1, z2 = chart.jet(0, 1), chart.jet(0, 2)
        invariant = z2 - z1 / u
        self.assert_expression_equal(chart.lie_action(chart.speeds, (invariant,))[0], 0)
        self.assert_expression_equal(chart.lie_action(chart.speeds, (z1,))[0], -2 * u)
        self.assert_expression_equal(chart.to_physical(invariant / (2 * u)**2),
                                     -2 / ((2 * u)**3 * p) - space.jet(0, 2) / ((2 * u)**2 * p**3))


if __name__ == "__main__":
    unittest.main()
