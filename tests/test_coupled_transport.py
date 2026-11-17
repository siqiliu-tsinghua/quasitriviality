import unittest

import sympy as sp

from quasitriviality.coupled_transport import (
    TransportUnresolved,
    first_primitive,
    jet_partial,
    jet_primitive,
    laurent_coefficients,
    solve_transport,
    transport_alpha,
    transport_operator,
)
from quasitriviality.hodograph import HodographChart
from quasitriviality.jets import JetSpace, is_zero


class CoupledTransportTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        self.r, self.s = self.space.states
        r, s = self.r, self.s
        self.chart = HodographChart(self.space, (2*r+s, r+2*s))

    def assert_expr(self, a, b):
        self.assertTrue(is_zero(a-b), f"{a} != {b}")

    def test_velocity_matches_chart_action(self):
        c = self.chart
        f = (self.r*c.jet(0, 2)*c.log(1) + c.jet(1, 2))/c.jet(1, 1)**2
        for i in range(2):
            expected = c.evolutionary(c.speeds, f) - c.speeds[i]*c.total(f) + transport_alpha(c,i)*f
            self.assert_expr(transport_operator(c, i, f), expected)

    def test_laurent_split_and_primitive(self):
        c = self.chart
        z, ell = c.jet(0, 1), c.log(0)
        f = ell**2 / z + 1 / (self.r*z**2 + self.s*z**2)
        terms = {tuple(sorted(exp.items(), key=lambda pair: str(pair[0]))): coeff
                 for exp, coeff in laurent_coefficients(c, f)}
        self.assertEqual(len(terms), 2)
        self.assert_expr(jet_partial(c, jet_primitive(c, f, 0, 1), 0, 1), f)
        self.assert_expr(first_primitive(c, 0, -2/z, 0), ell)

    def test_recover_order_two_and_logs(self):
        z, ell = self.chart.jet, self.chart.log
        for i, f, weight, order, alpha_factor in (
            (0, (ell(0)+self.r)/z(0,1) + (ell(1)-self.s)/z(1,1), 1, 1, 1),
            (0, z(0,2)**2/z(0,1)**3 + z(1,2)/(z(0,1)*z(1,1)), 1, 2, 1),
            (1, z(1,2)/z(1,1)**3 + z(0,2)/(z(0,1)**2*z(1,1)), 2, 2, -2),
        ):
            rhs = transport_operator(self.chart, i, f, alpha_factor)
            self.assert_expr(solve_transport(self.chart, i, rhs, weight, order, alpha_factor), f)

    def test_incompatible_source_is_unresolved(self):
        c = self.chart
        with self.assertRaises(TransportUnresolved):
            solve_transport(c, 0, 1/c.jet(1, 1)**2, 1, 1)
