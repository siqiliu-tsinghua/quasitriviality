import unittest
from types import SimpleNamespace

import sympy as sp

from quasitriviality.jets import JetSpace
from quasitriviality.polynomiality import analyze_polynomiality, PolynomialityStatus


def series(*vectors):
    return SimpleNamespace(dimension=len(vectors[0]), cutoff=len(vectors),
                           coefficient=lambda k: vectors[k])


class TestPolynomiality(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("u", "v"))
        self.u, self.v = self.space.states
        self.p = self.space.jet(0, 1)
        self.q = self.space.jet(1, 1)
        self.r = self.space.jet(1, 2)
        self.log_p = self.space.log(0)

    def test_exact_cancellation_and_state_only_denominators(self):
        p, q, r = self.p, self.q, self.r
        log = self.log_p
        expr = ((p + q) * (p + log) - q * (p + log)) / p - log + r / (1 + self.u)
        report = analyze_polynomiality(self.space, series((expr, q)))
        self.assertEqual(report.status, PolynomialityStatus.POLYNOMIAL)
        self.assertEqual(report.conditions, ())

    def test_surviving_log_and_negative_power_with_tags(self):
        report = analyze_polynomiality(
            self.space, series((self.p, self.q), (0, 3*self.log_p + 2*self.r/self.p)),
        )
        self.assertEqual(report.cutoff, 2)
        self.assertEqual(report.status, PolynomialityStatus.OBSTRUCTED)
        self.assertEqual(
            {(c.order, c.component, c.monomial, c.coefficient) for c in report.conditions},
            {(1, 1, self.log_p, 3), (1, 1, self.r/self.p, 2)},
        )
        self.assertTrue(all(c.certified_nonzero for c in report.conditions))

    def test_parameter_condition_and_free_function_split(self):
        a = sp.Symbol("a")
        h = sp.Function("h")(self.u)
        h1 = sp.diff(h, self.u)
        report = analyze_polynomiality(
            self.space, series((a*h1/self.p, 0)), independent=(h, h1),
        )
        self.assertEqual(report.status, PolynomialityStatus.CONDITIONAL)
        self.assertEqual(report.conditions[0].coefficient, a)
        self.assertEqual(report.conditions[0].speed_monomial, h1)
        self.assertFalse(report.conditions[0].certified_nonzero)
        report = analyze_polynomiality(
            self.space, series(((2*h1 + a*h)/self.p, 0)), independent=(h, h1),
        )
        self.assertEqual(report.status, PolynomialityStatus.OBSTRUCTED)
        self.assertEqual(
            {c.speed_monomial: c.coefficient for c in report.conditions}, {h: a, h1: 2},
        )

    def test_exact_reducer_and_outside_supported_class(self):
        b = sp.Symbol("b")
        report = analyze_polynomiality(self.space, series((b/self.p, 1/(1 + self.q))),
                                       coefficient_reducer=lambda expr: expr.subs(b, 0))
        self.assertEqual(report.status, PolynomialityStatus.UNDETERMINED)
        self.assertEqual(report.conditions, ())
        self.assertEqual(report.unresolved[0].component, 1)

    def test_unevaluated_primitive_is_not_polynomial_certified(self):
        f = sp.Function("f")
        primitive = sp.Integral(f(self.u), self.u)
        report = analyze_polynomiality(self.space, series((primitive*self.p, 0)))
        self.assertEqual(report.status, PolynomialityStatus.UNDETERMINED)


if __name__ == "__main__":
    unittest.main()
