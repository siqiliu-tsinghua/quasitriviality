import unittest

import sympy as sp

from quasitriviality.jets import JetSpace, is_zero
from quasitriviality.scalar import (
    ScalarPrimitiveError,
    laurent_log_primitive,
    scalar_primitive,
)


class ScalarPrimitiveTests(unittest.TestCase):
    def setUp(self):
        self.j = JetSpace(("u",))
        self.u = self.j.states[0]
        self.p = self.j.jet(0, 1)
        self.q = self.j.jet(0, 2)
        self.ell = self.j.log(0)

    def assertSolved(self, f, source):
        result = scalar_primitive(self.j, f, source)
        self.assertTrue(is_zero(self.j.bracket((f * self.p,), (result,))[0] - source))
        return result

    def test_non_hopf_quadratic_source_requires_log(self):
        f = self.u**2
        got = self.assertSolved(f, self.p**2)
        self.assertTrue(is_zero(got - self.p * self.ell / (2 * self.u)))
        self.assertTrue(self.j.is_homogeneous(got, 1))

    def test_non_hopf_higher_jet_and_unknown_coefficient(self):
        f = sp.Function("f")(self.u)
        a = sp.Function("a")(self.u)
        s = a * self.j.jet(0, 3) + self.q**2 / self.p
        got = self.assertSolved(f, s)
        self.assertTrue(self.j.is_homogeneous(got, 2))

    def test_input_log_with_resonant_and_nonresonant_powers(self):
        source = self.p**2 * self.ell + self.p**3 * self.ell**2
        self.assertSolved(self.u**2, source)

    def test_power_log_antiderivative(self):
        expr = self.p**3 * self.ell**2 + self.ell / self.p + 7 / self.p**3
        primitive = laurent_log_primitive(expr, self.p, self.ell)
        self.assertTrue(is_zero(self.j.partial(primitive, 0, 1) - expr))

    def test_rejects_unsupported_forms(self):
        cases = [
            (sp.Integer(1), self.p**2),
            (self.u**2, sp.sin(self.p)),
            (self.u**2, self.ell**-1),
            (self.u**2, 1 / self.q),
            (self.u**2, self.p**2 / (1 + self.p)),
            (self.u**2, sp.Float("1.1") * self.p**2),
        ]
        for f, source in cases:
            with self.subTest(f=f, source=source):
                with self.assertRaises(ScalarPrimitiveError):
                    scalar_primitive(self.j, f, source)

    def test_expression_text_is_not_evaluated(self):
        with self.assertRaises(TypeError):
            scalar_primitive(self.j, "u**2", self.p**2)
        with self.assertRaises(TypeError):
            laurent_log_primitive("u_1**2", self.p, self.ell)

    def test_zero_source(self):
        self.assertEqual(self.assertSolved(self.u**2, sp.S.Zero), 0)


if __name__ == "__main__":
    unittest.main()
