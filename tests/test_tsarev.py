import unittest

import sympy as sp

from quasitriviality.jets import JetSpace, is_zero
from quasitriviality.tsarev import TsarevFamily


class TsarevFamilyTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        r, s = self.r, self.s = self.space.states
        self.speeds = (2 * r + s, r + 2 * s)
        self.family = TsarevFamily(self.space, self.speeds)

    def test_coupled_family_matches_explicit_local_potential(self):
        r, s = self.r, self.s
        first, second = self.family.values
        a, b = sp.Function("A"), sp.Function("B")
        potential = (a(r) + b(s)) / (r - s)
        concrete = {first: potential.diff(r), second: potential.diff(s)}
        self.assertTrue(is_zero(first.diff(s).xreplace(concrete).doit()
                                - second.diff(r).xreplace(concrete).doit()))
        expr = first.diff(r, s, s) + second.diff(r, r, s)
        actual = self.family.reduce(expr)
        self.assertTrue(is_zero((actual - expr).xreplace(concrete).doit()))
        self.assertLessEqual(actual.atoms(sp.Derivative), set(self.family.free_jets(3)))

    def test_reduce_rejects_expression_text(self):
        for text in ("r + s", b"r + s"):
            with self.subTest(text=text):
                with self.assertRaisesRegex(TypeError, "not expression text"):
                    self.family.reduce(text)
        self.assertEqual(self.family.reduce(self.r + self.s), self.r + self.s)

    def test_arbitrary_leading_characteristic_commutes(self):
        p, q = (self.space.jet(i, 1) for i in range(2))
        first, second = self.family.values
        bracket = self.space.bracket((self.speeds[0] * p, self.speeds[1] * q),
                                     (first * p, second * q))
        self.assertTrue(all(self.family.reduce(expr) == 0 for expr in bracket))

    def test_three_component_transverse_rules_are_compatible(self):
        space = JetSpace(("a", "b", "c"))
        a, b, c = space.states
        v = a + b + c
        family = TsarevFamily(space, (v, v + 1, v + 2))
        w0, w1, w2 = family.values
        via_b = (family.gamma(0, 1) * (w1 - w0)).diff(c)
        via_c = (family.gamma(0, 2) * (w2 - w0)).diff(b)
        self.assertEqual(family.reduce(via_b - via_c), 0)

    def test_reject_equal_or_non_semihamiltonian_speeds(self):
        r, s = self.space.states
        with self.assertRaisesRegex(ValueError, "distinct"):
            TsarevFamily(self.space, (r + s, r + s))
        space = JetSpace(("a", "b", "c"))
        a, b, c = space.states
        with self.assertRaisesRegex(ValueError, "semi-Hamiltonian"):
            TsarevFamily(space, (a + b * c, b, c))

    def test_reject_jet_dependent_speed(self):
        with self.assertRaisesRegex(ValueError, "positive physical jets"):
            TsarevFamily(self.space, (self.space.jet(0, 1), self.s))

    def test_leave_pure_derivatives_and_unrelated_functions(self):
        first, _ = self.family.values
        f = sp.Function("f")
        expr = first.diff(self.r, 3) + f(self.r, self.s).diff(self.s)
        self.assertEqual(self.family.reduce(expr), expr)


if __name__ == "__main__":
    unittest.main()
