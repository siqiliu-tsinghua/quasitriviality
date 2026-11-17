import unittest

import sympy as sp

from quasitriviality.jets import JetSpace, canonical, is_zero


class JetCalculusTests(unittest.TestCase):
    def test_scalar_total_derivative_and_formal_log(self):
        jet = JetSpace("u")
        (u,) = jet.states
        p, q, r = (jet.jet(0, s) for s in (1, 2, 3))
        ell = jet.log(0)
        h = sp.Function("h")
        self.assertEqual(jet.active_jets(ell), ((0, 1),))
        self.assertEqual(jet.dx(p + q), q + r)
        self.assertEqual(jet.dx(ell), q / p)
        self.assertTrue(is_zero(jet.dx(h(u) * p**-1)
                                - sp.diff(h(u), u) + h(u) * q / p**2))
        self.assertTrue(is_zero(jet.dx(ell, 2) - r / p + q**2 / p**2))

    def test_joel_bracket_sign(self):
        jet = JetSpace("u")
        (u,) = jet.states
        p = jet.jet(0, 1)
        f = sp.Function("f")(u)
        self.assertTrue(is_zero(jet.bracket((f * p,), (p * jet.log(0),))[0]
                                - sp.diff(f, u) * p**2))

    def test_evolutionary_commutes_with_total_derivative(self):
        jet = JetSpace(("r", "s"))
        u, v = jet.states
        p, q = jet.jet(0, 1), jet.jet(1, 1)
        Q = ((2 * u + v) * p, (u + 2 * v) * q)
        expr = (u + v) * jet.log(0) + jet.jet(1, 2) / p
        self.assertTrue(is_zero(jet.dx(jet.evolutionary(Q, expr)) -
                                jet.evolutionary(Q, jet.dx(expr))))
        self.assertEqual(jet.bracket(Q, (p, q)), (0, 0))
        self.assertEqual(jet.bracket((p, 0), Q), (-p * q, p * q))

    def test_weight_keeps_logs_fixed(self):
        jet = JetSpace("u")
        p, q = jet.jet(0, 1), jet.jet(0, 2)
        ell = jet.log(0)
        self.assertTrue(jet.is_homogeneous(q / p + p * ell**2, 1))
        self.assertTrue(jet.is_homogeneous(1 / p, -1))
        self.assertFalse(jet.is_homogeneous(q + p, 1))
        self.assertTrue(jet.is_homogeneous(sp.S.Zero, 9))

    def test_exact_cancellation_and_namespace(self):
        jet = JetSpace("u")
        (u,) = jet.states
        p = jet.jet(0, 1)
        ell = jet.log(0)
        self.assertEqual(canonical((ell * p + p**2) / p - ell), p)
        self.assertEqual(jet.jet(0, 0), u)
        with self.assertRaises(ValueError):
            JetSpace(("u", "u_1"))
        with self.assertRaises(ValueError):
            JetSpace(("u", "L_u"))
        with self.assertRaises(ValueError):
            jet.jet(0, -1)
        with self.assertRaises(IndexError):
            jet.log(1)


if __name__ == "__main__":
    unittest.main()
