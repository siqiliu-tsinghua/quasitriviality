import unittest

from sympy import O, Symbol

from quasitriviality.jets import JetSpace, is_zero
from quasitriviality.series import Series, lie_exponential, series_bracket, substitution


class SeriesTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("u",))
        self.u, = self.space.states
        self.p = self.space.jet(0, 1)

    def test_cutoff_and_internal_zero(self):
        eps = Symbol("eps")
        s = Series(((self.u,), (0,), (self.p,)), 4)
        self.assertEqual(s.coefficient(1), (0,))
        self.assertEqual(s.coefficient(3), (0,))
        self.assertEqual(s.as_expr(eps)[0].getO(), O(eps**4))
        with self.assertRaises(IndexError):
            s.coefficient(4)
        with self.assertRaises(ValueError):
            Series(((0,), (1,)), 1)

    def test_logarithmic_generator_sign_and_inverse(self):
        s, u, p = self.space, self.u, self.p
        g1 = -p * s.log(0)
        g = Series(((0,), (g1,)), 3)
        x0 = Series(((u*p,),), 3)
        x = lie_exponential(s, g, x0)
        self.assertTrue(is_zero(x.coefficient(1)[0] - p**2))
        back = lie_exponential(s, g.map(lambda y: -y), x)
        for a, b in zip(back.coefficients, x0.coefficients):
            self.assertTrue(is_zero(a[0] - b[0]))
        self.assertEqual(substitution(s, g).coefficient(1), (g1,))
        self.assertEqual(substitution(s, g, inverse=True).coefficient(1), (-g1,))

    def test_bracket_is_cutoff_limited(self):
        s = self.space
        x = Series(((self.u*self.p,),), 2)
        y = Series(((self.p,),), 3)
        self.assertEqual(series_bracket(s, x, y).cutoff, 2)
        with self.assertRaises(ValueError):
            series_bracket(s, x, y, cutoff=3)
