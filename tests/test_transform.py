import unittest

from sympy import O, Float, Symbol

from quasitriviality.jets import JetSpace, is_zero
from quasitriviality.series import Series, lie_exponential
from quasitriviality.transform import UnresolvedTransformation, quasitrivialize


class TransformTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace("u")
        self.u, = self.space.states
        self.p, self.q, self.r = (self.space.jet(0, s) for s in (1, 2, 3))
        self.eps = Symbol("eps")

    def test_explicit_cutoff_is_retained(self):
        s, u, p, e = self.space, self.u, self.p, self.eps
        X = Series.from_expr(u*p + O(e), e)
        q = quasitrivialize(s, X)
        self.assertEqual(q.cutoff, 1)
        self.assertEqual(q.reducing.as_expr(e), (u+O(e),))
        self.assertEqual(q.generator.coefficient(0), (0,))
        with self.assertRaises(ValueError):
            Series.from_expr(u*p, e)

    def test_generic_non_hopf_log_feature(self):
        s, u, p, e = self.space, self.u, self.p, self.eps
        X = Series.from_expr(u**2*p + e*p**2 + O(e**2), e)
        q = quasitrivialize(s, X)
        self.assertEqual(q.speeds, (u**2,))
        self.assertTrue(is_zero(q.generator.coefficient(1)[0]+p*s.log(0)/(2*u)))
        self.assertEqual(q.reducing.as_expr(e)[0].getO(), O(e**2))

    def test_kdv_coefficient(self):
        s, u, p, q, r, e = self.space, self.u, self.p, self.q, self.r, self.eps
        X = Series.from_expr(u*p+e**2*r/12+O(e**3), e)
        t = quasitrivialize(s, X)
        self.assertEqual(t.reducing.coefficient(1), (0,))
        self.assertTrue(is_zero(t.reducing.coefficient(2)[0]-(q**2/p**2-r/p)/24))

    def test_burgers_reconstructs_higher_order(self):
        s, u, p, q, e = self.space, self.u, self.p, self.q, self.eps
        X = Series.from_expr(u*p+e*q+O(e**4), e)
        t = quasitrivialize(s, X)
        computed = lie_exponential(s, t.generator, Series(((u*p,),), 4))
        for a,b in zip(computed.coefficients,X.coefficients):
            self.assertTrue(is_zero(a[0]-b[0]))
        self.assertTrue(is_zero(t.reducing.coefficient(1)[0]+q/p))

    def test_bad_witness_and_bad_grading(self):
        s, u, p, e = self.space, self.u, self.p, self.eps
        X = Series.from_expr(u*p+e*p**2+O(e**2), e)
        with self.assertRaises(UnresolvedTransformation) as failure:
            quasitrivialize(s, X, solver=lambda *_: (0,))
        self.assertEqual(failure.exception.order, 1)
        with self.assertRaises(ValueError):
            quasitrivialize(s, Series(((u*p,), (Float("0.1")*p**2,)), 2))
        with self.assertRaises(ValueError):
            quasitrivialize(s, Series(((p,),), 1))
