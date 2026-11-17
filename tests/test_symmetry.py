import unittest

from sympy import Function, O, Symbol

from quasitriviality import JetSpace, Series, check_integrability, transport_symmetry
from quasitriviality.jets import is_zero
from quasitriviality.polynomiality import PolynomialityStatus


class SymmetryTests(unittest.TestCase):
    def setUp(self):
        self.s = JetSpace("u")
        self.u, = self.s.states
        self.p, self.q, self.r = (self.s.jet(0, k) for k in (1, 2, 3))
        self.e = Symbol("eps")

    def test_kdv_arbitrary_speed_and_reuse(self):
        s,u,p,r,e = self.s,self.u,self.p,self.r,self.e
        result = check_integrability(s, Series.from_expr(u*p+e**2*r/12+O(e**3), e))
        self.assertEqual(result.symmetry.polynomiality.status, PolynomialityStatus.POLYNOMIAL)
        self.assertTrue(result.symmetry.universal_family)
        self.assertEqual(result.symmetry.symmetry.coefficient(0), (Function("h")(u)*p,))
        w = transport_symmetry(result.transformation, u)
        self.assertFalse(w.universal_family)
        for a,b in zip(w.symmetry.coefficients, result.transformation.equation.coefficients):
            self.assertTrue(is_zero(a[0]-b[0]))
        self.assertEqual(w.symmetry.as_expr(e)[0].getO(), O(e**3))

    def test_burgers_arbitrary_speed(self):
        s,u,p,q,e = self.s,self.u,self.p,self.q,self.e
        result = check_integrability(s, Series.from_expr(u*p+e*q+O(e**4), e))
        self.assertEqual(result.symmetry.polynomiality.status, PolynomialityStatus.POLYNOMIAL)
        self.assertEqual(result.transformation.reducing.as_expr(e)[0].getO(), O(e**4))
        self.assertEqual(result.symmetry.symmetry.as_expr(e)[0].getO(), O(e**4))

    def test_integrability_obstruction_does_not_hide_transform(self):
        # Liu–Zhang Example 5.2, f1=1, f2=1, which violates (5.8).
        s,u,p,q,e = self.s,self.u,self.p,self.q,self.e
        result = check_integrability(s, Series.from_expr(u*p+e*(q+p**2)+O(e**4), e))
        report = result.symmetry.polynomiality
        self.assertEqual(report.status, PolynomialityStatus.OBSTRUCTED)
        self.assertTrue(any(c.order == 3 and c.certified_nonzero for c in report.conditions))
        self.assertEqual(result.transformation.cutoff, result.symmetry.cutoff)

    def test_log_in_transformation_is_allowed(self):
        s,u,p,e = self.s,self.u,self.p,self.e
        result = check_integrability(s, Series.from_expr(u**2*p+e*p**2+O(e**2), e))
        self.assertTrue(result.transformation.reducing.coefficient(1)[0].has(s.log(0)))
        self.assertEqual(result.symmetry.polynomiality.status, PolynomialityStatus.POLYNOMIAL)
