import unittest

from sympy import Function, O, Symbol

from quasitriviality import JetSpace, Series, check_integrability, quasitrivialize, transport_symmetry
from quasitriviality.coupled import CoupledSolver
from quasitriviality.jets import is_zero
from quasitriviality.polynomiality import PolynomialityStatus
from quasitriviality.series import lie_exponential
from quasitriviality.tsarev import TsarevFamily
from quasitriviality.transform import UnresolvedTransformation


class CoupledPipelineTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        self.r, self.s = self.space.states
        self.p, self.q = self.space.jet(0, 1), self.space.jet(1, 1)
        self.speeds = 2*self.r+self.s, self.r+2*self.s
        self.x0 = tuple(v*p for v,p in zip(self.speeds, (self.p, self.q)))

    def test_transform_first_coupled_generic_and_realized_family(self):
        space,r,s,p,q = self.space,self.r,self.s,self.p,self.q
        N = 3
        g = Series(((0,0),(p,0),(p**2,0)),N)
        leading = Series((self.x0,),N)
        equation = lie_exponential(space,g,leading)
        self.assertEqual(equation.coefficient(1), (-p*q,p*q))
        computed = check_integrability(space,equation)
        self.assertTrue(any(not is_zero(x) for x in computed.transformation.generator.coefficient(2)))
        self.assertTrue(computed.symmetry.universal_family)
        self.assertEqual(computed.symmetry.polynomiality.status,PolynomialityStatus.POLYNOMIAL)
        eps = Symbol("eps")
        self.assertEqual([x.getO() for x in computed.symmetry.symmetry.as_expr(eps)], [O(eps**N)]*2)
        self.assertEqual([x.getO() for x in computed.transformation.reducing.as_expr(eps)], [O(eps**N)]*2)

        a,b = Function("A")(r),Function("B")(s)
        P = (a+b)/(r-s)
        w = P.diff(r),P.diff(s)
        realized = transport_symmetry(computed.transformation,w)
        self.assertFalse(realized.universal_family)
        self.assertEqual(realized.polynomiality.status,PolynomialityStatus.POLYNOMIAL)

    def test_noncommuting_leader_and_unresolved_span(self):
        space,p,q = self.space,self.p,self.q
        eq = Series((self.x0,(-p*q,p*q)),2)
        with self.assertRaises(UnresolvedTransformation) as exc:
            quasitrivialize(space,eq,solver=CoupledSolver(state_basis=()))
        self.assertEqual(exc.exception.order,1)
        self.assertIn("empty search space",str(exc.exception))
        transform = quasitrivialize(space,eq)
        with self.assertRaises(ValueError):
            transport_symmetry(transform,(1,0))

    def test_named_symbolic_tsarev_family(self):
        x = quasitrivialize(self.space,Series((self.x0,),1))
        w = TsarevFamily(self.space,self.speeds,names=("fixture_one","fixture_two"))
        result = transport_symmetry(x,family=w)
        self.assertEqual(result.leading_speeds,w.values)
        self.assertTrue(result.universal_family)
