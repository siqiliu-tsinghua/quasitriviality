import unittest

import sympy as sp

from quasitriviality.cocycle import (
    CocyclePrimitiveError,
    DarbouxPrimitiveRequired,
    primitive_from_lift,
)
from quasitriviality.coupled import CoupledSolver
from quasitriviality.hodograph import HodographChart
from quasitriviality.jets import JetSpace, is_zero
from quasitriviality.tsarev import TsarevFamily


class CocyclePrimitiveTests(unittest.TestCase):
    def setUp(self):
        self.space = JetSpace(("r", "s"))
        r, s = self.space.states
        self.speeds = (2*r + s, r + 2*s)
        self.chart = HodographChart(self.space, self.speeds)
        self.family = TsarevFamily(self.space, self.speeds)
        self.p, self.q = (self.space.jet(i, 1) for i in range(2))
        self.X = tuple(v*p for v, p in zip(self.speeds, (self.p, self.q)))
        self.W = tuple(w*p for w, p in zip(self.family.values, (self.p, self.q)))

    def source_and_lift(self, F):
        return self.space.bracket(self.X, F), self.space.bracket(self.W, F)

    def test_weight_above_one_recovers_log_potential(self):
        F = (self.p**2 * (self.space.log(0)**2 + self.space.log(1)), sp.S.Zero)
        source, lift = self.source_and_lift(F)
        found = primitive_from_lift(self.chart, self.family, source, lift, 3)
        self.assertTrue(all(is_zero(a - b) for a, b in zip(found, F)))
        direct = CoupledSolver()(self.space, self.speeds, source)
        self.assertTrue(all(is_zero(a - b) for a, b in zip(direct, found)))

    def test_weight_one_uses_darboux_state_primitive(self):
        source, lift = self.source_and_lift((self.p, sp.S.Zero))
        with self.assertRaises(DarbouxPrimitiveRequired) as needed:
            primitive_from_lift(self.chart, self.family, source, lift, 2)
        self.assertTrue(needed.exception.equations)

        def solver(chart, rhs):
            self.assertEqual(rhs, needed.exception.equations)
            self.assertIs(chart, self.chart)
            return (sp.S.One, sp.S.Zero)

        found = primitive_from_lift(self.chart, self.family, source, lift, 2,
                                    darboux_solver=solver)
        self.assertEqual(found, (self.p, 0))

    def test_weight_one_vertical_log(self):
        F = (self.p * self.space.log(0), sp.S.Zero)
        source, lift = self.source_and_lift(F)
        found = primitive_from_lift(self.chart, self.family, source, lift, 2)
        self.assertTrue(all(is_zero(a - b) for a, b in zip(found, F)))

    def test_reject_unverified_generic_lift(self):
        source, _ = self.source_and_lift((self.p, sp.S.Zero))
        with self.assertRaisesRegex(CocyclePrimitiveError, "failed its physical bracket"):
            primitive_from_lift(self.chart, self.family, source, (0, 0), 2)


if __name__ == "__main__":
    unittest.main()
