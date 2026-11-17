"""Regression checks for the Toda / CP1 benchmark."""
import json
from pathlib import Path
import unittest
from contextlib import redirect_stdout
from io import StringIO

from examples.toda.compare import compare, unpack
from examples.toda.compute import compute, serialize
from examples.toda.reference import d, p, q


class TodaBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result = compute()
        cls.data = serialize(result)

    def test_recomputed_sixth_order_data(self):
        path = Path(__file__).resolve().parents[1] / "examples" / "toda" / "toda-epsilon6.json"
        self.assertEqual(self.data, json.loads(path.read_text(encoding="utf-8")))
        self.assertEqual(self.data["metadata"]["cutoff"], 7)
        self.assertEqual(len(self.data["free_energy"]["6"][0]), 290)

    def test_published_genera(self):
        with redirect_stdout(StringIO()):
            self.assertEqual(compare(self.data), [(0, 0), (0, 0)])

    def test_generic_hodograph_solver_agrees_at_first_order(self):
        import sympy as sp

        from quasitriviality import JetSpace, Series, quasitrivialize
        from quasitriviality.jets import is_zero

        space = JetSpace(("r","s"))
        r,s = space.states
        h,V = (r-s)/4,(r+s)/2
        a,b = space.dx(h**2,3),h*space.dx(V,3)
        flow = Series(((h*space.jet(0,1), -h*space.jet(1,1)), (0,0),
                       ((a+b)/24,(a-b)/24)), 3)
        transform = quasitrivialize(space, flow)
        expected_substitution = {d:r-s, **{p[i-1]:space.jet(0,i) for i in range(1,5)},
                                       **{q[i-1]:space.jet(1,i) for i in range(1,5)}}
        expected = [unpack(row).xreplace(expected_substitution) for row in self.data["dress"]["2"]]
        self.assertTrue(all(is_zero(sp.factor(a-b)) for a,b in zip(
            transform.inverse.coefficient(2), expected
        )))


if __name__ == "__main__":
    unittest.main()
