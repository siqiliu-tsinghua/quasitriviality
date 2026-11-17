"""Reconstruct A2 from its Boussinesq flow and check the DZ loop equation."""
import json
from fractions import Fraction as Q
from pathlib import Path
import unittest

from examples.a2 import compute as a2, loop_equation as loop, potentials


ARCHIVE = Path(a2.__file__).with_name("a2-epsilon6.json")


class A2BenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = a2.compute()
        cls.energies = loop.read_free_energies(ARCHIVE)

    def test_equation_only_sixth_order_archive(self):
        self.assertEqual(a2.serialize(self.data), json.loads(ARCHIVE.read_text(encoding="utf-8")))
        self.assertFalse(any(a2.input_flow(4)))
        self.assertFalse(any(a2.input_flow(6)))
        self.assertTrue(all(any(self.data["dress"][n]) for n in (2, 4, 6)))

        f1, f2, f3 = self.energies
        self.assertEqual(f1, (a2.log(0) + a2.log(1))/24)
        self.assertEqual(tuple(f.degree_order() for f in self.energies), (1, 4, 7))
        for i, sign in enumerate((1, -1)):
            self.assertEqual(f2.partial(a2.variable_index(i, 4)),
                             sign*a2.b/(192*a2.j(i, 1)**2))
            self.assertEqual(f3.partial(a2.variable_index(i, 7)),
                             a2.b**2/(2304*a2.j(i, 1)**3))
        self.assertEqual(tuple(len(f) for f in self.energies), (2, 32, 290))

    def test_independent_published_genera(self):
        self.assertEqual(potentials.verify_archive(ARCHIVE, algebra=a2), {1: 2, 2: 32})

    def test_generic_sympy_solver_at_first_dispersive_order(self):
        import sympy as sp
        from quasitriviality import JetSpace, Series, quasitrivialize
        from quasitriviality.jets import is_zero

        space = JetSpace(("r", "s"))
        r, s = space.states
        b = ((r-s)/12)**sp.Rational(1, 3)
        perturbation = space.dx(3*b**2, 3)/4
        leading = (3*b*space.jet(0, 1), -3*b*space.jet(1, 1))
        transform = quasitrivialize(space, Series((leading, (0, 0),
                                                    (perturbation, perturbation)), 3))
        indexed = (b,) + tuple(space.jet(i, k) for k in range(1, a2.N+1)
                               for i in range(2))
        expected = [sum(sp.Rational(c)*sp.prod(indexed[i]**e for i, e in powers)
                        for powers, c in vector.serial())
                    for vector in self.data["dress"][2]]
        self.assertTrue(all(is_zero(sp.factor(x-y)) for x, y in zip(
            transform.inverse.coefficient(2), expected
        )))

    def test_dz_loop_equation_symbolically_through_genus_three(self):
        self.assertEqual(loop.coefficients(self.energies), ((2, {}), (5, {}), (8, {})))

    def test_dz_loop_equation_detects_changed_genus_three_coefficient(self):
        f1, f2, f3 = self.energies
        changed = f3 + Q(1, 2304)*a2.b**2*a2.j(0, 7)/a2.j(0, 1)**3
        residuals = loop.coefficients((f1, f2, changed))
        self.assertTrue(all(not numerator for _, numerator in residuals[:2]))
        self.assertTrue(residuals[2][1])


if __name__ == "__main__":
    unittest.main()
