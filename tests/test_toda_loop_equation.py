"""The Toda free energies solve the corrected DZ loop equation in all jets."""
from pathlib import Path
import unittest

from examples.toda import loop_equation as loop
from examples.toda.compute import j


DATA = Path(loop.__file__).with_name("toda-epsilon6.json")


class TodaLoopEquationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.densities = loop.read_free_energies(DATA)

    def test_exact_rational_identity_through_genus_three(self):
        self.assertEqual(loop.coefficients(self.densities),
                         ((2, {}), (5, {}), (8, {})))

    def test_each_printed_error_breaks_the_equation(self):
        printed_source = loop.coefficients(self.densities, printed_source=True)
        self.assertEqual(printed_source[0], (3, -2*loop.w*loop.a**2))
        printed_connection = loop.coefficients(self.densities, printed_connection=True)
        self.assertFalse(printed_connection[0][1])
        self.assertTrue(printed_connection[1][1])
        self.assertTrue(printed_connection[2][1])

    def test_changing_genus_three_density_is_detected(self):
        f1, f2, f3 = self.densities
        changed = f3 + loop.d**2 * j(0, 7) / j(0, 1)**3
        residuals = loop.coefficients((f1, f2, changed))
        self.assertTrue(all(not numerator for _, numerator in residuals[:2]))
        self.assertTrue(residuals[2][1])


if __name__ == "__main__":
    unittest.main()
