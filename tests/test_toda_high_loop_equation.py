"""Regression of the packed period equation and its two errata."""
from pathlib import Path
import unittest

from examples.toda import high_loop_equation as loop
from examples.toda.high_loop_canonical import CanonicalLoop


DATA = Path(loop.__file__).with_name("toda-free-energies-g1-g5.json")


class HighTodaLoopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chart, cls.densities = loop.read_archive(DATA)

    def test_exact_first_three_genera(self):
        rows = list(loop.residuals(self.chart, self.densities[:3]))
        self.assertEqual(rows, [(1, 2, {}), (2, 5, {}), (3, 8, {})])
        self.assertEqual(tuple(map(len, self.densities)), (3, 30, 290, 1824, 9199))

    def test_independent_laurent_residual_agrees(self):
        rows = list(CanonicalLoop(self.chart).residuals(self.densities[:3]))
        self.assertEqual(rows, [(1, {}), (2, {}), (3, {})])

    def test_breaking_either_erratum_is_detected(self):
        chart, first = self.chart, self.densities[:2]
        rows = tuple(loop.residuals(chart, first, printed_source=True))
        self.assertEqual(rows[0], (1, 3, -chart.d**2 * chart.V**2/8))
        self.assertFalse(rows[1][2])
        rows = tuple(loop.residuals(chart, first, printed_connection=True))
        self.assertFalse(rows[0][2])
        self.assertTrue(rows[1][2])
