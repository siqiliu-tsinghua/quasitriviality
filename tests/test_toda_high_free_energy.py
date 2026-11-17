"""Regression check for tau extraction from the centered dressing map."""
import json
from pathlib import Path
import unittest

from examples.toda._high_order_backend import make_chart
from examples.toda import compute
from examples.toda.high_free_energy import read_sparse, recover


class HighTodaFreeEnergyTests(unittest.TestCase):
    def test_recovers_existing_densities_and_all_tau_identities(self):
        data = json.loads(Path(compute.__file__).with_name("toda-epsilon6.json").read_text())
        source = {
            "cutoff": 7, "basis": {"jet_order": compute.N},
            "coefficients": {order: {"dress": row} for order, row in data["dress"].items()},
        }
        chart = make_chart()
        densities = recover(chart, source)
        for genus in (2, 3):
            expected = read_sparse(chart, data["free_energy"][str(2*genus)][0],
                                   source_bound=compute.N)
            self.assertEqual(densities[genus], expected)
