"""Cross-check transported Toda symmetries against the published Lax flows."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from examples.toda.symmetries import DATA, verify


class TodaSymmetriesTests(unittest.TestCase):
    def test_ordinary_and_logarithmic_flows_at_all_known_orders(self):
        with TemporaryDirectory() as folder:
            output = Path(folder)/"result.json"
            counts = verify(output=output)
            self.assertEqual({flow: rows[6] for flow, rows in counts.items()},
                             {"t21":(16,54), "t22":(56,122), "t11":(112,126)})
            self.assertEqual(output.read_bytes(), DATA.with_name(
                "toda-symmetries-epsilon6.json"
            ).read_bytes())

    def test_a_different_generator_is_not_silently_replaced(self):
        changed = json.loads(DATA.read_text(encoding="utf-8"))
        changed["generator"]["2"] = [[],[]]
        with TemporaryDirectory() as folder:
            bad = Path(folder)/"bad.json"
            bad.write_text(json.dumps(changed),encoding="utf-8")
            with self.assertRaisesRegex(ArithmeticError, "t21, order 2"):
                verify(bad,("t21",))


if __name__=="__main__":
    unittest.main()
