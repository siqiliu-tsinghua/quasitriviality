"""Keep the private high-order chart and centered map convention aligned."""
import unittest
from contextlib import redirect_stdout
from io import StringIO

from examples.toda import compute
from examples.toda._high_order_backend import make_chart
from examples.toda.high_order import benchmark, verify_archived_low_orders, verify_reducing


class TodaHighOrderTests(unittest.TestCase):
    def test_extended_chart_matches_archived_maps(self):
        chart = make_chart()
        original_indices = (compute.N, compute.log_index(0), compute.V_INDEX)
        self.assertGreater(chart.N, compute.N)
        rows = dict((order, row) for order, row, _ in benchmark(chart, 7))
        verify_archived_low_orders(chart, rows)
        with redirect_stdout(StringIO()):
            verify_reducing(chart, rows)
        self.assertEqual((compute.N, compute.log_index(0), compute.V_INDEX),
                         original_indices)
