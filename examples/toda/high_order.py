"""Profile exact centered Toda quasi-Miura maps through O(epsilon**11).

Run with `python -m examples.toda.high_order --verify --output toda-epsilon10.json`.
"""
import argparse
import json
import platform
import resource
import sys
from fractions import Fraction
from math import factorial
from pathlib import Path
from time import perf_counter

from . import compute
from ._high_order_backend import make_chart

NAMES = ("generator", "reduce", "dress")
Q = Fraction


def sum_vectors(chart, terms):
    rows = ({}, {})
    for vector, scale in terms:
        for result, component in zip(rows, vector, strict=True):
            for key, value in component.items():
                result[key] = result.get(key, 0) + scale*value
    return tuple(chart.L(row) for row in rows)


def benchmark(chart, cutoff=11):
    """Yield each timed coefficient from one run with warm lower orders."""
    if cutoff not in (3, 5, 7, 9, 11):
        raise ValueError("unsupported cutoff; expected 3, 5, 7, 9, or 11")
    adjoint, evolution, generators = {}, {}, {}
    start = perf_counter()
    for order in range(2, cutoff, 2):
        stage = perf_counter()
        target = chart.input_flow(order)
        # [epsilon^order] ad_G^k X0. For k=1, brackets already solved at
        # lower orders can be reused exactly.
        for k in range(2, order//2 + 1):
            adjoint[order, k] = sum_vectors(chart, (
                (chart.bracket(g, adjoint[order-j, k-1]), 1)
                for j, g in generators.items() if order-j >= 2*(k-1)))
        rhs = sum_vectors(chart, [(target, 1)] + [
            (adjoint[order, k], -Q(1, factorial(k))) for k in range(2, order//2 + 1)])
        generator = generators[order] = chart.solve(chart.mul(rhs, -1), order)
        adjoint[order, 1] = rhs
        gen_time, gen_peak = perf_counter()-stage, peak_mib()
        # [epsilon^order] E_G^k (r,s); exp(±E_G) gives inverse directions.
        evolution[order, 1] = generator
        for k in range(2, order//2 + 1):
            evolution[order, k] = sum_vectors(chart, (
                (chart.evol(g, evolution[order-j, k-1]), 1)
                for j, g in generators.items() if order-j >= 2*(k-1)))
        reduce, dress = (sum_vectors(chart, (
            (evolution[order, k], Q(sign**k, factorial(k)))
            for k in range(1, order//2 + 1))) for sign in (1, -1))
        finished = perf_counter()
        result = dict(zip(NAMES, (generator, reduce, dress), strict=True))
        yield order, result, {
            "order": order, "stage_seconds": finished-stage,
            "cumulative_seconds": finished-start, "peak_rss_MiB": peak_mib(),
            "generator_stage_seconds": gen_time, "generator_peak_MiB": gen_peak,
            "term_counts": {key: list(map(len, value)) for key, value in result.items()},
        }


def peak_mib():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw/(2**20 if sys.platform == "darwin" else 2**10)


def verify_archived_low_orders(chart, maps):
    """Compare each available low-order coefficient to the tracked exact archive."""
    archive = json.loads(Path(compute.__file__).with_name("toda-epsilon6.json").read_text())
    for order in (2, 4, 6):
        if order not in maps:
            continue
        for name in NAMES:
            before = archive[name][str(order)]
            for row in before:
                for exponents, _ in row:
                    for exponent in exponents:
                        if exponent[0] > 2*compute.N:
                            exponent[0] += 2*(chart.N-compute.N)
            chart.require([a.serial() for a in maps[order][name]] == before,
                          f"archived {name} mismatch at epsilon^{order}")


def verify_reducing(chart, maps):
    """Independently check E_X of the reducing map against X0 of its outputs."""
    flow = {order: chart.input_flow(order) for order in (0, *maps)}
    for order in flow:
        start = perf_counter()
        left = sum_vectors(chart, [(flow[order], 1)] + [
            (chart.evol(flow[k], maps[order-k]["reduce"]), 1) for k in range(0, order, 2)])
        right = []
        for component in range(2):
            value = chart.L()
            for k in range(0, order+1, 2):
                gap = chart.d if k == 0 else maps[k]["reduce"][0] - maps[k]["reduce"][1]
                derivative = (chart.j(component, 1) if k == order else
                              chart.dx(maps[order-k]["reduce"][component]))
                value += gap*derivative*Q(1 if component == 0 else -1, 4)
            right.append(value)
        chart.require(left == tuple(right), f"direct reducing PDE residual at epsilon^{order}")
        print(f"Direct reducing PDE, epsilon^{order}: exact zero ({perf_counter()-start:.3f}s)",
              flush=True)


def save(path, chart, cutoff, metrics, maps):
    """Write coefficient rows as [[exponent-index/value pairs], exact rational string]."""
    with path.open("w", encoding="utf-8") as out:
        out.write(json.dumps({
            "cutoff": cutoff, "python": platform.python_version(),
            "chart": "centered Toda, r=V+2exp(U/2), s=V-2exp(U/2)",
            "basis": {"gap_index": 0, "gap": "d=r-s", "jet_order": chart.N,
                      "jet_index": "1+2*(order-1)+component, component 0=r, 1=s",
                      "log_indices": [chart.log_index(0), chart.log_index(1)],
                      "log": "log(abs(first derivative))"},
            "benchmark": metrics,
        })[:-1])
        out.write(', "coefficients": {\n')
        for i, (order, row) in enumerate(maps.items()):
            if i:
                out.write(",\n")
            out.write(json.dumps(str(order)) + ": ")
            json.dump({name: [a.serial() for a in value] for name, value in row.items()},
                      out, separators=(",", ":"))
        out.write("\n}}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff", type=int, default=11, choices=(3, 5, 7, 9, 11))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", action="store_true", help="also check the direct reducing PDE")
    args = parser.parse_args(argv)
    chart, maps, metrics = make_chart(), {}, {}
    for order, row, metric in benchmark(chart, args.cutoff):
        maps[order], metrics[order] = row, metric
        print(json.dumps(metric), flush=True)
    verify_archived_low_orders(chart, maps)
    if args.output:
        save(args.output, chart, args.cutoff, metrics, maps)
    if args.verify:
        verify_reducing(chart, maps)


if __name__ == "__main__":
    main()
