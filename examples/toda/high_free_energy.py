"""Recover the CP1 free energies from the centered high-order Toda archive."""
import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
from time import perf_counter

from ._high_order_backend import make_chart

Q = Fraction


def read_sparse(chart, row, *, source_bound):
    """Decode the archive's sparse variables into a private Toda chart."""
    terms = {}
    for powers, coefficient in row:
        key = chart.ZERO_KEY
        for index, exponent in powers:
            if index > 2*source_bound:
                index += 2*(chart.N-source_bound)
            key = chart.move_power(key, index, exponent)
        terms[key] = Q(coefficient)
    return chart.L(terms)


def recover(chart, archive, *, verbose=False):
    """Return F1..Fg, verifying the coordinate and t,t tau identities."""
    cutoff = archive["cutoff"]
    if cutoff not in (3, 5, 7, 9, 11):
        raise ValueError("expected an odd archive cutoff from 3 through 11")
    source_bound = archive["basis"]["jet_order"]
    if not 0 < source_bound <= chart.N:
        raise ValueError("archived jet basis is not supported by the chart")
    d = chart.d
    w = d**2/16
    leading = chart.input_flow(0)
    v, u, gaps = {}, {}, {}
    b, b2 = {0: Q(1)}, {0: Q(1)}
    densities = {}
    for g in range(1, cutoff//2+1):
        start = perf_counter()
        record = archive["coefficients"][str(2*g)]
        r, s = (read_sparse(chart, row, source_bound=source_bound) for row in record["dress"])
        gaps[g] = (r-s)/d
        v[g] = (r+s)/2
        # U-U0 = 2 log(1+gap). Compute in z=epsilon**2, using
        # (1+gap) d(log(1+gap))/dz = d(gap)/dz.
        u[g] = 2*gaps[g] - sum(
            (h*u[h]*gaps[g-h] for h in range(1, g)), chart.L())/g
        exp_coeff = w*(2*gaps[g] + sum(
            (gaps[h]*gaps[g-h] for h in range(1, g)), chart.L()))
        b[g] = -sum((b[g-h]*Q(1, 2**(2*h)*factorial(2*h+1))
                     for h in range(1, g+1)), Q(0))
        b2[g] = sum(b[h]*b[g-h] for h in range(g+1))
        u_tau = b2[g]*chart.x_derivative_leading_flat(0, 2*g)
        v_tau = b[g]*chart.x_derivative_leading_flat(1, 2*g)
        for h in range(1, g+1):
            u_tau += b2[g-h]*chart.dx(u[h], 2*(g-h))
            v_tau += b[g-h]*chart.dx(v[h], 2*(g-h))
        if g == 1:
            # The known genus-one density supplies the state-dependent log;
            # check it against the same archived map as every higher genus.
            f = (chart.log(0)+chart.log(1)-chart.U)/24
        else:
            f = chart.total_primitive(chart.total_primitive(u_tau))
        chart.require(chart.dx(f, 2) == u_tau, f"flat U tau identity at genus {g}")
        first, = chart.evol(leading, (f,))
        chart.require(chart.dx(first) == v_tau, f"flat V tau identity at genus {g}")
        second, = chart.evol(leading, (first,))
        chart.require(second == exp_coeff, f"t,t tau identity at genus {g}")
        densities[g] = f
        if verbose:
            print(f"genus {g}: {len(f)} terms, jet order {f.degree_order()}, "
                  f"three exact tau identities ({perf_counter()-start:.3f}s)", flush=True)
    return densities


def serialize(chart, archive, densities):
    basis = {**archive["basis"], "jet_order": chart.N,
             "log_indices": [chart.log_index(0), chart.log_index(1)],
             "V_index": chart.V_INDEX, "U_index": chart.U_INDEX}
    return {"cutoff": archive["cutoff"], "basis": basis,
            "gauge": "CP1 DZ flat coordinates: V_tau=B V, U_tau=B^2 U",
            "free_energy": {str(2*g): [f.serial()] for g, f in densities.items()}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args(argv)
    archive = json.loads(args.archive.read_text(encoding="utf-8"))
    chart = make_chart()
    densities = recover(chart, archive, verbose=True)
    if args.output:
        out = serialize(chart, archive, densities)
        args.output.write_text(json.dumps(out, separators=(",", ":"))+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
