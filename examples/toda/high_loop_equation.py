"""Exact Dubrovin–Zhang CP1 loop checks for high-order Toda free energies.

The corrected period equation is the same as in loop_equation.py
and LOOP_EQUATION.md. This version uses the packed high-order chart
and clears the spectral denominator without any numerical evaluation.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import json
from math import comb
from pathlib import Path
from time import perf_counter

from ._high_order_backend import make_chart
from .high_free_energy import read_sparse


def read_archive(path):
    """Read independent equation-derived free energies in the packed basis."""
    archive = json.loads(Path(path).read_text(encoding="utf-8"))
    basis = archive["basis"]
    chart = make_chart()
    if not 0 < basis["jet_order"] <= chart.N:
        raise ValueError("unsupported archived jet basis")
    components = archive["free_energy"]
    orders = list(range(2, archive["cutoff"], 2))
    if set(components) != {str(order) for order in orders}:
        raise ValueError("expected consecutive even free-energy coefficients")
    rows = (components[str(order)][0] for order in orders)
    return chart, tuple(read_sparse(chart, row, source_bound=basis["jet_order"])
                        for row in rows)


def _accumulate(output, terms):
    for key, value in terms.items():
        if value := output.get(key, 0) + value:
            output[key] = value
        else:
            output.pop(key, None)


def _grad(chart, expr, order, jacobian):
    """Compute flat (V,U) jet derivatives in the canonical Laurent chart."""
    L = chart.L
    vd = [expr.partial(chart.V_INDEX)]
    ud = [expr.partial(chart.U_INDEX) + expr.partial(0)*chart.d/2]
    ud.extend(L() for _ in range(order))
    for n in range(1, order+1):
        sides = [expr.partial(chart.variable_index(alpha, n)) for alpha in (0, 1)]
        if n == 1:
            for alpha in (0, 1):
                log_partial = expr.partial(chart.log_index(alpha))
                if log_partial:
                    sides[alpha] += log_partial/chart.j(alpha, 1)
        vd.append(sides[0]+sides[1])
        diff = sides[0]-sides[1]
        if diff:
            for k in range(n+1):
                ud[k] += diff*jacobian[k, n]
    return {(alpha, n): value for alpha, row in enumerate((vd, ud))
            for n, value in enumerate(row) if value}


def _derivatives(chart, delta, initial, twice_power, order):
    """Numerators of derivatives of initial / delta**(twice_power/2)."""
    delta_x = chart.dx(delta)
    result = [initial]
    for i in range(order):
        previous = result[-1]
        result.append(chart.dx(previous)*delta - previous*delta_x*Q(twice_power+2*i, 2))
    return result


def _linear_term(terms, denominator, grad, kernel):
    if grad and kernel:
        _accumulate(terms[denominator], grad * kernel)


def residuals(chart, free_energies, *, printed_source=False, printed_connection=False):
    """Yield (genus, delta's power, exact cleared residual) for each F_g.

    For DeltaF=F1+epsilon**2 F2+..., genus g checks epsilon**(2*g-2).
    All tests are identities in arbitrary jets and the spectral parameter.
    Use the two switches only for diagnostic restoration of printed errors.
    """
    if not free_energies:
        raise ValueError("at least F1 is required")
    L, d, a, j = chart.L, chart.d, chart.V, chart.j
    w = d**2/16
    delta = a**2 - 4*w
    orders = [f.degree_order() for f in free_energies]
    order = max(orders)
    # Jacobian dr_n/dU_k = -ds_n/dU_k. V derivatives simply add rows.
    jacobian = {(k, n): comb(n, k)*(d if n == k else j(0, n-k)-j(1, n-k))/4
                for n in range(1, order+1) for k in range(n+1)}
    grads = [_grad(chart, f, n, jacobian) for f, n in zip(free_energies, orders)]

    units = (a, L(-2))
    resolvent = [_derivatives(chart, delta, n, 2, order) for n in units]
    period_order = max(order, max(orders[:-1], default=0)+1)
    periods = [_derivatives(chart, delta, n, 1, period_order) for n in units]
    p0half = [p / (-2) for p in periods[1]]

    vx = (j(0, 1)+j(1, 1))/2
    ux = 2*(j(0, 1)-j(1, 1))/d
    norm = a*a+4*w
    conn = (w*(norm*vx-4*w*a*ux), w*(norm*ux-4*a*vx))
    conn = [_derivatives(chart, delta, n, 6, max(orders[:-1], default=0)+1)
            for n in conn]

    @lru_cache(maxsize=None)
    def linear(alpha, k):
        lin = L(resolvent[alpha][k])
        for m in range(1, k+1):
            _accumulate(lin, p0half[m-1]*periods[alpha][k-m+1]*comb(k, m))
        return lin

    for g, grad in enumerate(grads):
        terms = defaultdict(dict)
        for (alpha, k), part in grad.items():
            _linear_term(terms, k+1, part, linear(alpha, k))
        if g == 0:
            if printed_source:
                _accumulate(terms[3], -w*norm)
            else:
                _accumulate(terms[2], w)
        else:
            quadratic = {}
            for left, part in grads[g-1].items():
                alpha, k = left
                scale = Q(-1, 2) if printed_connection else 1
                _linear_term(terms, k+4, part*scale, conn[alpha][k+1])
                for right, hess in _grad(chart, part, orders[g-1], jacobian).items():
                    # The symmetric Hessian gives identical contributions.
                    if left <= right:
                        quadratic[left, right] = hess*(2 if left < right else 1)
            # Pair terms related by exchanging both the jet and genus labels.
            for h in range(g//2+g%2):
                other = g-1-h
                for left, part in grads[h].items():
                    for right, partner in grads[other].items():
                        if h == other and left > right:
                            continue
                        pair = (left, right) if left <= right else (right, left)
                        weight = 1 if h == other and left == right else 2
                        _accumulate(quadratic.setdefault(pair, L()), part*partner*weight)
            for ((alpha, k), (beta, l)), coeff in quadratic.items():
                if not coeff:
                    continue
                # Multiplying smallest supports first reduces intermediate terms.
                pieces = sorted((coeff/4, periods[alpha][k+1], periods[beta][l+1]), key=len)
                _accumulate(terms[k+l+3], pieces[0]*pieces[1]*pieces[2])

        power = max(terms, default=0)
        # Horner clearing instead of multiplying every numerator by its
        # own expanded power, to give cancellations at each spectral level.
        numerator = L()
        for k in range(1, power+1):
            numerator = numerator * delta
            if k in terms:
                _accumulate(numerator, terms[k])
        yield g+1, power, numerator


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path,
                        default=Path(__file__).with_name("toda-free-energies-g1-g5.json"))
    parser.add_argument("--genus", type=int, help="check only through this genus")
    args = parser.parse_args(argv)
    chart, free_energies = read_archive(args.archive)
    if args.genus is not None:
        if not 1 <= args.genus <= len(free_energies):
            parser.error("--genus must be available in the archive")
        free_energies = free_energies[:args.genus]
    start = perf_counter()
    for g, p, numerator in residuals(chart, free_energies):
        print(f"genus {g}, epsilon^{2*g-2}: Delta^{p} times residual has "
              f"{len(numerator)} nonzero monomials ({perf_counter()-start:.2f}s elapsed)",
              flush=True)
        if numerator:
            raise ArithmeticError(f"the archived Toda F{g} violates the DZ loop equation")


if __name__ == "__main__":
    main()
