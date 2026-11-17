"""Verify the Toda specialization of the DZ jet-space loop equation exactly.

Use the period formulation (3.10.95), not the misprinted specialization
(3.10.102) of arXiv:math/0108160. All operations are over rational Laurent
polynomials in the canonical jets, their gap d, and a = V - lambda.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path

from .compute import (
    L, U, U_INDEX, V, V_INDEX, ZERO_KEY, d, dx, j, log, log_index, variable_index,
)


# V's sparse index here denotes a=V-lambda. Their total derivatives agree.
a = V
w = d**2 / 16
delta = a**2 - 4*w
delta_x = dx(delta)
unit_numerators = (a, L(-2))


def read_sparse(rows):
    terms = {}
    for powers, coefficient in rows:
        key = list(ZERO_KEY)
        for index, exponent in powers:
            key[index] = exponent
        terms[tuple(key)] = Fraction(coefficient)
    return L(terms)


def read_free_energies(archive):
    data = json.loads(Path(archive).read_text(encoding="utf-8"))
    f1 = (log(0) + log(1) - U) / 24
    return (f1,) + tuple(read_sparse(data["free_energy"][str(2*g)][0]) for g in (2, 3))


@lru_cache(maxsize=None)
def flat_u_jacobian(k, n):
    """Differentiation of r^(n) with respect to the flat U^(k)."""
    return comb(n, k) * (d if n == k else j(0, n-k)-j(1, n-k)) / 4


def flat_partial(expr, component, k):
    if component == 0:
        if k == 0:
            return expr.partial(V_INDEX)
        result = expr.partial(variable_index(0, k)) + expr.partial(variable_index(1, k))
        if k == 1:
            result += expr.partial(log_index(0)) / j(0, 1) + expr.partial(log_index(1)) / j(1, 1)
        return result
    result = expr.partial(U_INDEX) + expr.partial(0)*d/2 if k == 0 else L()
    for n in range(max(k, 1), expr.degree_order()+1):
        jac = expr.partial(variable_index(0, n)) - expr.partial(variable_index(1, n))
        if n == 1:
            jac += expr.partial(log_index(0))/j(0, 1) - expr.partial(log_index(1))/j(1, 1)
        if jac:
            result += jac * flat_u_jacobian(k, n)
    return result


def derivative_numerators(numerator, twice_power, order):
    """Build numerators of D_x^k(numerator/delta**(twice_power/2))."""
    rows = [numerator]
    for k in range(order):
        rows.append(dx(rows[-1])*delta - Fraction(twice_power+2*k, 2)*rows[-1]*delta_x)
    return rows


def add_product(terms, denominator_exponent, *factors):
    if any(not item for item in factors):
        return
    term = L(1)
    for factor in sorted(factors, key=len):
        term = term * factor
    if term:
        acc = terms[denominator_exponent]
        for key, coeff in term.items():
            acc[key] = acc.get(key, 0) + coeff


def coefficients(free_energies, *, printed_source=False, printed_connection=False):
    """Return cleared numerator of LHS minus RHS for epsilon^0,2,4.

    Optional switches reintroduce separate errors of (3.10.102); they are
    diagnostic mutations, not alternative loop equations.
    """
    if len(free_energies) != 3:
        raise ValueError("expected F1, F2, F3")
    orders = [f.degree_order() for f in free_energies]
    order = max(orders)
    indices = [[(alpha,k) for alpha in range(2) for k in range(n+1)] for n in orders]
    grads = [{(alpha,k): flat_partial(f, alpha, k) for alpha,k in idx}
             for f,idx in zip(free_energies,indices)]

    resolvent = [derivative_numerators(n, 2, order) for n in unit_numerators]
    periods = [derivative_numerators(n, 1, order+1) for n in unit_numerators]
    p0half = [p / (-2) for p in periods[1]]

    vx = (j(0,1)+j(1,1))/2
    ux = 2*(j(0,1)-j(1,1))/d
    norm = a*a+4*w
    # Numerator of -(G11/2) (grad p_{1,lambda}) · (grad p_{1,lambda}) · v_x
    # in the Frobenius product, used when subtracting the corrected RHS.
    connection = (w*(norm*vx-4*w*a*ux), w*(norm*ux-4*a*vx))
    connection = [derivative_numerators(n, 6, order+1) for n in connection]
    residuals = []
    for g in range(3):
        terms = defaultdict(dict)
        for (alpha,k), grad in grads[g].items():
            add_product(terms, k+1, grad, resolvent[alpha][k])
            if k and grad:
                product = sum((comb(k,m) * p0half[m-1] * periods[alpha][k-m+1]
                               for m in range(1,k+1)), L())
                add_product(terms, k+1, grad, product)

        if g == 0:
            if printed_source:
                add_product(terms, 3, -w*norm)
            else:
                add_product(terms, 2, w)
        else:
            for (alpha,k), grad in grads[g-1].items():
                c = Fraction(-1, 2) if printed_connection else 1
                add_product(terms, 4+k, c*grad, connection[alpha][k+1])
                if grad:
                    for beta,l in indices[g-1]:
                        hess = flat_partial(grad, beta, l)
                        add_product(terms, k+l+3, hess/4, periods[alpha][k+1],
                                    periods[beta][l+1])
            for left_genus in range(g):
                right_genus = g-1-left_genus
                for (alpha,k), grad in grads[left_genus].items():
                    if grad:
                        for (beta,l), other in grads[right_genus].items():
                            add_product(terms, k+l+3, grad/4, other,
                                        periods[alpha][k+1], periods[beta][l+1])

        power = max(terms, default=0)
        cleared = sum((L(terms[p])*delta**(power-p) for p in sorted(terms)), L())
        residuals.append((power, cleared))
    return tuple(residuals)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path, default=Path(__file__).with_name(
        "toda-epsilon6.json"))
    args = parser.parse_args(argv)
    residuals = coefficients(read_free_energies(args.archive))
    for g, (p, residual) in enumerate(residuals, start=1):
        print(f"genus {g}, epsilon^{2*g-2}: D^{p} times residual has "
              f"{len(residual)} nonzero monomials", flush=True)
    if any(residual for _, residual in residuals):
        raise ArithmeticError("the archived Toda free energies violate the DZ loop equation")


if __name__ == "__main__":
    main()
