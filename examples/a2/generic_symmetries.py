"""Exact arbitrary-Tsarev-leader screen for the integrable A2 power branch.

For a fixed order the generic symmetry is a finite linear form in the free
symbols w_i, D_i w_i, D_i**2 w_i, ... . The coefficients are sparse Laurent
expressions in b and the physical jets. No numerical choices of w are made.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import json
from pathlib import Path

from examples.a2.compute import L, b, dx, input_flow, j, jet_coeff, value_coeff
from examples.toda.compute import N, VARIABLE_COUNT, log_index, variable_index


DATA = Path(__file__).with_name("a2-epsilon6.json")


@lru_cache(None)
def speed_partial(label, direction):
    """Differentiate a free Tsarev jet and eliminate transverse derivatives."""
    i, n = label
    if direction == i:
        return {(i, n + 1): L(1)}
    a = value_coeff(i, direction, n)
    result = {(direction, 0): a, (i, 0): -a}
    result.update({(i, t + 1): jet_coeff(i, direction, n, t) for t in range(n)})
    return {k: v for k, v in result.items() if v}


def add(*forms):
    out = {}
    for form in forms:
        for k, v in form.items():
            out[k] = out.get(k, L()) + v
            if not out[k]:
                del out[k]
    return out


def mul(form, factor):
    return {k: w for k, v in form.items() if (w := v * factor)}


def total_derivative(form):
    result = {name: dx(v) for name, v in form.items()}
    for name, v in form.items():
        for direction in range(2):
            result = add(result, mul(speed_partial(name, direction), v * j(direction, 1)))
    return {k: v for k, v in result.items() if v}


def evolve_coefficients(vector, target):
    """Apply a physical characteristic (independent of w) to a generic target."""
    order = max(
        (v.degree_order() for component in target for v in component.values()), default=0
    )
    derivatives = [[v] for v in vector]
    for _ in range(order):
        for entry in derivatives:
            entry.append(dx(entry[-1]))
    delta = (vector[0] - vector[1]) / (36 * b**2)
    out = []
    for component in target:
        result = {}
        for name, v in component.items():
            coefficient = v.partial(0) * delta
            for i in range(2):
                for n in range(1, v.degree_order() + 1):
                    if a := v.partial(variable_index(i, n)):
                        coefficient += a * derivatives[i][n]
                if a := v.partial(log_index(i)):
                    coefficient += a * derivatives[i][1] / j(i, 1)
                result = add(result, mul(speed_partial(name, i), v * vector[i]))
            result = add(result, {name: coefficient})
        out.append(result)
    return tuple(out)


def evolve_generic(vector, target):
    """Apply a generic characteristic to a physical target independent of w."""
    order = max((component.degree_order() for component in target), default=0)
    derivatives = [[v] for v in vector]
    for _ in range(order):
        for entry in derivatives:
            entry.append(total_derivative(entry[-1]))
    delta = {
        k: v / (36 * b**2) for k, v in add(vector[0], mul(vector[1], -1)).items()
    }
    out = []
    for component in target:
        result = mul(delta, component.partial(0))
        for i in range(2):
            for n in range(1, component.degree_order() + 1):
                if a := component.partial(variable_index(i, n)):
                    result = add(result, mul(derivatives[i][n], a))
            if a := component.partial(log_index(i)):
                result = add(result, mul(derivatives[i][1], a / j(i, 1)))
        out.append(result)
    return tuple(out)


def bracket(physical, generic):
    return tuple(add(a, mul(b, -1)) for a, b in zip(
        evolve_coefficients(physical, generic), evolve_generic(generic, physical)
    ))


def vadd(*vectors):
    return tuple(add(*components) for components in zip(*vectors))


def vmul(vector, factor):
    return tuple(mul(c, factor) for c in vector)


def transport(generator):
    """Return coefficients of exp(ad_G) W at epsilon orders 0,2,4,6."""
    g2, g4, g6 = generator
    y0 = ({(0, 0): j(0, 1)}, {(1, 0): j(1, 1)})
    y2 = bracket(g2, y0)
    y22 = bracket(g2, y2)
    y40 = bracket(g4, y0)
    y4 = vadd(y40, vmul(y22, Q(1, 2)))
    y6 = vadd(
        bracket(g6, y0),
        vmul(vadd(bracket(g4, y2), bracket(g2, y40)), Q(1, 2)),
        vmul(bracket(g2, y22), Q(1, 6)),
    )
    return (y0, y2, y4, y6)


def decode(items):
    result = {}
    for key, coeff in items:
        exponents = [0] * VARIABLE_COUNT
        for i, power in key:
            exponents[i] = power
        result[tuple(exponents)] = Q(coeff)
    return L(result)


def obstruction_terms(vector):
    """Extract all forbidden physical-jet powers and formal first-jet logs."""
    return tuple({
        label: bad
        for label, terms in component.items()
        if (bad := L({
            key: c for key, c in terms.items()
            if any(e < 0 for e in key[1:1 + 2*N])
            or key[log_index(0)] or key[log_index(1)]
        }))
    } for component in vector)


def verify(archive=DATA, dispersion=Q(1, 4)):
    """Require exact generic polynomiality and commutation through O(eps**7).

    The archived generator belongs to dispersion coefficient 1/4. For any
    constant coefficient A, homogeneity in the expansion parameter multiplies
    its epsilon**(2*g) term by (4*A)**g. Negative or zero rational coefficient
    values therefore require no complex or singular square-root substitution.
    """
    factor = 4 * Q(dispersion)
    data = json.loads(Path(archive).read_text())
    generator = tuple(
        tuple(decode(c) * factor**(k // 2) for c in data["generator"][str(k)])
        for k in (2, 4, 6)
    )
    symmetries = transport(generator)
    x0 = input_flow(0)
    x2 = tuple(factor * coefficient for coefficient in input_flow(2))
    report = {}
    for idx, symmetry in enumerate(symmetries):
        order = 2 * idx
        obstructions = obstruction_terms(symmetry)
        if any(obstructions):
            raise ArithmeticError(f"epsilon^{order}: transported generic symmetry is not jet polynomial")
        commutation = bracket(x0, symmetry)
        if idx:
            commutation = vadd(commutation, bracket(x2, symmetries[idx - 1]))
        if any(commutation):
            raise ArithmeticError(f"epsilon^{order}: nonzero generic commutator")
        report[order] = {
            "free_speed_jet_order": max(
                (n for component in symmetry for _, n in component), default=0
            ),
            "monomial_counts": [sum(map(len, component.values())) for component in symmetry],
            "obstructions": [0, 0],
            "commutator_terms": [0, 0],
        }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DATA)
    parser.add_argument("--dispersion-coefficient", type=Q, default=Q(1, 4))
    args = parser.parse_args()
    report = verify(args.archive, args.dispersion_coefficient)
    for order, result in report.items():
        print(f"epsilon^{order}: {result}", flush=True)


if __name__ == "__main__":
    main()
