"""Recover a homological primitive from a verified generic Tsarev lift.

This optional checker follows the cocycle construction in the Lean reference.
A generic lift is explicit input; this is not an inverse on arbitrary sources.
"""

from __future__ import annotations

from itertools import combinations
from operator import index

import sympy as sp

from .jets import canonical, is_zero
from .coupled_transport import TransportUnresolved, jet_partial, jet_primitive, jet_order, laurent_coefficients


class CocyclePrimitiveError(ValueError):
    """The generic lift or reconstruction could not be verified."""


class DarbouxPrimitiveRequired(CocyclePrimitiveError):
    """A nonzero weight-one value cocycle needs a state-space primitive."""

    def __init__(self, equations):
        self.equations = dict(equations)
        super().__init__("inhomogeneous Darboux primitive required for remaining value cocycle")


def _vertical_primitive(chart, form, weight, order):
    """Integrate a closed scalar positive-jet form; no state derivative is taken."""
    for coefficient in form.values():
        if jet_order(chart, coefficient) > order:
            raise CocyclePrimitiveError("cocycle exceeds the certified inverse-jet order")
        try:
            for _ in laurent_coefficients(chart, coefficient):
                pass
        except TransportUnresolved as error:
            raise CocyclePrimitiveError(str(error)) from error
    if any(not is_zero(jet_partial(chart, form[a], *b) - jet_partial(chart, form[b], *a))
           for a, b in combinations(form, 2)):
        raise CocyclePrimitiveError("the extracted vertical one-form is not closed")

    if weight > 0:
        # In fixed state, sum z*∂_z = (−weight + sum ∂_log_z) on homogeneous
        # expressions. The log derivative sum is nilpotent on polynomials.
        log_symbols = [chart.log(i) for i in range(chart.dimension)]
        term = canonical(sum(chart.jet(*a) * rhs for a, rhs in form.items()))
        completed = sp.S.Zero
        while term != 0:
            completed -= term / weight
            term = canonical(sum(term.diff(ell) for ell in log_symbols) / weight)
    else:
        # Weight zero has a nontrivial state-only ambiguity. Deterministically
        # integrate the residual in each jet direction with zero constants.
        completed = sp.S.Zero
        for (i, s), rhs in form.items():
            residual = canonical(rhs - jet_partial(chart, completed, i, s))
            try:
                q = jet_primitive(chart, residual, i, s)
            except TransportUnresolved as error:
                raise CocyclePrimitiveError("could not integrate the vertical Laurent–log form") from error
            completed = canonical(completed + q)

    completed = canonical(completed)
    if any(not is_zero(jet_partial(chart, completed, *a) - rhs) for a, rhs in form.items()):
        raise CocyclePrimitiveError("the computed vertical potential failed its gradient check")
    return completed


def primitive_from_lift(chart, family, source, lift, degree, *, darboux_solver=None):
    """Return F satisfying [v_i*u_i,x, F] = source, via a generic lift.

    The physical, homogeneous `lift` must use `family.values` and satisfy
    [v_i*u_i,x, lift] = [w_i*u_i,x, source] for every Tsarev speed.
    `degree` is the physical differential degree of source (at least two).
    The returned primitive has degree `degree - 1`. If the vertical
    weight-zero primitive leaves a value cocycle, caller may supply
    `darboux_solver(chart, rhs)`, with rhs keyed by (i,j), i != j. Its state
    vector h must satisfy ∂_j h_i − Γ_ij*(h_j−h_i) = rhs[i,j].
    """
    degree = index(degree)
    if degree < 2:
        raise CocyclePrimitiveError("source must be homogeneous of degree at least two")
    space = chart.space
    if tuple(space.states) != tuple(family.space.states) or any(
        not is_zero(a - b) for a, b in zip(chart.speeds, family.speeds)
    ) or len(chart.speeds) != len(family.speeds):
        raise CocyclePrimitiveError("Tsarev family must belong to the same leading system")
    if any(is_zero(v.diff(u)) for v, u in zip(chart.speeds, space.states)):
        raise CocyclePrimitiveError("leading speeds need nonvanishing self-derivatives")
    source, lift = tuple(map(canonical, source)), tuple(map(family.reduce, lift))
    if len(source) != space.dimension or len(lift) != space.dimension:
        raise CocyclePrimitiveError("source and lift must be characteristic vectors")
    if any(not space.is_homogeneous(expr, degree) for expr in source + lift):
        raise CocyclePrimitiveError("source and lift must have the given physical degree")
    if any(expr.has(*family.values) for expr in source):
        raise CocyclePrimitiveError("the homological source must not depend on the generic speed")
    order = max(1, *(space.max_order(expr) for expr in source))
    if any(space.max_order(expr) > order for expr in lift):
        raise CocyclePrimitiveError("lift exceeds the physical jet-order bound of the source")
    p = tuple(space.jet(i, 1) for i in range(space.dimension))
    X = tuple(v * px for v, px in zip(chart.speeds, p))
    W = tuple(w * px for w, px in zip(family.values, p))
    if any(family.reduce(a - b) != 0 for a, b in zip(
        space.bracket(X, lift), space.bracket(W, source)
    )):
        raise CocyclePrimitiveError("supplied generic lift failed its physical bracket identity")

    free = family.free_jets(order)
    basis = tuple(sp.Dummy(f"cocycle_free_{k}") for k in range(len(free)))
    normalize_free = dict(zip(free, basis))
    indices = tuple((j, s) for j in range(space.dimension) for s in range(1, order + 1))
    cocycle = chart.normalize(lift)
    vertical = []
    for expr in cocycle:
        try:
            polynomial = sp.Poly(expr.xreplace(normalize_free), *basis)
        except sp.PolynomialError as error:
            raise CocyclePrimitiveError("normalized lift is not polynomial in free speed jets") from error
        if any(sum(monomial) != 1 for monomial, coefficient in polynomial.terms()
               if coefficient != 0) or any(
            coefficient.has(*family.values) for coefficient in polynomial.coeffs()
        ):
            raise CocyclePrimitiveError("normalized lift is not linear in the certified free jets")
        form = {(j, s): canonical(-polynomial.coeff_monomial(basis[j * (order + 1) + s]))
                for j, s in indices}
        vertical.append(_vertical_primitive(chart, form, degree - 2, order))

    residual = tuple(family.reduce(c - a) for c, a in zip(
        cocycle, chart.lie_action(family.values, vertical)
    ))
    if any(not is_zero(expr) for expr in residual):
        if degree != 2:
            raise CocyclePrimitiveError("unexpected value cocycle above weight one")
        rhs = {}
        for i, expr in enumerate(residual):
            terms = []
            for j, wj in enumerate(family.values):
                if j == i:
                    continue
                z = chart.jet(j, 1)
                coefficient = canonical(z * sp.diff(expr, wj))
                if chart.active_symbols(coefficient) or coefficient.has(*family.values):
                    raise CocyclePrimitiveError("weight-one remainder is not a Darboux value cocycle")
                rhs[i, j] = coefficient
                terms.append((wj - family.values[i]) * coefficient / z)
            if not is_zero(expr - sum(terms)):
                raise CocyclePrimitiveError("weight-one remainder is not a Darboux value cocycle")
        if darboux_solver is None:
            raise DarbouxPrimitiveRequired(rhs)
        state = tuple(map(canonical, darboux_solver(chart, rhs)))
        if len(state) != space.dimension or any(
            chart.active_symbols(expr)
            or any(s > 0 for _, s in space.active_jets(expr))
            or expr.has(*family.values)
            for expr in state
        ):
            raise CocyclePrimitiveError("Darboux primitive must be a state-only vector")
        if any(not is_zero(
            state[i].diff(space.states[j]) - chart.tsarev(i, j) * (state[j] - state[i]) - f
        ) for (i, j), f in rhs.items()):
            raise CocyclePrimitiveError("supplied Darboux primitive failed its differential equations")
        vertical = [canonical(h + base) for h, base in zip(vertical, state)]

    if any(family.reduce(c - a) != 0 for c, a in zip(
        cocycle, chart.lie_action(family.values, vertical)
    )):
        raise CocyclePrimitiveError("potential does not reconstruct the full normalized lift")
    F = chart.denormalize(vertical)
    if any(not space.is_homogeneous(expr, degree - 1) for expr in F) or any(
        not is_zero(a - b) for a, b in zip(space.bracket(X, F), source)
    ):
        raise CocyclePrimitiveError("reconstructed primitive failed the final physical check")
    return F
