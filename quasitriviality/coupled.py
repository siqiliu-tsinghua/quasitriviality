"""Constructive finite-jet inverse of the coupled homological operator."""
from __future__ import annotations

import sympy as sp

from .coupled_transport import (
    first_primitive, jet_order, solve_transport, transverse_primitive, TransportUnresolved,
)
from .darboux import DarbouxUnresolved, solve_darboux
from .hodograph import HodographChart
from .jets import canonical, is_zero


class CoupledSolveUnresolved(ValueError):
    """The source did not yield a verified local Laurent–log primitive."""


def _weight(space, source):
    nonzero = next(expr for expr in source if not is_zero(expr))
    euler = sum(s * space.jet(i, s) * nonzero.diff(space.jet(i, s))
                for i, s in space.active_jets(nonzero) if s)
    weight = canonical(euler / nonzero)
    if not isinstance(weight, sp.Integer) or weight < 2:
        raise ValueError("source must have a homogeneous differential weight at least two")
    degree = int(weight)
    if not all(space.is_homogeneous(expr, degree) for expr in source):
        raise ValueError("all source components must have the same differential weight")
    return degree


def _state_forcing(chart, residual):
    """Read the z_j,1^(-1) residual and check its exact Darboux shape."""
    forcing = {}
    dimension = chart.dimension
    jets = tuple(chart.jet(j, 1) for j in range(dimension))
    for i, component in enumerate(residual):
        if is_zero(component):
            for j in range(dimension):
                if j != i:
                    forcing[i, j] = sp.S.Zero
            continue
        # If the shape is right, multiplying by the first jets gives a
        # polynomial with exactly one distinct monomial for each j != i.
        reduced = canonical(component * sp.prod(jets))
        try:
            poly = sp.Poly(reduced, *jets, domain="EX")
        except sp.PolynomialError as exc:
            raise CoupledSolveUnresolved("residual is not of Darboux form") from exc
        for j in range(dimension):
            if j == i:
                continue
            powers = tuple(0 if a == j else 1 for a in range(dimension))
            coeff = poly.coeff_monomial(powers)
            if chart.active_symbols(coeff):
                raise CoupledSolveUnresolved("Darboux forcing still contains jets or logarithms")
            forcing[i, j] = canonical(coeff / (chart.speeds[j] - chart.speeds[i]))
        predicted = sum((chart.speeds[j] - chart.speeds[i]) * forcing[i, j] / jets[j]
                        for j in range(dimension) if i != j)
        if not is_zero(component - predicted):
            raise CoupledSolveUnresolved("residual is not of Darboux form")
    return forcing


class CoupledSolver:
    """Integrate a coupled homological equation in generalized hodograph jets.

    No candidate generators are searched. All finite jet orders and polynomial
    degrees are dictated by the source. At degree one the residual state part
    is a locally compatible affine Darboux problem. Optional `darboux_solver`
    supplies an exact state-space primitive, not a jet-space basis.
    """

    def __init__(self, *, darboux_solver=None, basepoint=None):
        self.darboux_solver = darboux_solver
        self.basepoint = basepoint

    def __call__(self, space, speeds, source):
        speeds, source = tuple(map(canonical, speeds)), tuple(map(canonical, source))
        if space.dimension < 2:
            raise ValueError("use the scalar solver for a one-component equation")
        if len(speeds) != space.dimension or len(source) != space.dimension:
            raise ValueError("leading speeds and source must match the jet space dimension")
        if any(expr.has(sp.Float) for expr in (*speeds, *source)):
            raise ValueError("the exact solver does not accept floating-point coefficients")
        if all(is_zero(expr) for expr in source):
            return (sp.S.Zero,) * space.dimension
        m = _weight(space, source) - 2
        chart = HodographChart(space, speeds)
        if any(not is_zero(v) for v in chart.semi_hamiltonian_residual().values()):
            raise ValueError("leading speeds must satisfy the semi-Hamiltonian identities")
        if any(is_zero(v.diff(u)) for v, u in zip(speeds, chart.states)):
            raise ValueError("leading self-derivatives must be nonzero")
        residual = chart.normalize(source)
        answer = (sp.S.Zero,) * chart.dimension
        nmax = max(1, *(jet_order(chart, item) for item in residual))

        def subtract(parts):
            nonlocal residual, answer
            if all(is_zero(h) for h in parts):
                return
            action = chart.lie_action(speeds, parts)
            residual = tuple(canonical(r - a) for r, a in zip(residual, action))
            answer = tuple(canonical(h + x) for h, x in zip(answer, parts))

        try:
            for n in range(nmax, 0, -1):
                if any(jet_order(chart, item) > n + 1
                       for item in residual):
                    raise CoupledSolveUnresolved("residual exceeds the homological order bound")
                subtract(tuple(transverse_primitive(chart, i, r, n)
                               for i, r in enumerate(residual)))
                if any(jet_order(chart, item) > n for item in residual):
                    raise CoupledSolveUnresolved("top transverse derivative did not cancel")
                if n == 1:
                    subtract(tuple(first_primitive(chart, i, r, m)
                                   for i, r in enumerate(residual)))
                    continue
                # In row i no lower-order term or other unknown component
                # contains z[i,n]. Determine all positive own powers.
                for i in range(chart.dimension):
                    z = chart.jet(i, n)
                    while residual[i].has(z):
                        polynomial = sp.Poly(residual[i], z, domain="EX")
                        k = int(polynomial.degree())
                        if k <= 0:
                            raise CoupledSolveUnresolved("unsupported dependence on an own top jet")
                        term = z**k * solve_transport(
                            chart, i, polynomial.LC(), m + k, n - 1, alpha_factor=1-k
                        )
                        subtract(tuple(term if j == i else sp.S.Zero
                                       for j in range(chart.dimension)))
                        if residual[i].has(z) and sp.Poly(
                                residual[i], z, domain="EX").degree() >= k:
                            raise CoupledSolveUnresolved("own leading jet power did not cancel")
        except (TransportUnresolved, sp.PolynomialError) as exc:
            raise CoupledSolveUnresolved(str(exc)) from exc

        if m == 0 and any(not is_zero(r) for r in residual):
            forcing = _state_forcing(chart, residual)
            try:
                h = solve_darboux(chart, forcing, solver=self.darboux_solver, basepoint=self.basepoint)
            except DarbouxUnresolved as exc:
                raise CoupledSolveUnresolved(str(exc)) from exc
            subtract(h)
        if any(not is_zero(r) for r in residual):
            raise CoupledSolveUnresolved("normalized bracket has a nonzero residual")
        result = chart.denormalize(answer)
        leading = tuple(speed * space.jet(i, 1) for i, speed in enumerate(speeds))
        if (any(not space.is_homogeneous(x, m+1) for x in result)
                or any(not is_zero(a-b) for a, b in zip(space.bracket(leading, result), source))):
            raise CoupledSolveUnresolved("primitive did not verify in physical jets")
        return result
