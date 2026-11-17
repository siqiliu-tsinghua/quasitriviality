"""Verified finite-span solvers for coupled homological equations.

This is not a general inverse for ad(X0), or the full smooth Darboux
construction. An unsuccessful finite search is not a nonexistence proof.
"""

from __future__ import annotations

import operator
from collections.abc import Iterable, Iterator, Sequence

import sympy as sp

from .jets import JetSpace, canonical, is_zero


class CoupledSolveUnresolved(ValueError):
    """No verified primitive was constructed in the selected search space."""


def _bound(value: int, label: str, minimum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be an integer")
    value = operator.index(value)
    if value < minimum:
        raise ValueError(f"{label} must be at least {minimum}")
    return value


def _weight(space: JetSpace, source: tuple[sp.Expr, ...]) -> int:
    nonzero = next(expr for expr in source if not is_zero(expr))
    euler = sum(
        s * space.jet(i, s) * sp.diff(nonzero, space.jet(i, s))
        for i, s in space.active_jets(nonzero)
        if s
    )
    weight = canonical(euler / nonzero)
    if not isinstance(weight, sp.Integer) or weight < 2:
        raise ValueError("source must have a homogeneous differential weight at least two")
    degree = int(weight)
    if not all(space.is_homogeneous(expr, degree) for expr in source):
        raise ValueError("all source components must have the same differential weight")
    return degree


def _jet_monomials(
    space: JetSpace, weight: int, max_order: int, max_inverse_power: int
) -> Iterator[sp.Expr]:
    # Lower bounds permit positive higher jet powers to be balanced by
    # bounded negative powers of first jets.
    entries = [
        (space.jet(i, order), order, -max_inverse_power if order == 1 else 0)
        for order in range(max_order, 0, -1)
        for i in range(space.dimension)
    ]
    minima = [0] * (len(entries) + 1)
    for k in range(len(entries) - 1, -1, -1):
        _, degree, lower = entries[k]
        minima[k] = minima[k + 1] + degree * lower

    def visit(k: int, remaining: int, powers: list[sp.Expr]) -> Iterator[sp.Expr]:
        if k == len(entries):
            if remaining == 0:
                yield sp.Mul(*powers)
            return
        symbol, degree, lower = entries[k]
        upper = (remaining - minima[k + 1]) // degree
        for exponent in range(lower, upper + 1):
            if exponent:
                powers.append(symbol**exponent)
            yield from visit(k + 1, remaining - degree * exponent, powers)
            if exponent:
                powers.pop()

    yield from visit(0, weight, [])


def _log_monomials(space: JetSpace, limit: int) -> Iterator[sp.Expr]:
    def visit(i: int, remaining: int, powers: list[sp.Expr]) -> Iterator[sp.Expr]:
        if i == space.dimension:
            yield sp.Mul(*powers)
            return
        for exponent in range(remaining + 1):
            if exponent:
                powers.append(space.log(i)**exponent)
            yield from visit(i + 1, remaining - exponent, powers)
            if exponent:
                powers.pop()

    yield from visit(0, limit, [])


def _split_identity(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> list[sp.Expr]:
    """Split a denominator-cleared identity by variable-dependent factors.

    Exact coefficient comparison works in a rational function ring.
    Nonrational factors can be treated as distinct, so the search may
    miss a solution even in the supplied span. Verification prevents
    false positives.
    """
    numerator = sp.fraction(sp.together(expr))[0]
    coefficients: dict[sp.Expr, sp.Expr] = {}
    for term in sp.Add.make_args(sp.expand(numerator)):
        coefficient, factor = term.as_independent(*variables, as_Add=False)
        coefficients[factor] = coefficients.get(factor, sp.S.Zero) + coefficient
    return [value for value in map(canonical, coefficients.values()) if not is_zero(value)]


class CoupledSolver:
    """Solve [X0, F] = source in an explicit finite homogeneous span.

    Default candidates are constant-coefficient jet polynomials, with
    derivative order bounded by the requested generator weight. Override
    max_jet_order for a different bound; max_inverse_power bounds the
    negative exponent of each first jet, and max_log_degree bounds the
    total formal-log exponent. Each allowed jet monomial is multiplied
    by all entries of state_basis. extra_candidates supplies particular
    vectors, including vectors whose components share a single coefficient.
    Candidates of another weight are used only at that weight.

    By default compare identities in the generalized hodograph chart
    of the supplied semihamiltonian diagonal leading system. All
    results are verified in physical jets. Work on a domain where the
    speed gaps, self-derivatives and result denominators are nonzero;
    a symbolic parameter in a solved denominator must be specialized
    away from its zeros.
    """

    def __init__(
        self,
        *,
        state_basis: Iterable[object] = (1,),
        max_jet_order: int | None = None,
        max_inverse_power: int = 0,
        max_log_degree: int = 0,
        extra_candidates: Iterable[Sequence[object]] = (),
        max_candidates: int = 400,
        use_hodograph: bool = True,
    ) -> None:
        self.state_basis = tuple(canonical(term) for term in state_basis)
        self.max_jet_order = (
            None if max_jet_order is None else _bound(max_jet_order, "max_jet_order", 1)
        )
        self.max_inverse_power = _bound(max_inverse_power, "max_inverse_power", 0)
        self.max_log_degree = _bound(max_log_degree, "max_log_degree", 0)
        self.extra_candidates = tuple(tuple(map(canonical, row)) for row in extra_candidates)
        self.max_candidates = _bound(max_candidates, "max_candidates", 1)
        self.use_hodograph = use_hodograph

    def _candidates(
        self, space: JetSpace, weight: int
    ) -> Iterator[tuple[sp.Expr, ...]]:
        if any(space.max_order(coefficient) for coefficient in self.state_basis):
            raise ValueError("state_basis entries may not depend on positive-order jets")
        zero = (sp.S.Zero,) * space.dimension
        for candidate in self.extra_candidates:
            if len(candidate) != space.dimension:
                raise ValueError("extra candidate has the wrong number of components")
            if candidate != zero and all(
                space.is_homogeneous(expr, weight) for expr in candidate
            ):
                yield candidate

        basis = tuple(term for term in self.state_basis if not is_zero(term))
        if not basis:
            return
        max_order = self.max_jet_order or weight
        for monomial in _jet_monomials(space, weight, max_order, self.max_inverse_power):
            for logs in _log_monomials(space, self.max_log_degree):
                for coefficient in basis:
                    expr = canonical(monomial * logs * coefficient)
                    for i in range(space.dimension):
                        yield zero[:i] + (expr,) + zero[i + 1 :]

    def __call__(
        self, space: JetSpace, speeds: Sequence[object], source: Sequence[object]
    ) -> tuple[sp.Expr, ...]:
        if space.dimension < 2:
            raise ValueError("use the scalar solver for a one-component equation")
        speeds = tuple(map(canonical, speeds))
        source = tuple(map(canonical, source))
        if len(speeds) != space.dimension or len(source) != space.dimension:
            raise ValueError("leading speeds and source must match the jet space dimension")
        if any(space.max_order(speed) for speed in speeds):
            raise ValueError("leading speeds may only depend on the state")
        if any(expr.has(sp.Float) for expr in (*speeds, *source, *self.state_basis)):
            raise ValueError("the exact solver does not accept floating-point coefficients")
        if all(is_zero(expr) for expr in source):
            return (sp.S.Zero,) * space.dimension
        weight = _weight(space, source) - 1
        candidates = []
        seen = set()
        for candidate in self._candidates(space, weight):
            if candidate in seen:
                continue
            if any(expr.has(sp.Float) for expr in candidate):
                raise ValueError("the exact solver does not accept floating-point candidates")
            seen.add(candidate)
            candidates.append(candidate)
            if len(candidates) > self.max_candidates:
                raise CoupledSolveUnresolved(
                    f"weight {weight} search exceeds its {self.max_candidates}-candidate bound"
                )
        if not candidates:
            raise CoupledSolveUnresolved(f"empty search space at generator weight {weight}")

        parameters = set()
        for expr in (*speeds, *source, *(x for row in candidates for x in row)):
            parameters.update(
                expr.free_symbols
                - {space.jet(i, s) for i, s in space.active_jets(expr)}
                - {space.log(i) for i in space.active_logs(expr)}
            )
        parameters.difference_update(space.states)
        leading = tuple(speed * space.jet(i, 1) for i, speed in enumerate(speeds))
        if self.use_hodograph:
            from .hodograph import HodographChart

            chart = HodographChart(space, speeds)
            target = chart.normalize(source)
            # Compute matrix columns by bracket conjugacy. Applying the
            # transported derivation to each rational inverse-jet expression
            # separately creates much larger intermediate numerators; both
            # forms give the same normalized operator.
            columns = [chart.normalize(space.bracket(leading, row)) for row in candidates]
        else:
            target = source
            columns = [space.bracket(leading, candidate) for candidate in candidates]
        unknowns = tuple(sp.Dummy(f"c_{k}") for k in range(len(candidates)))
        unknown_set = set(unknowns)
        equations = []
        for i, component in enumerate(target):
            residual = sum(a * column[i] for a, column in zip(unknowns, columns)) - component
            variables = tuple(sorted(
                residual.free_symbols - unknown_set - parameters, key=sp.default_sort_key
            ))
            equations.extend(_split_identity(residual, variables))

        if not equations:
            # Even a vacuous coefficient solve gets verified independently.
            coefficients = (sp.S.Zero,) * len(candidates)
        else:
            matrix, values = sp.linear_eq_to_matrix(equations, unknowns)
            solution = sp.linsolve((matrix, values), unknowns)
            if solution == sp.EmptySet:
                raise CoupledSolveUnresolved(
                    f"no exact witness found in the selected weight {weight} span"
                )
            coefficients = next(iter(solution))
            free_zero = dict.fromkeys(unknowns, sp.S.Zero)
            coefficients = tuple(term.xreplace(free_zero) for term in coefficients)

        answer = tuple(
            canonical(sum(c * candidate[i] for c, candidate in zip(coefficients, candidates)))
            for i in range(space.dimension)
        )
        if any(
            not is_zero(actual - expected)
            for actual, expected in zip(space.bracket(leading, answer), source)
        ):
            raise CoupledSolveUnresolved("candidate solution did not verify in physical jets")
        return answer
