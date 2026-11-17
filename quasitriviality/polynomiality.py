"""Finite-order polynomiality in physical jets, not in inverse jets.

Supported inputs are rational expressions that reduce to finite Laurent
polynomials in the first jets, polynomial in higher jets and the jet space's
formal logarithms, with exact jet-independent coefficients. General rational
or nonalgebraic jet dependence is reported as unresolved, not as a proof of
nonintegrability.
"""
from dataclasses import dataclass
from enum import Enum
from typing import Callable, Iterable

import sympy as sp

from .jets import canonical, is_zero


class PolynomialityStatus(str, Enum):
    POLYNOMIAL = "polynomial"
    OBSTRUCTED = "obstructed"
    CONDITIONAL = "conditional"
    UNDETERMINED = "undetermined"


@dataclass(frozen=True)
class JetCondition:
    """The coefficient of a bad monomial must vanish."""

    order: int
    component: int
    monomial: sp.Expr
    coefficient: sp.Expr
    speed_monomial: sp.Expr = sp.S.One
    certified_nonzero: bool = False


@dataclass(frozen=True)
class UnresolvedCoefficient:
    order: int
    component: int
    expression: sp.Expr
    reason: str


@dataclass(frozen=True)
class PolynomialityReport:
    status: PolynomialityStatus
    cutoff: int
    conditions: tuple[JetCondition, ...]
    unresolved: tuple[UnresolvedCoefficient, ...]


# These objects might be jet-independent as represented, but we must not
# certify expressions carrying an unevaluated primitive or approximation.
_UNEVALUATED = (sp.Integral, sp.Sum, sp.Product, sp.Limit, sp.Order)
_BAD_NUMBERS = (sp.nan, sp.zoo, sp.oo, -sp.oo)


def _split_free_data(coefficient, independent, dummies):
    """Coefficient splitting in caller-certified independent speed data."""
    if not independent:
        return [(sp.S.One, coefficient)]
    expr = coefficient.xreplace(dict(zip(independent, dummies)))
    poly = sp.Poly(expr, *dummies, domain="EX")
    return [
        (sp.prod(x**p for x, p in zip(independent, powers)), canonical(value))
        for powers, value in poly.terms()
    ]


def analyze_polynomiality(
    space,
    series,
    *,
    independent: Iterable[sp.Expr] = (),
    coefficient_reducer: Callable[[sp.Expr], sp.Expr] | None = None,
) -> PolynomialityReport:
    """Analyze each component of a physical-jet Series through its cutoff.

    ``coefficient_reducer`` is an optional exact rewrite (for example,
    elimination of transverse Tsarev derivatives) applied before testing.
    ``independent`` lists free speed data, such as ``h(u)`` and its pure
    derivatives. Splitting assumes its entries can take independent values.
    Without that declaration, a symbolic bad coefficient gives a conditional
    identity, not a certified nonzero obstruction. Component and order numbers
    are zero-based.
    """
    if series.dimension != space.dimension:
        raise ValueError("series and jet space dimensions must agree")
    independent = tuple(sp.sympify(x, strict=True) for x in independent)
    if len(set(independent)) != len(independent):
        raise ValueError("independent speed data must be distinct")
    dummies = tuple(sp.Dummy("free_speed") for _ in independent)
    conditions = []
    unresolved = []

    def unclear(order, component, expression, reason):
        unresolved.append(UnresolvedCoefficient(order, component, expression, reason))

    for order in range(series.cutoff):
        for component, coefficient in enumerate(series.coefficient(order)):
            expr = sp.sympify(coefficient, strict=True)
            if coefficient_reducer is not None:
                expr = sp.sympify(coefficient_reducer(expr), strict=True)
            if expr.has(*_UNEVALUATED, *_BAD_NUMBERS, sp.Float):
                unclear(order, component, expr, "unevaluated, singular, or inexact coefficient")
                continue
            expr = canonical(expr)
            if expr.has(*_UNEVALUATED, *_BAD_NUMBERS, sp.Float):
                unclear(order, component, expr, "unevaluated, singular, or inexact reduction")
                continue
            if is_zero(expr):
                continue
            jets = tuple(
                space.jet(i, s) for i, s in sorted(set(space.active_jets(expr))) if s > 0
            )
            logs = tuple(space.log(i) for i in range(space.dimension) if expr.has(space.log(i)))
            variables = jets + logs
            if not variables:
                continue
            first = {space.jet(i, 1) for i in range(space.dimension)}
            log_count = len(logs)
            try:
                num, den = sp.fraction(sp.cancel(expr))
                numpoly = sp.Poly(num, *variables, domain="EX")
                denpoly = sp.Poly(den, *variables, domain="EX")
            except (sp.PolynomialError, TypeError, NotImplementedError, ValueError):
                unclear(order, component, expr, "not in the supported Laurent–log class")
                continue
            den_terms = denpoly.terms()
            if len(den_terms) != 1:
                unclear(order, component, expr, "denominator is not a first-jet monomial")
                continue
            den_powers, den_coeff = den_terms[0]
            if any(p and x not in first for x, p in zip(variables, den_powers)):
                unclear(order, component, expr, "denominator contains a higher jet or a logarithm")
                continue

            for num_powers, num_coeff in numpoly.terms():
                powers = tuple(a - b for a, b in zip(num_powers, den_powers))
                if min(powers) >= 0 and (not log_count or not any(powers[-log_count:])):
                    continue
                monomial = sp.prod(x**p for x, p in zip(variables, powers))
                value = canonical(num_coeff / den_coeff)
                try:
                    split = _split_free_data(value, independent, dummies)
                except (sp.PolynomialError, TypeError, NotImplementedError, ValueError):
                    unclear(order, component, value * monomial,
                            "bad coefficient is not polynomial in the declared free speed data")
                    continue
                for speed_monomial, scalar in split:
                    if not is_zero(scalar):
                        conditions.append(JetCondition(
                            order, component, monomial, scalar, speed_monomial,
                            certified_nonzero=scalar.is_zero is False,
                        ))

    if any(c.certified_nonzero for c in conditions):
        status = PolynomialityStatus.OBSTRUCTED
    elif unresolved:
        status = PolynomialityStatus.UNDETERMINED
    elif conditions:
        status = PolynomialityStatus.CONDITIONAL
    else:
        status = PolynomialityStatus.POLYNOMIAL
    return PolynomialityReport(status, series.cutoff, tuple(conditions), tuple(unresolved))
