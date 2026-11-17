"""Deterministic normalized scalar transport inversion in hodograph jets."""

from functools import lru_cache
from operator import index

import sympy as sp

from .jets import _expr, canonical, is_zero


class TransportUnresolved(ValueError):
    """A required finite Laurent–log transport identity was not solved."""


def laurent_coefficients(chart, expr):
    """Yield (nonzero exponents keyed by auxiliary symbol, scalar coefficient)."""
    expr = canonical(expr)
    if expr == 0:
        return
    variables = chart.active_symbols(expr)
    if not variables:
        yield {}, expr
        return
    allowed_inverse = {chart.jet(i, 1) for i in range(chart.dimension)}
    try:
        numerator, denominator = sp.fraction(sp.together(expr))
        denom = sp.Poly(denominator, *variables, domain=sp.EX)
        if len(denom.terms()) != 1:
            raise TransportUnresolved("nonmonomial jet denominator")
        shift, factor = denom.terms()[0]
        if any(e and q not in allowed_inverse for q, e in zip(variables, shift)):
            raise TransportUnresolved("non-first-order jet denominator")
        for powers, c in sp.Poly(numerator, *variables, domain=sp.EX).terms():
            coefficient = canonical(c / factor)
            if coefficient != 0:
                yield {q: int(p - d) for q, p, d in zip(variables, powers, shift) if p != d}, coefficient
    except sp.PolynomialError as exc:
        raise TransportUnresolved("expected finite Laurent–log support") from exc


def jet_order(chart, expr):
    """Largest inverse jet order; a formal log counts as first order."""
    return max((chart._coordinate_indices.get(q, (None, 1))[1]
                for q in chart.active_symbols(expr)), default=0)


def jet_partial(chart, expr, j, order):
    """Jet derivative including the first-jet formal log rule."""
    expr = _expr(expr)
    z = chart.jet(j, order)
    derivative = expr.diff(z)
    if order == 1:
        derivative += expr.diff(chart.log(j)) / z
    return canonical(derivative)


def jet_primitive(chart, expr, j, order):
    """Integrate a finite Laurent–log expression; set the constant to zero."""
    z = chart.jet(j, order)
    ell = chart.log(j) if order == 1 else None
    result = sp.S.Zero
    for powers, c in laurent_coefficients(chart, expr):
        n = powers.pop(z, 0)
        if order > 1:
            result += c * z**(n + 1) * sp.prod(q**e for q, e in powers.items()) / (n + 1)
            continue
        k = powers.pop(ell, 0)
        c *= sp.prod(q**e for q, e in powers.items())
        if n == -1:
            result += c * ell**(k + 1) / (k + 1)
        else:
            m = sp.Integer(n + 1)
            polynomial = sum((-1)**d * sp.factorial(k) * ell**(k - d)
                             / sp.factorial(k - d) / m**(d + 1)
                             for d in range(k + 1))
            result += c * z**m * polynomial
    return canonical(result)


def own_projection(chart, expr, i, exponent):
    """Coefficient of own first jet to exponent, allowing its log powers only."""
    z, ell = chart.jet(i, 1), chart.log(i)
    return canonical(sum(c * ell**powers.get(ell, 0)
                         for powers, c in laurent_coefficients(chart, expr)
                         if powers.get(z, 0) == exponent and powers.keys() <= {z, ell}))


def transport_alpha(chart, i):
    """The transverse logarithmic slope in the ith normalized action."""
    return canonical(sum(chart.speeds[i].diff(chart.states[j]) / chart.jet(j, 1)
                         for j in range(chart.dimension) if j != i))


def first_primitive(chart, i, rhs, weight):
    """Recover the forced own first-jet part; remove the weight-zero gauge."""
    if weight < 0:
        raise ValueError("require nonnegative weight")
    z, ell = chart.jet(i, 1), chart.log(i)
    slope = chart.speeds[i].diff(chart.states[i])
    if is_zero(slope):
        raise TransportUnresolved("zero self derivative")
    pure = own_projection(chart, rhs, i, -weight - 1)
    if weight == 0:
        return canonical(-sp.integrate(pure, ell) / slope)
    poly = sp.Poly(pure, ell, domain=sp.EX)
    result = sum(c * sp.factorial(k) * ell**(k - d)
                 / sp.factorial(k - d) / sp.Integer(weight)**(d + 1)
                 for (k,), c in poly.terms() for d in range(k + 1))
    return canonical(z**(-weight) * result / slope)


def transverse_primitive(chart, i, rhs, n):
    """Recover forced dependence on order-n transverse inverse jets."""
    forced = {}
    for j in range(chart.dimension):
        if j != i:
            gap = chart.speeds[j] - chart.speeds[i]
            forced[j] = canonical(rhs.diff(chart.jet(j, n + 1)) * chart.jet(j, 1) / gap)
            if jet_order(chart, forced[j]) > n:
                raise TransportUnresolved("top transverse derivative is nonlinear in the next jet")
    result = sp.S.Zero
    for j, f in forced.items():
        result += jet_primitive(chart, f - jet_partial(chart, result, j, n), j, n)
    for j, f in forced.items():
        if not is_zero(jet_partial(chart, result, j, n) - f):
            raise TransportUnresolved("incompatible top transverse jet derivatives")
    return canonical(result)


def transport_operator(chart, i, expr, alpha_factor=1):
    """Evaluate (E_v - v_i D_x + alpha_factor*alpha_i) on inverse jets."""
    return _Transport(chart, i).apply(_expr(expr), alpha_factor)


class _Transport:
    def __init__(self, chart, i):
        self.chart = chart
        self.i = i
        self.states = chart.states
        self.speeds = chart.speeds
        self.n = chart.dimension
        self.gaps = tuple(canonical(v - self.speeds[i]) for v in self.speeds)
        self.alpha = transport_alpha(chart, i)
        self.logs = tuple(chart.log(j) for j in range(self.n))

    @lru_cache(maxsize=256)
    def _velocity(self, j, s):
        c = self.chart
        return canonical(
            self.gaps[j] * c.jet(j, s + 1) / c.jet(j, 1)
            + sum(
                self.gaps[l] / c.jet(l, 1)
                * sum(c.jet_coefficient(j, l, s, t) * c.jet(j, t + 1)
                      for t in range(s))
                for l in range(self.n) if l != j and l != self.i
            )
            - self.speeds[j].diff(self.states[j], s)
        )

    def apply(self, expr, factor):
        c = self.chart
        symbols = expr.free_symbols
        terms = [factor * self.alpha * expr]
        for j in range(self.n):
            if j != self.i and self.states[j] in symbols:
                terms.append(self.gaps[j] * expr.diff(self.states[j]) / c.jet(j, 1))
            if self.logs[j] in symbols:
                terms.append(self._velocity(j, 1) * expr.diff(self.logs[j]) / c.jet(j, 1))
        for symbol in c.active_symbols(expr):
            if symbol not in self.logs:
                j, s = c._coordinate_indices[symbol]
                terms.append(self._velocity(j, s) * expr.diff(symbol))
        return canonical(sum(terms))

    def solve(self, rhs, weight, order, factor):
        rhs = canonical(rhs)
        if rhs == 0:
            return sp.S.Zero
        if order < 1:
            raise TransportUnresolved("nonzero source has no positive-weight jet-free primitive")
        if jet_order(self.chart, rhs) > order + 1:
            raise TransportUnresolved("transport source exceeds allowed order")
        top = transverse_primitive(self.chart, self.i, rhs, order)
        rest = canonical(rhs - self.apply(top, factor))
        if jet_order(self.chart, rest) > order:
            raise TransportUnresolved("uncancelled transverse highest jets")
        if order == 1:
            result = canonical(top + first_primitive(self.chart, self.i, rest, weight))
        else:
            y = self.chart.jet(self.i, order)
            try:
                poly = sp.Poly(rest, y, domain=sp.EX)
            except sp.PolynomialError as exc:
                raise TransportUnresolved("nonpolynomial own higher inverse jet") from exc
            for power in range(int(poly.degree()) if not poly.is_zero else 0, 0, -1):
                leading = sp.Poly(rest, y, domain=sp.EX).nth(power)
                lower = self.solve(leading, weight + power, order - 1, factor - power)
                term = lower * y**power
                rest = canonical(rest - self.apply(term, factor))
                if not is_zero(sp.Poly(rest, y, domain=sp.EX).nth(power)):
                    raise TransportUnresolved("own-jet coefficient did not cancel")
                top += term
            if rest.has(y):
                raise TransportUnresolved("uncancelled highest own jet")
            result = canonical(top + self.solve(rest, weight, order - 1, factor))
        if not is_zero(self.apply(result, factor) - rhs):
            raise TransportUnresolved("candidate did not satisfy the complete transport equation")
        return result


def solve_transport(chart, i, rhs, weight, max_order, alpha_factor=1):
    """Solve a finite Laurent–log transport equation of positive q-degree.

    Each inverse jet has weight -1. max_order bounds the unknown, not the
    source, which may contain one additional transverse derivative.
    """
    weight, max_order = index(weight), index(max_order)
    if weight <= 0 or max_order < 1:
        raise ValueError("require positive transport weight and jet order")
    return _Transport(chart, i).solve(_expr(rhs), weight, max_order, alpha_factor)
