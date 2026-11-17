"""Exact finite jet calculus for local evolutionary characteristics."""

from __future__ import annotations

import operator
import re
from collections.abc import Iterable, Sequence

import sympy as sp


def _expr(value: object) -> sp.Expr:
    if isinstance(value, (str, bytes)):
        raise TypeError("pass an explicit SymPy expression, not expression text")
    expr = sp.sympify(value, strict=True)
    if not isinstance(expr, sp.Expr) or expr.is_commutative is not True:
        raise TypeError("a jet coefficient must be a commutative scalar expression")
    return expr


def canonical(expr: object) -> sp.Expr:
    """Expand and cancel exact algebraic factors without changing log branches."""
    expr = sp.expand(_expr(expr))
    if expr.is_number or expr.is_polynomial(*expr.free_symbols):
        return expr
    return sp.expand(sp.cancel(expr))


def is_zero(expr: object) -> bool:
    """Return True only when exact symbolic manipulation establishes zero."""
    expr = canonical(expr)
    return expr == 0 or expr.is_zero is True


def _index(value: int, label: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{label} must be an integer") from exc
    if result < 0:
        raise ValueError(f"{label} must be nonnegative")
    return result


class JetSpace:
    """Jet symbols of finitely many scalar fields and their formal first-jet logs.

    A state can be a name or an explicit commutative SymPy Symbol.
    Use the factories rather than hand-creating jet and formal-log symbols.
    """

    def __init__(self, names: Iterable[str | sp.Symbol] | str | sp.Symbol):
        if isinstance(names, (str, sp.Symbol)):
            names = (names,)
        states = tuple(sp.Symbol(s) if isinstance(s, str) else s for s in names)
        if not states or any(
            not isinstance(s, sp.Symbol) or s.is_commutative is not True or not s.name
            for s in states
        ):
            raise ValueError("provide at least one nonempty commutative state symbol")
        labels = {s.name for s in states}
        if len(labels) != len(states):
            raise ValueError("state names must be distinct")
        for s in states:
            if f"L_{s.name}" in labels or any(
                re.fullmatch(re.escape(s.name) + r"_[1-9][0-9]*", label)
                for label in labels
            ):
                raise ValueError("state name collides with a jet or formal-log name")
        self.states: tuple[sp.Symbol, ...] = states
        self.dimension = len(states)
        self._jet_symbols = {(i, 0): s for i, s in enumerate(states)}
        self._log_symbols = tuple(sp.Symbol(f"L_{s.name}") for s in states)
        self._from_symbol = {s: (i, 0) for i, s in enumerate(states)}

    def _component(self, i: int) -> int:
        i = _index(i, "component")
        if i >= self.dimension:
            raise IndexError(f"no component {i} in dimension {self.dimension}")
        return i

    def jet(self, i: int, s: int) -> sp.Symbol:
        """Return the symbol for the order-s jet of zero-indexed component i."""
        i, s = self._component(i), _index(s, "jet order")
        pair = (i, s)
        if pair not in self._jet_symbols:
            symbol = sp.Symbol(f"{self.states[i].name}_{s}")
            self._jet_symbols[pair] = symbol
            self._from_symbol[symbol] = pair
        return self._jet_symbols[pair]

    def log(self, i: int) -> sp.Symbol:
        """Return the formal generator representing log(abs(first jet))."""
        return self._log_symbols[self._component(i)]

    def active_logs(self, expr: object) -> tuple[int, ...]:
        """Return components whose formal log symbol occurs in an expression."""
        free = _expr(expr).free_symbols
        return tuple(i for i, symbol in enumerate(self._log_symbols) if symbol in free)

    def active_jets(self, expr: object) -> tuple[tuple[int, int], ...]:
        """Return occurring (component, order) pairs, including log dependencies."""
        free = _expr(expr).free_symbols
        pairs = {self._from_symbol[s] for s in free if s in self._from_symbol}
        pairs.update((i, 1) for i in self.active_logs(expr))
        return tuple(sorted(pairs))

    def max_order(self, expr: object) -> int:
        return max((s for _, s in self.active_jets(expr)), default=0)

    def partial(self, expr: object, i: int, s: int) -> sp.Expr:
        """Jet partial derivative, with the formal-log chain rule at order one."""
        expr = _expr(expr)
        target = self.jet(i, s)
        result = sp.diff(expr, target)
        if s == 1 and self.log(i) in expr.free_symbols:
            result += sp.diff(expr, self.log(i)) / target
        return result

    def dx(self, expr: object, times: int = 1) -> sp.Expr:
        """Total spatial derivative, applied times times."""
        result = _expr(expr)
        for _ in range(_index(times, "derivative count")):
            result = canonical(
                sum(self.jet(i, s + 1) * self.partial(result, i, s)
                    for i, s in self.active_jets(result))
            )
        return result

    def _characteristic(self, vector: Sequence[object]) -> tuple[sp.Expr, ...]:
        result = tuple(_expr(f) for f in vector)
        if len(result) != self.dimension:
            raise ValueError(f"expected a characteristic with {self.dimension} components")
        return result

    def _evolve(
        self,
        vector: tuple[sp.Expr, ...],
        expr: object,
        derivatives: dict[tuple[int, int], sp.Expr],
    ) -> sp.Expr:
        expr = _expr(expr)
        terms = []
        for i, s in self.active_jets(expr):
            key = (i, s)
            if key not in derivatives:
                for order in range(1, s + 1):
                    if (i, order) not in derivatives:
                        derivatives[i, order] = self.dx(derivatives[i, order - 1])
            terms.append(derivatives[key] * self.partial(expr, i, s))
        return sum(terms, sp.S.Zero)

    def evolutionary(self, Q: Sequence[object], expr: object) -> sp.Expr:
        """Apply the evolutionary derivation induced by Q."""
        Q = self._characteristic(Q)
        derivatives = {(i, 0): q for i, q in enumerate(Q)}
        return canonical(self._evolve(Q, expr, derivatives))

    def bracket(self, Q: Sequence[object], P: Sequence[object]) -> tuple[sp.Expr, ...]:
        """Return E_Q(P) - E_P(Q)."""
        Q, P = self._characteristic(Q), self._characteristic(P)
        q_cache = {(i, 0): q for i, q in enumerate(Q)}
        p_cache = {(i, 0): p for i, p in enumerate(P)}
        return tuple(
            canonical(self._evolve(Q, p, q_cache) - self._evolve(P, q, p_cache))
            for q, p in zip(Q, P)
        )

    def is_homogeneous(self, expr: object, weight: int) -> bool:
        """Check physical differential weight by scale covariance, fixing formal logs."""
        expr = canonical(expr)
        if isinstance(weight, bool):
            raise TypeError("weight must be an integer")
        weight = operator.index(weight)
        t = sp.Dummy("jet_scale", positive=True)
        subs = {self.jet(i, s): t**s * self.jet(i, s)
                for i, s in self.active_jets(expr) if s}
        return is_zero(expr.xreplace(subs) - t**weight * expr)
