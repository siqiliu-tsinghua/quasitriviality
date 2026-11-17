"""Generic commuting speeds on a local semi-Hamiltonian Riemann chart."""

from __future__ import annotations

from itertools import permutations
from typing import Iterable

import sympy as sp
from sympy.core.function import AppliedUndef

from .jets import JetSpace, canonical, is_zero


class TsarevFamily:
    """Represent arbitrary commuting speeds with independent pure derivatives.

    For each component i, keep w_i and all pure state derivatives in u^i.
    The transverse rule is

        d_j w_i = Gamma_(i,j) * (w_j - w_i),    j != i.

    Every rewrite of a derivative of order m uses derivatives of order at
    most m-1. Exact semi-Hamiltonian identities ensure the result does not
    depend on the order of transverse differentiations. Function symbols
    represent a local family on the domain where the speed differences are
    nonzero; this class does not construct global speed functions.
    """

    def __init__(
        self,
        space: JetSpace,
        speeds: Iterable[sp.Expr],
        names: Iterable[str] | None = None,
    ) -> None:
        self.space = space
        self.speeds = tuple(map(canonical, speeds))
        if len(self.speeds) != space.dimension:
            raise ValueError("expected one leading speed per component")
        if any(
            any(order for _, order in space.active_jets(speed))
            for speed in self.speeds
        ):
            raise ValueError("leading speeds must not depend on positive physical jets")
        self._gamma = {}
        for i, j in permutations(range(space.dimension), 2):
            gap = canonical(self.speeds[j] - self.speeds[i])
            if is_zero(gap):
                raise ValueError("leading speeds must be distinct on the local domain")
            self._gamma[i, j] = canonical(self.speeds[i].diff(space.states[j]) / gap)
        for i, j, k in permutations(range(space.dimension), 3):
            if j < k and not is_zero(
                self.gamma(i, j).diff(space.states[k])
                - self.gamma(i, k).diff(space.states[j])
            ):
                raise ValueError(
                    "could not verify the semi-Hamiltonian identity "
                    f"for components {(i, j, k)}"
                )
        if names is None:
            names = (f"w_{i}" for i in range(space.dimension))
        names = tuple(names)
        if (
            len(names) != space.dimension
            or any(not isinstance(name, str) or not name for name in names)
            or len(set(names)) != len(names)
        ):
            raise ValueError("provide one distinct nonempty function name per component")
        existing = {str(term.func) for speed in self.speeds for term in speed.atoms(AppliedUndef)}
        if set(names) & existing:
            raise ValueError("generic speed names conflict with leading-speed functions")
        self.values = tuple(sp.Function(name)(*space.states) for name in names)
        self._components = {speed: i for i, speed in enumerate(self.values)}
        self._state_indices = {state: i for i, state in enumerate(space.states)}
        self._rules: dict[sp.Derivative, sp.Expr] = {}

    def gamma(self, i: int, j: int) -> sp.Expr:
        """Return the off-diagonal Tsarev coefficient (zero on the diagonal)."""
        if not (0 <= i < self.space.dimension and 0 <= j < self.space.dimension):
            raise IndexError("component index out of range")
        return sp.S.Zero if i == j else self._gamma[i, j]

    def free_jets(self, max_order: int) -> tuple[sp.Expr, ...]:
        """List the independent speed values and pure derivatives, by component."""
        if not isinstance(max_order, int) or max_order < 0:
            raise ValueError("max_order must be a nonnegative integer")
        return tuple(
            w.diff(self.space.states[i], k)
            for i, w in enumerate(self.values)
            for k in range(max_order + 1)
        )

    def _rule(self, derivative: sp.Derivative) -> sp.Expr:
        if derivative in self._rules:
            return self._rules[derivative]
        i = self._components.get(derivative.expr)
        if i is None:
            return derivative
        variables = list(derivative.variables)
        # Non-state derivative arguments are outside this local speed model.
        if any(var not in self._state_indices for var in variables):
            return derivative
        transverse = [
            (self._state_indices[var], pos)
            for pos, var in enumerate(variables)
            if self._state_indices[var] != i
        ]
        if not transverse:
            return derivative
        j, pos = min(transverse)
        del variables[pos]
        rhs = self.gamma(i, j) * (self.values[j] - self.values[i])
        rhs = rhs.diff(*variables) if variables else rhs
        result = self._rewrite(rhs)
        self._rules[derivative] = result
        return result

    def _rewrite(self, expr: sp.Expr) -> sp.Expr:
        replacements = {}
        for derivative in expr.atoms(sp.Derivative):
            replacement = self._rule(derivative)
            if replacement != derivative:
                replacements[derivative] = replacement
        return expr.xreplace(replacements)

    def reduce(self, expr: sp.Expr) -> sp.Expr:
        """Rewrite direct derivatives of these speed symbols in free speed jets.

        Unrelated undefined functions or unevaluated derivatives stay intact.
        For a formally differentiated Tsarev expression, evaluate the state
        differentiations normally (e.g. with SymPy diff) before reducing.
        """
        return canonical(self._rewrite(sp.sympify(expr)))
