"""Triangular generalized hodograph jet charts over Tsarev systems.

Here z[i,s] is the s-th pure auxiliary Tsarev derivative on the hodograph
section; z[i,1] = 1/u[i,1]. The two algebras share state variables but have
distinct positive-order jets. Speeds and their differences are understood
as local expressions on the supplied nonvanishing domain.
"""

from functools import lru_cache
from operator import index

import sympy as sp

from .jets import _expr, canonical, is_zero


class HodographChart:
    """A diagonal leading system and its normalized hodograph calculus.

    The transport identities assume a semihamiltonian leading system, and
    the evolutionary operator assumes its speed argument solves Tsarev's
    equations. The chart exposes the equations for checking these hypotheses.
    """

    def __init__(self, space, speeds):
        self.space = space
        self.dimension = space.dimension
        self.states = space.states
        self.speeds = self._vector(speeds)
        for i, speed in enumerate(self.speeds):
            if any(order for _, order in space.active_jets(speed)):
                raise ValueError("leading speeds must depend only on the states")
            if any(speed.has(space.log(j)) for j in range(self.dimension)):
                raise ValueError("leading speeds cannot depend on jet logarithms")
            for j in range(i):
                if is_zero(speed - self.speeds[j]):
                    raise ValueError("leading speeds must have nonzero differences")
        self._coordinate_indices = {}
        self._log_indices = {}
        self._coordinates = {}
        self._logs = {}

    def _vector(self, values):
        values = tuple(map(_expr, values))
        if len(values) != self.dimension:
            raise ValueError("component count does not match the hodograph chart")
        return values

    def _component(self, i):
        i = index(i)
        if not 0 <= i < self.dimension:
            raise IndexError("component index out of range")
        return i

    def jet(self, i, s):
        """Auxiliary inverse-jet symbol, with s starting at one."""
        i, s = self._component(i), index(s)
        if s < 1:
            raise ValueError("hodograph jet indices start at one")
        if (i, s) not in self._coordinates:
            z = sp.Dummy(f"z{i}_{s}", real=True)
            self._coordinates[i, s] = z
            self._coordinate_indices[z] = i, s
        return self._coordinates[i, s]

    def log(self, i):
        """Formal log of abs(z[i,1]); its derivative is 1/z[i,1]."""
        i = self._component(i)
        if i not in self._logs:
            ell = sp.Dummy(f"logz{i}", real=True)
            self._logs[i] = ell
            self._log_indices[ell] = i
        return self._logs[i]

    @lru_cache(maxsize=256)
    def tsarev(self, i, j):
        i, j = self._component(i), self._component(j)
        if i == j:
            raise ValueError("Tsarev coefficients require distinct indices")
        return canonical(self.speeds[i].diff(self.states[j]) / (self.speeds[j] - self.speeds[i]))

    def commuting_residual(self, speeds):
        """Return all transverse Tsarev equations, keyed by (i,j), i != j."""
        w = self._vector(speeds)
        return {
            (i, j): canonical(w[i].diff(self.states[j]) - self.tsarev(i, j) * (w[j] - w[i]))
            for i in range(self.dimension)
            for j in range(self.dimension)
            if i != j
        }

    def semi_hamiltonian_residual(self):
        """Return the Tsarev compatibility equations for distinct triples."""
        return {
            (i, j, k): canonical(self.tsarev(i, j).diff(self.states[k])
                                  - self.tsarev(i, k).diff(self.states[j]))
            for i in range(self.dimension)
            for j in range(self.dimension)
            for k in range(self.dimension)
            if len({i, j, k}) == 3
        }

    @lru_cache(maxsize=512)
    def value_coefficient(self, i, j, n):
        """Coefficient of w[j]-w[i] in d_i^n d_j w[i] (Lean ValueCoeff)."""
        if n < 0:
            raise ValueError("coefficient order cannot be negative")
        if n == 0:
            return self.tsarev(i, j)
        previous = self.value_coefficient(i, j, n - 1)
        return canonical(previous.diff(self.states[i]) - self.tsarev(j, i) * previous)

    @lru_cache(maxsize=2048)
    def jet_coefficient(self, i, j, n, t):
        """Coefficient of d_i^(t+1) w[i] in d_i^n d_j w[i]."""
        if n < 0 or t < 0:
            raise ValueError("coefficient orders cannot be negative")
        if n == 0 or t >= n:
            return sp.S.Zero
        previous = self.jet_coefficient(i, j, n - 1, t).diff(self.states[i])
        if t == 0:
            return canonical(previous - self.value_coefficient(i, j, n - 1))
        return canonical(previous + self.jet_coefficient(i, j, n - 1, t - 1))

    @lru_cache(maxsize=512)
    def inverse_jet(self, i, s):
        """Physical expression for z[i,s], triangular in u[i,s]."""
        i, s = self._component(i), index(s)
        if s < 1:
            raise ValueError("hodograph jet indices start at one")
        p = self.space.jet(i, 1)
        if s == 1:
            return 1 / p
        transverse = sum(
            self.space.jet(j, 1)
            * sum(self.jet_coefficient(i, j, s - 1, t) * self.inverse_jet(i, t + 1)
                  for t in range(s - 1))
            for j in range(self.dimension)
            if j != i
        )
        return canonical((self.space.dx(self.inverse_jet(i, s - 1)) - transverse) / p)

    @lru_cache(maxsize=512)
    def physical_jet(self, i, s):
        """Inverse triangular substitution: expression for u[i,s] in z."""
        i, s = self._component(i), index(s)
        if s < 0:
            raise ValueError("physical jet indices must be nonnegative")
        if s == 0:
            return self.states[i]
        z1 = self.jet(i, 1)
        if s == 1:
            return 1 / z1
        # The top term is -u[i,s]/p_i**(s+1); all remaining jets have
        # strictly lower order. Remove the top term before recursing.
        lower = self.inverse_jet(i, s).xreplace({self.space.jet(i, s): sp.S.Zero})
        active = tuple((j, k) for j, k in self.space.active_jets(lower) if k > 0)
        if any(k >= s for _, k in active):
            raise ArithmeticError("hodograph inverse has an unexpected top-order jet")
        converted = lower.xreplace({self.space.jet(j, k): self.physical_jet(j, k)
                                    for j, k in active})
        return canonical((converted - self.jet(i, s)) / z1 ** (s + 1))

    def to_hodograph(self, expr):
        """Substitute physical jets and log|u_x| into the auxiliary chart."""
        expr = _expr(expr)
        replacements = {self.space.jet(i, s): self.physical_jet(i, s)
                        for i, s in self.space.active_jets(expr) if s > 0}
        replacements.update({self.space.log(i): -self.log(i)
                             for i in range(self.dimension) if expr.has(self.space.log(i))})
        return canonical(expr.xreplace(replacements))

    def active_symbols(self, expr):
        """Return occurring auxiliary positive jets and logs in stable order."""
        symbols = _expr(expr).free_symbols
        coordinates = sorted(
            ((i, k), symbol)
            for symbol, (i, k) in self._coordinate_indices.items()
            if symbol in symbols
        )
        logs = sorted(
            (i, symbol) for symbol, i in self._log_indices.items() if symbol in symbols
        )
        return tuple(symbol for _, symbol in coordinates + logs)

    def to_physical(self, expr):
        """Substitute inverse jets and their formal logs back to physical jets."""
        expr = _expr(expr)
        replacements = {}
        for a in expr.free_symbols:
            if a in self._coordinate_indices:
                replacements[a] = self.inverse_jet(*self._coordinate_indices[a])
            elif a in self._log_indices:
                replacements[a] = -self.space.log(self._log_indices[a])
        return canonical(expr.xreplace(replacements))

    def normalize(self, characteristic):
        return tuple(self.to_hodograph(q / self.space.jet(i, 1))
                     for i, q in enumerate(self._vector(characteristic)))

    def denormalize(self, characteristic):
        return tuple(canonical(self.space.jet(i, 1) * self.to_physical(h))
                     for i, h in enumerate(self._vector(characteristic)))

    @lru_cache(maxsize=2048)
    def _aux_velocity(self, w, i, s):
        """Transported evolutionary action on z[i,s] (s >= 1)."""
        own = w[i] * self.jet(i, s + 1) / self.jet(i, 1)
        transverse = sum(
            w[j] / self.jet(j, 1)
            * sum(self.jet_coefficient(i, j, s, t) * self.jet(i, t + 1)
                  for t in range(s))
            for j in range(self.dimension)
            if j != i
        )
        return canonical(own + transverse - w[i].diff(self.states[i], s))

    def evolutionary(self, speeds, expr):
        """Explicit transported derivative of a function in hodograph jets."""
        w, expr = self._vector(speeds), _expr(expr)
        terms = [w[i] * expr.diff(state) / self.jet(i, 1)
                 for i, state in enumerate(self.states) if state in expr.free_symbols]
        for a in expr.free_symbols:
            if a in self._coordinate_indices:
                i, s = self._coordinate_indices[a]
                terms.append(self._aux_velocity(w, i, s) * expr.diff(a))
            elif a in self._log_indices:
                i = self._log_indices[a]
                terms.append(self._aux_velocity(w, i, 1) * expr.diff(a) / self.jet(i, 1))
        return canonical(sum(terms))

    def total(self, expr):
        """Transported total x derivative (the common unit speed)."""
        return self.evolutionary((sp.S.One,) * self.dimension, expr)

    def lie_action(self, speeds, characteristic):
        """Normalized bracket corresponding to [w[i] u[i,1], F]."""
        w = self._vector(speeds)
        h = self._vector(characteristic)
        return tuple(
            canonical(
                self.evolutionary(w, h[i])
                - w[i] * self.total(h[i])
                + sum(w[i].diff(self.states[j]) * (h[i] - h[j]) / self.jet(j, 1)
                      for j in range(self.dimension) if j != i)
            )
            for i in range(self.dimension)
        )
