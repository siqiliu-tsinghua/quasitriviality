"""Constructive local primitives of the inhomogeneous Tsarev equations."""

from __future__ import annotations

from hashlib import sha256
from itertools import permutations
from types import MappingProxyType

import sympy as sp
from sympy.core.function import ArgumentIndexError

from .jets import canonical, is_zero


class DarbouxUnresolved(ValueError):
    """The required exact local Darboux primitive could not be certified."""


def _subst(expr, states, points):
    if tuple(points) == tuple(states):
        return expr
    # Simultaneous substitution also handles derivatives of smooth
    # coefficient functions when a radial argument is not a symbol.
    return expr.subs(dict(zip(states, points)), simultaneous=True)


def _quadrature(expr, variable, lower):
    """Exact definite quadrature; do not guess elementary branches."""
    if is_zero(expr):
        return sp.S.Zero
    if not expr.has(variable):
        return canonical(expr * (variable - lower))
    try:
        poly = sp.Poly(expr, variable)
    except sp.PolynomialError:
        t = sp.Dummy("darboux_t", real=True)
        return sp.Integral(expr.subs(variable, t), (t, lower, variable))
    return canonical(sum(
        coeff * (variable**(k + 1) - lower**(k + 1)) / (k + 1)
        for (k,), coeff in poly.terms()
    ))


class DarbouxProblem:
    """Certified local zero-axis Volterra problem on a small state cube.

    Under the reference theorem's hypotheses (smooth coefficients and
    forcing on a cube with sufficiently small radius), h is the unique
    fixed point of h = A + T(h), equivalently the convergent Neumann sum
    of T**m(A) for m >= 0. The components may be retained as exact
    named coefficient functions, without choosing a finite truncation.
    """

    def __init__(self, chart, forcing, *, basepoint=None):
        self.chart = chart
        self.states = chart.states
        self.dimension = chart.dimension
        n = self.dimension
        self._pairs = tuple(permutations(range(n), 2))
        if any(key not in self._pairs for key in forcing):
            raise ValueError("forcing keys must be valid off-diagonal pairs (i,j)")
        self.forcing = MappingProxyType({key: canonical(forcing.get(key, 0)) for key in self._pairs})
        self.gamma = MappingProxyType({key: chart.tsarev(*key) for key in self._pairs})
        if any(
            value.has(sp.Float)
            or chart.space.max_order(value)
            or chart.active_symbols(value)
            for value in self.forcing.values()
        ):
            raise DarbouxUnresolved("Darboux forcing must consist of exact state coefficients")
        if basepoint is None:
            basepoint = tuple(sp.Dummy(f"darboux_a{i}", real=True) for i in range(n))
        self.basepoint = tuple(map(canonical, basepoint))
        if len(self.basepoint) != n or any(
            value.free_symbols.intersection(self.states) or value.has(sp.Float)
            for value in self.basepoint
        ):
            raise ValueError("basepoint must have one exact state-independent entry per component")
        self._check_basepoint()
        if any(not is_zero(v) for v in chart.semi_hamiltonian_residual().values()):
            raise DarbouxUnresolved("could not verify homogeneous Tsarev compatibility")
        f, g, x = self.forcing, self.gamma, self.states
        self.compatibility = {}
        for i, j, k in permutations(range(n), 3):
            coefficient_residual = canonical(
                g[i, j].diff(x[k]) - g[i, j] * g[j, k]
                - g[i, k] * g[k, j] + g[i, j] * g[i, k]
            )
            if not is_zero(coefficient_residual):
                raise DarbouxUnresolved("could not verify the homogeneous Darboux identity")
            if j > k:
                continue
            self.compatibility[i, j, k] = canonical(
                f[i, j].diff(x[k]) + g[i, k] * (f[i, j] - f[k, j])
                - f[i, k].diff(x[j]) - g[i, j] * (f[i, k] - f[j, k])
            )
        if any(not is_zero(v) for v in self.compatibility.values()):
            raise DarbouxUnresolved("could not verify inhomogeneous Darboux compatibility")
        self._component_types = None

    def _check_basepoint(self):
        speeds = self.chart.speeds
        n = self.dimension
        values = {
            **{f"speed {i}": speeds[i] for i in range(n)},
            **{f"speed gap ({i},{j})": canonical(speeds[j] - speeds[i])
               for i in range(n) for j in range(i + 1, n)},
            **{f"self derivative {i}": speeds[i].diff(self.states[i]) for i in range(n)},
            **{f"Gamma{key}": value for key, value in self.gamma.items()},
            **{f"forcing{key}": value for key, value in self.forcing.items()},
        }
        for label, expr in values.items():
            value = _subst(expr, self.states, self.basepoint)
            if (value.has(sp.nan, sp.zoo, sp.oo, -sp.oo)
                    or value.is_finite is False or value.is_real is False):
                raise DarbouxUnresolved(f"singular or non-real Darboux basepoint: {label}")
            if (label.startswith(("speed gap", "self derivative"))
                    and is_zero(value)):
                raise DarbouxUnresolved(f"singular Darboux basepoint: {label} vanishes")

    def verify(self, vector):
        vector = tuple(map(canonical, vector))
        if len(vector) != self.dimension:
            raise DarbouxUnresolved("Darboux solver returned the wrong component count")
        if any(value.has(sp.Float) or self.chart.space.max_order(value)
               or self.chart.active_symbols(value) for value in vector):
            raise DarbouxUnresolved("Darboux primitive must consist of exact state coefficients")
        for i, j in self._pairs:
            residual = canonical(
                vector[i].diff(self.states[j])
                - self.gamma[i, j] * (vector[j] - vector[i])
                - self.forcing[i, j]
            )
            if not is_zero(residual) and not is_zero(sp.simplify(residual)):
                raise DarbouxUnresolved(f"Darboux primitive did not verify for ({i},{j})")
        return vector

    def radial(self, row, t):
        return tuple(
            x if j == row else a + t * (x - a)
            for j, (x, a) in enumerate(zip(self.states, self.basepoint))
        )

    def _row_integral(self, rows):
        t = sp.Dummy("t", real=True)
        x, a = self.states, self.basepoint
        result = [sp.S.Zero] * self.dimension
        for i, j in self._pairs:
            term = rows[i, j]
            if not is_zero(term):
                result[i] += (x[j] - a[j]) * sp.Integral(
                    _subst(term, x, self.radial(i, t)), (t, 0, 1)
                )
        return tuple(result)

    def inhomogeneity(self):
        """The radial forcing integral A, with zero initial axis data."""
        return self._row_integral(self.forcing)

    def apply_volterra(self, vector):
        """One application of the homogeneous radial integral operator T."""
        vector = tuple(vector)
        if len(vector) != self.dimension:
            raise ValueError("expected one expression per component")
        return self._row_integral({
            (i, j): self.gamma[i, j] * (vector[j] - vector[i])
            for i, j in self._pairs
        })

    def fixed_point_rhs(self):
        """Display the exact right side A + T(h) defining the components."""
        h = self.zero_axis_solution()
        return tuple(a + t for a, t in zip(self.inhomogeneity(), self.apply_volterra(h)))

    def contraction_condition(self, radius, sup_bound):
        """Sufficient common-cube bound from the smooth Lean theorem.

        Supply a radius and a joint sup bound for Gamma and its first
        state derivative on the cube (which must lie in the smooth domain).
        """
        radius, sup_bound = map(sp.sympify, (radius, sup_bound))
        return sp.And(radius > 0, radius <= 1, sup_bound >= 0,
                      8 * self.dimension * radius * sup_bound < 1)

    def zero_axis_solution(self):
        """Return components of the exact local fixed point, not approximations."""
        if self._component_types is None:
            problem = self
            n = self.dimension
            center = self.basepoint
            key = (self.states, center, tuple(self.gamma.values()), tuple(self.forcing.values()))
            tag = sha256(repr(key).encode("utf-8")).hexdigest()[:10]

            def build(index):
                def eval_at_axis(cls, *args):
                    if all(args[j] == center[j] for j in range(n) if j != index):
                        return sp.S.Zero
                    return None

                def derivative(obj, argindex=1):
                    j = argindex - 1
                    if not 0 <= j < n or j == index:
                        # The characteristic free derivative is left as a
                        # genuine SymPy derivative (including after substitution).
                        raise ArgumentIndexError(obj, argindex)
                    gamma = _subst(problem.gamma[index, j], problem.states, obj.args)
                    force = _subst(problem.forcing[index, j], problem.states, obj.args)
                    other = problem._component_types[j](*obj.args)
                    return gamma * (other - obj) + force

                return type(
                    f"Darboux_{tag}_{index}", (sp.Function,),
                    {
                        "nargs": n,
                        "eval": classmethod(eval_at_axis),
                        "fdiff": derivative,
                        "darboux_problem": problem,
                        "darboux_component": index,
                    }
                )

            self._component_types = tuple(build(i) for i in range(n))
        return self.verify(cls(*self.states) for cls in self._component_types)

    def decoupled_solution(self):
        # Integrate each closed row along a coordinate path with the own
        # coordinate held fixed, starting from its coordinate axis.
        x, a, n = self.states, self.basepoint, self.dimension
        result = []
        for i in range(n):
            row = sp.S.Zero
            for j in range(n):
                if j != i:
                    remaining = {x[k]: a[k] for k in range(j + 1, n) if k != i}
                    part = self.forcing[i, j].subs(remaining, simultaneous=True)
                    row += _quadrature(part, x[j], a[j])
            result.append(row)
        return self.verify(result)

    def difference_solution(self):
        # This EPD class satisfies Gamma_(0,1) = -Gamma_(1,0)
        # = -1/(r-s). Adding I_s = f01 - f10 gives h0 = P_r + I,
        # h1 = P_s, and ( (r-s)*P )_rs = I+(r-s)*f10.
        if self.dimension != 2:
            raise DarbouxUnresolved("inverse-difference quadrature requires two components")
        r, s = self.states
        a, b = self.basepoint
        gap = r - s
        k = -1 / gap
        if (not is_zero(self.gamma[0, 1] - k)
                or not is_zero(self.gamma[1, 0] + k)):
            raise DarbouxUnresolved("inverse-difference quadrature requires its special Tsarev coefficients")
        i = _quadrature(self.forcing[0, 1] - self.forcing[1, 0], s, b)
        numerator = _quadrature(_quadrature(i + gap * self.forcing[1, 0], r, a), s, b)
        potential = numerator / gap
        return self.verify((potential.diff(r) + i, potential.diff(s)))


def solve_darboux(chart, forcing, *, solver=None, basepoint=None):
    """Recover h with d_j h_i = Gamma_ij*(h_j-h_i)+forcing_ij.

    A custom solver, if supplied, is called as solver(chart, forcing).
    Otherwise use zero, decoupled quadrature, or the exact local zero-axis
    Volterra solution. A special inverse-difference two-component
    quadrature is available on DarbouxProblem, with a rectangular gauge.
    The basepoint must lie in a common nonsingular coordinate box.
    With no basepoint, real symbolic anchors encode that choice.
    """
    problem = DarbouxProblem(chart, forcing, basepoint=basepoint)
    if all(is_zero(v) for v in problem.forcing.values()):
        return (sp.S.Zero,) * chart.dimension
    if solver is not None:
        return problem.verify(solver(chart, problem.forcing.copy()))
    if all(is_zero(v) for v in problem.gamma.values()):
        return problem.decoupled_solution()
    return problem.zero_axis_solution()


def volterra_definitions(*expressions):
    """Find exact defining Darboux problems referenced in SymPy expressions."""
    found = {}
    for value in expressions:
        value = sp.sympify(value)
        for function in sorted(value.atoms(sp.Function), key=sp.default_sort_key):
            problem = getattr(function.func, "darboux_problem", None)
            if isinstance(problem, DarbouxProblem):
                found.setdefault(id(problem), problem)
    return tuple(found.values())
