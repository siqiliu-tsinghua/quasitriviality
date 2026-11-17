"""Compute and verify the first-class quasi-Miura transformation."""
from dataclasses import dataclass

from sympy import Float, Function, S

from .jets import JetSpace, canonical, is_zero
from .scalar import ScalarPrimitiveError, scalar_primitive
from .series import Series, evolutionary_action, evolution_exponential, lie_exponential, substitution


class UnresolvedTransformation(ValueError):
    """The configured exact solver could not produce a verified witness."""

    def __init__(self, order, residual, message):
        self.order = order
        self.residual = residual
        super().__init__(f"order {order}: {message}")


def _series_equal(left, right):
    return (
        left.cutoff == right.cutoff and left.dimension == right.dimension
        and all(is_zero(a-b)
                for l, r in zip(left.coefficients, right.coefficients)
                for a, b in zip(l, r))
    )


def _function_names(expressions):
    return {str(fn.func) for expr in expressions for fn in expr.atoms(Function)}


def _fresh_function_names(stems, used):
    occupied = set(used)
    names = []
    for name in stems:
        if not isinstance(name, str) or not name:
            raise ValueError("symbolic speed names must be nonempty strings")
        while name in occupied:
            name += "_"
        names.append(name)
        occupied.add(name)
    return tuple(names)


def _leading_speeds(space, equation):
    if equation.dimension != space.dimension:
        raise ValueError("equation dimension does not match the jet space")
    speeds = []
    for i, coefficient in enumerate(equation.coefficient(0)):
        p = space.jet(i, 1)
        v = canonical(coefficient/p)
        if any(s > 0 for _, s in space.active_jets(v)):
            raise ValueError("leading term must be diagonal hydrodynamic in the input coordinates")
        if is_zero(v.diff(space.states[i])):
            raise ValueError("the leading speed must have a nonzero self-derivative")
        speeds.append(v)
    if space.dimension > 1:
        from .tsarev import TsarevFamily

        names = _fresh_function_names(
            (f"w_{i}" for i in range(space.dimension)), _function_names(speeds)
        )
        TsarevFamily(space, speeds, names=names)
    for k, coefficient in enumerate(equation.coefficients):
        for x in coefficient:
            if x.has(Float) or not space.is_homogeneous(x, k+1):
                raise ValueError(f"flow coefficient {k} needs exact differential weight {k+1}")
    return tuple(speeds)


@dataclass(frozen=True)
class QuasiMiura:
    """Verified witness; all series use the input's first unknown epsilon order."""

    space: JetSpace
    equation: Series
    speeds: tuple
    generator: Series
    reducing: Series
    inverse: Series

    @property
    def cutoff(self):
        return self.equation.cutoff


def quasitrivialize(space, equation: Series, *, solver=None) -> QuasiMiura:
    """Find G with exp(ad_G) X0 = X and return G, T and T inverse.

    The domain must have nonvanishing speed gaps and self-derivatives.
    Coupled input is integrated by triangular elimination in generalized
    hodograph jets, with an exact local Darboux integral at weight one.
    An unresolved computation is not an integrability obstruction.
    """
    speeds = _leading_speeds(space, equation)
    equation = equation.map(canonical)
    n = equation.cutoff
    x0 = Series((equation.coefficient(0),), n)
    if solver is None and space.dimension > 1 and n > 1:
        from .coupled import CoupledSolver

        solver = CoupledSolver()
    zeros = (S.Zero,) * space.dimension
    cs = [zeros]
    for k in range(1, n):
        if all(is_zero(x) for vector in cs for x in vector):
            approximate = zeros
        else:
            g_partial = Series(cs, k+1)
            approximate = lie_exponential(
                space, g_partial, Series((x0.coefficient(0),), k+1)
            ).coefficient(k)
        residual = tuple(canonical(x-y) for x, y in zip(equation.coefficient(k), approximate))
        source = tuple(-x for x in residual)
        try:
            if all(is_zero(y) for y in residual):
                gk = zeros
            elif solver is not None:
                gk = tuple(solver(space, speeds, source))
            else:
                gk = (scalar_primitive(space, speeds[0], source[0]),)
        except (ScalarPrimitiveError, ValueError, NotImplementedError) as exc:
            raise UnresolvedTransformation(k, residual, str(exc)) from exc
        gk = tuple(canonical(x) for x in gk)
        if (len(gk) != space.dimension
                or any(x.has(Float) or not space.is_homogeneous(x, k) for x in gk)
                or any(not is_zero(a-b) for a, b in zip(
                    space.bracket(x0.coefficient(0), gk), source
                ))):
            raise UnresolvedTransformation(k, residual, "solver witness failed exact verification")
        cs.append(gk)
    g = Series(cs, n)
    if not _series_equal(lie_exponential(space, g, x0), equation):
        raise ArithmeticError("Lie reconstruction failed")
    reducing = substitution(space, g)
    inverse = substitution(space, g, inverse=True)
    if not _series_equal(
        evolutionary_action(space, equation, reducing), evolution_exponential(space, g, x0)
    ):
        raise ArithmeticError("reducing substitution failed the chain rule verification")
    identity = Series((space.states,), n)
    minus_g = g.map(lambda x: -x)
    if not (_series_equal(evolution_exponential(space, g, inverse), identity)
            and _series_equal(evolution_exponential(space, minus_g, reducing), identity)):
        raise ArithmeticError("reducing and inverse substitutions failed composition")
    return QuasiMiura(space, equation, speeds, g, reducing, inverse)
