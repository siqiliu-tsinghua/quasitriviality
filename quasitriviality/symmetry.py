"""Transport dispersionless symmetries by a previously verified transform."""
from dataclasses import dataclass

from sympy import Derivative, Float, Function

from .jets import canonical, is_zero
from .polynomiality import PolynomialityReport, analyze_polynomiality
from .series import Series, lie_exponential, series_bracket
from .transform import QuasiMiura, quasitrivialize


@dataclass(frozen=True)
class SymmetryResult:
    symmetry: Series
    leading_speeds: tuple
    polynomiality: PolynomialityReport
    universal_family: bool

    @property
    def cutoff(self):
        return self.symmetry.cutoff


@dataclass(frozen=True)
class IntegrabilityResult:
    transformation: QuasiMiura
    symmetry: SymmetryResult


def transport_symmetry(
    transform: QuasiMiura,
    leading_speeds=None,
    *,
    family=None,
    scalar_name="h",
) -> SymmetryResult:
    """Transport and check arbitrary leaders, or a selected admissible speed.

    If omitted, a symbolic generic family is used. A concrete leading
    speed reports polynomiality only for that speed. An explicit
    TsarevFamily may choose the symbols used in a universal calculation.
    """
    if family is not None and leading_speeds is not None:
        raise ValueError("choose an explicit speed or a generic family, not both")
    space = transform.space
    universal = leading_speeds is None
    reducer = canonical
    if family is not None or (universal and space.dimension > 1):
        from .tsarev import TsarevFamily

        if family is None:
            family = TsarevFamily(space, transform.speeds)
        if (not isinstance(family, TsarevFamily) or family.space is not space
                or family.speeds != transform.speeds):
            raise ValueError("Tsarev family belongs to a different leading system")
        leading_speeds = family.values
        reducer = family.reduce
    elif universal:
        name = scalar_name
        used = {str(fn.func) for c in transform.equation.coefficients for x in c
                for fn in x.atoms(Function)}
        while name in used:
            name += "_"
        leading_speeds = (Function(name)(space.states[0]),)
    elif space.dimension == 1 and not isinstance(leading_speeds, (tuple, list)):
        leading_speeds = (leading_speeds,)

    speeds = tuple(canonical(w) for w in leading_speeds)
    if len(speeds) != space.dimension or any(
        w.has(Float) or any(s for _, s in space.active_jets(w)) for w in speeds
    ):
        raise ValueError("leader needs an exact state-dependent speed per component")
    leader = tuple(w*space.jet(i, 1) for i, w in enumerate(speeds))
    x0 = transform.equation.coefficient(0)
    if any(not is_zero(reducer(x)) for x in space.bracket(x0, leader)):
        raise ValueError("the supplied leading symmetry does not commute with X0")
    seed = Series((leader,), transform.cutoff)
    symmetry = lie_exponential(space, transform.generator, seed).map(reducer)
    if any(not is_zero(reducer(a-b)) for a,b in zip(symmetry.coefficient(0), leader)):
        raise ArithmeticError("transport changed the leading term")
    bracket = series_bracket(space, transform.equation, symmetry)
    if any(not is_zero(reducer(x)) for c in bracket.coefficients for x in c):
        raise ArithmeticError("could not verify commutation in physical jets")

    independent = ()
    if universal:
        active_orders = (len(d.variables)
                         for c in symmetry.coefficients for x in c
                         for d in x.atoms(Derivative) if d.expr in speeds)
        n = max(active_orders, default=0)
        if family is not None:
            independent = family.free_jets(n)
        else:
            independent = tuple(speeds[0].diff(space.states[0], k) for k in range(n+1))
    result = analyze_polynomiality(
        space, symmetry, independent=independent, coefficient_reducer=reducer
    )
    return SymmetryResult(symmetry, speeds, result, universal)


def check_integrability(space, equation, leading_speeds=None, *, solver=None):
    """Run the transform-first finite-order test; never assert all-order integrability."""
    transform = quasitrivialize(space, equation, solver=solver)
    return IntegrabilityResult(transform, transport_symmetry(transform, leading_speeds))
