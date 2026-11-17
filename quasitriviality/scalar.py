"""Joel's scalar homological inverse in the original dependent variable."""

from __future__ import annotations

from math import factorial

import sympy as sp

from .jets import JetSpace, canonical, is_zero


class ScalarPrimitiveError(ValueError):
    """The inputs or a symbolic step fall outside the verified scalar solver."""


def _terms(expr: sp.Expr):
    """Expose additive Laurent monomials after cancelling rational factors."""
    return sp.Add.make_args(sp.expand(canonical(expr), log=False))


def _split_monomial(term: sp.Expr, generators: tuple[sp.Symbol, ...]):
    powers = term.as_powers_dict()
    exponents = tuple(powers.get(g, sp.S.Zero) for g in generators)
    if any(exponent.is_Integer is not True for exponent in exponents):
        raise ScalarPrimitiveError("noninteger power of a jet or formal logarithm")
    coefficient = term / sp.prod(g**n for g, n in zip(generators, exponents))
    if coefficient.has(*generators):
        raise ScalarPrimitiveError("input has non-Laurent–log jet dependence")
    return coefficient, tuple(map(int, exponents))


def laurent_log_primitive(expr: sp.Expr, p: sp.Symbol, ell: sp.Symbol) -> sp.Expr:
    """Integrate a finite Laurent polynomial in p and log|p|, with zero constant.

    The relation used is d(ell)/dp = 1/p. Other symbols are constants
    for this operation. It deliberately does not invoke indefinite integration
    or branch-dependent simplification.
    """
    if p == ell:
        raise ScalarPrimitiveError("the jet and its formal logarithm must be distinct")
    result = []
    for term in _terms(expr):
        c, (k, m) = _split_monomial(term, (p, ell))
        if m < 0:
            raise ScalarPrimitiveError("negative powers of a formal logarithm")
        if k == -1:
            result.append(c * ell**(m + 1) / (m + 1))
        else:
            a = sp.Integer(k + 1)
            correction = sum(
                (-1)**j * (factorial(m) // factorial(m - j)) * ell**(m - j)
                / a**(j + 1)
                for j in range(m + 1)
            )
            result.append(c * p**a * correction)
    return sp.Add(*result)


def _inverse_jet_maps(space: JetSpace, speed: sp.Expr, order: int):
    """Produce simultaneous maps to/from invariants of E_(f p) - f D_x."""
    u, p = space.states[0], space.jet(0, 1)
    fprime = sp.diff(speed, u)

    # The adapted inverse jets are b_s = z_s - (f^(s)/f') z_1, s >= 2,
    # where z_s = d^s x/du^s and z_1 = 1/p. In particular, b_2/f'^2
    # is Joel's GGD[2]. Dummies avoid collisions with input parameters.
    z = {s: sp.Dummy(f"z{s}") for s in range(1, order + 1)}
    b = {s: sp.Dummy(f"b{s}") for s in range(2, order + 1)}
    ratios = {s: canonical(sp.diff(speed, u, s) / fprime) for s in b}
    z_to_invariants = {z[1]: 1 / p}
    z_to_invariants.update({z[s]: b[s] + ratios[s] / p for s in b})

    jet_of_z = 1 / z[1]
    p_to_inverse = {}
    z_of_p = 1 / p
    inverse_to_p = {}
    for s in range(2, order + 1):
        jet_of_z = sp.expand(sum(
            z[j + 1] * sp.diff(jet_of_z, z[j]) for j in range(1, s)
        ) / z[1])
        p_to_inverse[space.jet(0, s)] = jet_of_z.xreplace(z_to_invariants)

        z_of_p = canonical(space.dx(z_of_p) / p)
        inverse_to_p[b[s]] = z_of_p - ratios[s] / p
    return p_to_inverse, inverse_to_p


def scalar_primitive(space: JetSpace, speed: sp.Expr, source: sp.Expr) -> sp.Expr:
    """Return a verified Laurent–log F solving [speed(u)*u_x, F] = source.

    This is local on a domain where speed'(u) and u_x are nonzero. The symbolic
    condition speed' != 0 is assumed on that domain; an identically zero
    derivative is rejected. A zero primitive constant is chosen in the
    invariant coordinates. The method never requires a Hopf point change.
    """
    if space.dimension != 1:
        raise ScalarPrimitiveError("Joel's primitive is a scalar solver")
    speed, source = canonical(speed), canonical(source)
    if speed.atoms(sp.Float) or source.atoms(sp.Float):
        raise ScalarPrimitiveError("inexact floats are outside the exact solver")
    u, p, ell = space.states[0], space.jet(0, 1), space.log(0)
    if speed.has(ell) or any(n > 0 for _, n in space.active_jets(speed)):
        raise ScalarPrimitiveError("the speed must depend only on the state")
    fprime = canonical(sp.diff(speed, u))
    if is_zero(fprime):
        raise ScalarPrimitiveError("the speed derivative vanishes identically")

    order = max((n for _, n in space.active_jets(source)), default=1)
    order = max(order, 1)
    physical_generators = tuple(space.jet(0, s) for s in range(1, order + 1)) + (ell,)
    for term in _terms(source):
        _, exponents = _split_monomial(term, physical_generators)
        if any(n < 0 for n in exponents[1:]):
            raise ScalarPrimitiveError("only the first physical jet is invertible")

    to_inverse, to_physical = _inverse_jet_maps(space, speed, order)
    source_in_invariants = source.xreplace(to_inverse)
    primitive_in_invariants = p * laurent_log_primitive(
        source_in_invariants / p**3, p, ell
    ) / fprime
    result = canonical(primitive_in_invariants.xreplace(to_physical))
    if not is_zero(space.bracket((speed * p,), (result,))[0] - source):
        raise ScalarPrimitiveError("could not verify the primitive in physical jets")
    return result
