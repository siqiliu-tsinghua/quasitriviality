"""Truncated vector characteristics in the physical epsilon convention."""
from dataclasses import dataclass
from typing import Iterable

from sympy import Poly, PolynomialError, S, Symbol

from .jets import _expr

from .jets import canonical, is_zero


@dataclass(frozen=True, init=False)
class Series:
    """A vector power series with coefficient N and higher unknown."""

    coefficients: tuple
    cutoff: int
    dimension: int

    def __init__(self, coefficients: Iterable[Iterable], cutoff: int):
        if not isinstance(cutoff, int) or isinstance(cutoff, bool) or cutoff < 1:
            raise ValueError("cutoff must be a positive integer")
        cs = tuple(tuple(_expr(x) for x in c) for c in coefficients)
        if not cs or not cs[0] or len(cs) > cutoff:
            raise ValueError("supply a nonempty vector series within the cutoff")
        dim = len(cs[0])
        if any(len(c) != dim for c in cs):
            raise ValueError("all coefficients must have the same dimension")
        cs += ((S.Zero,) * dim,) * (cutoff - len(cs))
        object.__setattr__(self, "coefficients", cs)
        object.__setattr__(self, "cutoff", cutoff)
        object.__setattr__(self, "dimension", dim)

    @classmethod
    def from_expr(cls, expressions, epsilon):
        """Read component expressions with an explicit common O(epsilon**N)."""
        if not isinstance(epsilon, Symbol):
            raise TypeError("epsilon must be a SymPy symbol")
        if not isinstance(expressions, (tuple, list)):
            expressions = (expressions,)
        if not expressions:
            raise ValueError("provide at least one component")
        columns = []
        cutoff = None
        for value in expressions:
            value = _expr(value)
            tail = value.getO()
            if tail is None or tail.variables != (epsilon,) or tail.point != (0,):
                raise ValueError("every component needs O(epsilon**N) about zero")
            n = tail.getn()
            if n.is_Integer is not True or n < 1:
                raise ValueError("the first unknown order must be a positive integer")
            n = int(n)
            if cutoff is not None and cutoff != n:
                raise ValueError("all components must share the same cutoff")
            cutoff = n
            try:
                polynomial = Poly(value.removeO(), epsilon)
            except PolynomialError as exc:
                raise ValueError("coefficients must form a nonnegative epsilon polynomial") from exc
            if polynomial.degree() >= cutoff:
                raise ValueError("coefficient lies beyond the explicit cutoff")
            columns.append(tuple(polynomial.nth(k) for k in range(cutoff)))
        return cls(zip(*columns), cutoff)

    def coefficient(self, k: int) -> tuple:
        if not 0 <= k < self.cutoff:
            raise IndexError(f"coefficient {k} is outside the known range")
        return self.coefficients[k]

    def map(self, operation):
        return Series(
            (tuple(operation(x) for x in c) for c in self.coefficients), self.cutoff
        )

    def as_expr(self, epsilon) -> tuple:
        """Return one SymPy expression per component with the same O term."""
        from sympy import Order

        tail = Order(epsilon**self.cutoff, (epsilon, 0))
        return tuple(
            sum(c[i] * epsilon**k for k, c in enumerate(self.coefficients)) + tail
            for i in range(self.dimension)
        )


def _compatible(space, *series):
    if any(s.dimension != space.dimension for s in series):
        raise ValueError("series dimension does not match the jet space")


def series_bracket(space, left: Series, right: Series, *, cutoff=None) -> Series:
    """Lie bracket with no use of coefficients beyond either input cutoff."""
    _compatible(space, left, right)
    known = min(left.cutoff, right.cutoff)
    if cutoff is None:
        cutoff = known
    if not 1 <= cutoff <= known:
        raise ValueError("cannot request bracket coefficients beyond the known range")
    cs = []
    zero = (S.Zero,) * space.dimension
    active_left = [not all(is_zero(x) for x in c) for c in left.coefficients[:cutoff]]
    active_right = [not all(is_zero(x) for x in c) for c in right.coefficients[:cutoff]]
    for n in range(cutoff):
        current = zero
        for i in range(n + 1):
            j = n - i
            if active_left[i] and active_right[j]:
                term = space.bracket(left.coefficient(i), right.coefficient(j))
                current = tuple(a + b for a, b in zip(current, term))
        cs.append(tuple(canonical(x) for x in current))
    return Series(cs, cutoff)


def _series_action(space, generator, target, action) -> Series:
    _compatible(space, generator, target)
    if not all(is_zero(x) for x in generator.coefficient(0)):
        raise ValueError("the generator must be O(epsilon)")
    n = min(generator.cutoff, target.cutoff)
    lead = next(
        (k for k in range(1, n) if not all(is_zero(x) for x in generator.coefficient(k))),
        n,
    )
    accum = [list(c) for c in target.coefficients[:n]]
    current = tuple(tuple(c) for c in accum)
    for j in range(1, (n - 1) // lead + 1):
        terms = []
        for k in range(n):
            value = [S.Zero] * space.dimension
            for i in range(lead, k + 1):
                q = generator.coefficient(i)
                p = current[k - i]
                if all(is_zero(x) for x in q) or all(is_zero(x) for x in p):
                    continue
                value = [a + b / j for a, b in zip(value, action(q, p))]
            terms.append(tuple(canonical(x) for x in value))
        current = tuple(terms)
        accum = [[a + b for a, b in zip(v, t)] for v, t in zip(accum, terms)]
    return Series((tuple(canonical(x) for x in row) for row in accum), n)


def lie_exponential(space, generator: Series, characteristic: Series) -> Series:
    """Compute exp(ad_G) of a vector characteristic."""
    return _series_action(space, generator, characteristic, space.bracket)


def evolution_exponential(space, generator: Series, expressions: Series) -> Series:
    """Apply exp(E_G) to an epsilon series of expressions."""
    return _series_action(
        space, generator, expressions,
        lambda q, p: tuple(space.evolutionary(q, x) for x in p),
    )


def evolutionary_action(space, characteristic: Series, expressions: Series) -> Series:
    """Apply E_X to an epsilon series of expressions."""
    _compatible(space, characteristic, expressions)
    n = min(characteristic.cutoff, expressions.cutoff)
    cs = []
    for k in range(n):
        value = [S.Zero] * space.dimension
        for i in range(k+1):
            q, p = characteristic.coefficient(i), expressions.coefficient(k-i)
            if all(is_zero(x) for x in q) or all(is_zero(x) for x in p):
                continue
            value = [a + space.evolutionary(q, b) for a, b in zip(value, p)]
        cs.append(tuple(canonical(x) for x in value))
    return Series(cs, n)


def substitution(space, generator: Series, *, inverse=False) -> Series:
    """Coordinates exp(E_G)u (or exp(-E_G)u) with the generator's cutoff."""
    g = generator.map(lambda x: -x) if inverse else generator
    return evolution_exponential(space, g, Series((space.states,), g.cutoff))
