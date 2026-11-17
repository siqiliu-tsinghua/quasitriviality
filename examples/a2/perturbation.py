"""Transform-first residuals for an A₂ flow with one epsilon-square correction.

The three coefficient functions depend on U alone; their U-derivatives
are polynomial tags over the Toda example's exact sparse Laurent algebra.
All higher input corrections are zero. The routine exposes the **physical
bracket residual** of a deterministic triangular primitive, not a fitted
candidate vector. See the example notes for the independent polynomial
symmetry obstruction that rules out nonzero residual branches.
"""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from types import FunctionType, SimpleNamespace

import sympy as sp

from examples.a2 import compute as a2
from examples.toda import compute as toda


U = sp.Symbol("U", positive=True)
BROOT = sp.Symbol("b", positive=True)
FUNCTIONS = tuple(sp.Function(f"a{i}")(U) for i in (1, 2, 3))


def _copy_function(source, namespace):
    target = FunctionType(
        source.__code__, namespace, source.__name__, source.__defaults__, source.__closure__
    )
    target.__kwdefaults__ = source.__kwdefaults__
    return target


def _private_namespace(module, replacements=None):
    """Clone only locally declared functions, reusing algorithms without globals mutation."""
    scope = {}
    replacements = replacements or {}
    copies = {}
    for name, item in vars(module).items():
        wrapped = getattr(item, "__wrapped__", item)
        if (
            isinstance(wrapped, FunctionType)
            and wrapped.__globals__ is vars(module)
        ):
            local = _copy_function(wrapped, scope)
            if wrapped is not item:
                local = lru_cache(**item.cache_parameters())(local)
            scope[name] = copies[id(item)] = local
        else:
            scope[name] = replacements.get(id(item), item)
    return scope, copies


class _FormalCoefficients:
    """Dedicated Laurent algebra/chart for parameters and their U-derivatives."""

    def __init__(self, *, constant):
        self.constant = constant
        # e2 uses only the first few coefficient derivatives; e4 constants
        # never differentiate. A guard prevents silent omission at higher orders.
        self.max_derivative = 8
        self.tag_start = toda.VARIABLE_COUNT
        self.tag_count = 3 * (self.max_derivative + 1)
        self.tags = tuple(sp.Dummy(f"a{i+1}_{k}") for i in range(3)
                          for k in range(self.max_derivative + 1))
        self.tag_values = tuple(f.diff(U, k) for f in FUNCTIONS
                                for k in range(self.max_derivative + 1))
        t, copies = _private_namespace(toda)
        t["VARIABLE_COUNT"] += self.tag_count
        t["ZERO_KEY"] = (0,) * t["VARIABLE_COUNT"]
        # The inherited __class__ closure in Laurent.__init__ still points
        # to Toda's Laurent base, so its super() correctly resolves to dict.
        methods = {name: _copy_function(func, t)
                   for name, func in vars(toda.L).items() if isinstance(func, FunctionType)}
        t["L"] = type("_TaggedLaurent", (toda.L,), methods)
        copies[id(toda.L)] = t["L"]
        c, _ = _private_namespace(a2, copies)
        self.t = SimpleNamespace(**t)
        self.c = SimpleNamespace(**c)
        self.L = self.t.L
        self.b = self.t.mono(0)
        c["b"] = self.b
        old_hodo, old_phys = c["to_hodo_var"], c["to_physical_var"]
        c["to_hodo_var"] = lambda i: self.t.mono(i) if i >= self.tag_start else old_hodo(i)
        c["to_physical_var"] = lambda i: self.t.mono(i) if i >= self.tag_start else old_phys(i)
        c.update(dstate=self.dstate, dx=self.dx, evol=self.evol, firstpart=self.firstpart)

    def parameter(self, family, derivative=0):
        if not 0 <= derivative <= self.max_derivative:
            raise ValueError("coefficient derivative exceeds the configured jet bound")
        k = self.tag_start + family * (self.max_derivative + 1) + derivative
        return self.t.mono(k)

    def coefficient_derivative(self, expr):
        result = self.L()
        if self.constant:
            return result
        for index in expr.active():
            if index >= self.tag_start:
                family, order = divmod(index-self.tag_start, self.max_derivative+1)
                result += expr.partial(index)*self.parameter(family, order+1)
        return result

    def dstate(self, expr, i):
        # r-s=12*b**3; dU/dr=1/(6b), dU/ds=-1/(6b).
        difference = expr.partial(0)/(36*self.b**2)
        difference += self.coefficient_derivative(expr)/(6*self.b)
        return difference*(1 if i == 0 else -1) + expr.partial(self.t.V_INDEX)/2

    def dx(self, expr, times=1):
        t = self.t
        for _ in range(times):
            result = sum((self.dstate(expr, i)*t.j(i,1) for i in range(2)), self.L())
            for i in range(2):
                for order in range(1, expr.degree_order()+1):
                    term = expr.partial(t.variable_index(i, order))
                    if term:
                        result += term*t.j(i, order+1)
                term = expr.partial(t.log_index(i))
                if term:
                    result += term*t.j(i,2)/t.j(i,1)
            expr = result
        return expr

    def evol(self, vector, target):
        t = self.t
        maxorder = max((a.degree_order() for a in target), default=0)
        prolonged = [tuple(self.dx(a, k) for k in range(maxorder+1)) for a in vector]
        result = []
        for a in target:
            value = sum((self.dstate(a, i)*vector[i] for i in range(2)), self.L())
            for i in range(2):
                for order in range(1, a.degree_order()+1):
                    if term := a.partial(t.variable_index(i,order)):
                        value += term*prolonged[i][order]
                if term := a.partial(t.log_index(i)):
                    value += term*prolonged[i][1]/t.j(i,1)
            result.append(value)
        return tuple(result)

    def firstpart(self, row, rhs, weight):
        """Known own first-jet coefficient, leaving coefficient tags untouched."""
        t = self.t
        index = t.variable_index(row,1)
        ell = t.log_index(row)
        slope = self.dstate(self.c.vspeed(row),row)
        result = defaultdict(Q)
        for key, value in rhs.items():
            if key[index] != -weight-1:
                continue
            degree = key[ell]
            base = t.move_power(t.move_power(key,index,weight+1),ell,-degree)
            if any(base[1:self.tag_start]):
                continue
            for k in range(degree+1):
                new = t.move_power(t.move_power(base,index,-weight),ell,degree-k)
                result[new] += value*Q(factorial(degree),factorial(degree-k)*weight**(k+1))
        return self.L(result)/slope

    def _jet_filter(self, expr, order):
        # Physical and auxiliary higher-jet slots have the same indexing.
        return self.L({k:v for k,v in expr.items() if not any(k[1+2*order:1+2*toda.N])})

    def _transverse(self, row, rhs, order):
        t, c = self.t, self.c
        other = 1-row
        f = rhs.partial(t.variable_index(other,order+1))*t.j(other,1)
        f /= c.vspeed(other)-c.vspeed(row)
        return t.primitive(self._jet_filter(f,order),other,order)

    def _transport_projection(self, row, rhs, weight, order, factor):
        """Triangular candidate; compare the full original equation afterwards."""
        if not rhs or order < 1:
            return self.L()
        t, c = self.t, self.c
        top = self._transverse(row,rhs,order)
        rest = self._jet_filter(rhs-c.transport(row,top,factor),order)
        if order == 1:
            return top+self.firstpart(row,rest,weight)
        ind = t.variable_index(row,order)
        for k in range(rest.max_exponent(ind),0,-1):
            coeff = rest.coef(ind,k)
            term = t.j(row,order)**k*self._transport_projection(
                row,coeff,weight+k,order-1,factor-k)
            rest -= c.transport(row,term,factor)
            top += term
        return top+self._transport_projection(row,rest.coef(ind,0),weight,order-1,factor)

    def homological_projection(self, rhs, weight):
        """Recover a deterministic candidate and compute its unprojected bracket residual.

        At an inconsistent intermediate step, project away uncancelled top
        terms to continue the recurrence. The final residual is recomputed
        against every monomial of the original source in physical jets.
        """
        if weight < 2:
            raise ValueError("this screening starts with even generator orders >=2")
        t, c, m = self.t, self.c, weight-1
        rest = c.normalize(rhs)
        answer = self.L(),self.L()

        def subtract(parts):
            nonlocal rest, answer
            rest = t.plus(rest,t.mul(c.hodo_action(parts),-1))
            answer = t.plus(answer,parts)

        for n in range(max(1,*(a.degree_order() for a in rest)),0,-1):
            subtract(tuple(self._transverse(i,r,n) for i,r in enumerate(rest)))
            rest = tuple(self._jet_filter(a,n) for a in rest)
            if n == 1:
                subtract(tuple(self.firstpart(i,r,m) for i,r in enumerate(rest)))
                continue
            for i in range(2):
                ind = t.variable_index(i,n)
                for k in range(rest[i].max_exponent(ind),0,-1):
                    term = t.j(i,n)**k*self._transport_projection(
                        i,rest[i].coef(ind,k),m+k,n-1,1-k)
                    subtract((term,self.L()) if i == 0 else (self.L(),term))
                rest = tuple(a.coef(ind,0) if h == i else a for h,a in enumerate(rest))
        answer = c.denormalize(answer)
        residual = t.plus(c.bracket(c.input_flow(0),answer),t.mul(rhs,-1))
        return answer,residual

    def perturbation(self):
        field = 3*self.b**2
        ux,uxx,uxxx = (self.dx(field,k) for k in range(1,4))
        x2 = self.parameter(0)*uxxx+self.parameter(1)*uxx*ux
        if not self.constant:
            x2 += self.parameter(2)*ux**3
        return x2,x2

    def coefficients(self, vector):
        """Return exact state coefficient expressions grouped by physical jet monomial."""
        grouped = defaultdict(lambda: sp.S.Zero)
        for i, row in enumerate(vector):
            for key, value in row.items():
                monomial = tuple((j,e) for j,e in enumerate(key[:self.tag_start]) if j and e)
                coeff = sp.Rational(value)*BROOT**key[0]
                coeff *= sp.prod(x**e for x,e in zip(self.tags,key[self.tag_start:]) if e)
                grouped[i,monomial] += coeff
        return {key:sp.factor(c) for key,c in grouped.items() if c}

    def readable(self, expr):
        return expr.xreplace(dict(zip(self.tags,self.tag_values)))

    def ideal(self, coefficients):
        if not coefficients:
            return ()
        ordered = (self.tags[0],self.tags[self.max_derivative+1]) if self.constant else tuple(
            reversed(self.tags))
        field = sp.QQ.frac_field(BROOT)
        basis = sp.groebner(tuple(coefficients.values()),*ordered,domain=field)
        return tuple(self.readable(p.as_expr()) for p in basis.polys)


@dataclass(frozen=True)
class ScreeningResiduals:
    epsilon2_coefficients: dict
    epsilon2_ideal: tuple
    epsilon4_coefficients: dict
    epsilon4_ideal: tuple


def transform_residuals():
    """Compute epsilon² and epsilon⁴ transform residuals by graded recurrence.

    Epsilon⁴ is evaluated only after the differential conditions at
    epsilon²: a1,a2 constant and a3=0. Dict keys are (zero-based component,
    tuple of (Laurent jet variable index, exponent) pairs), using the
    same r_x,s_x,... basis as the A₂ archive; b is kept in coefficients.
    """
    arbitrary = _FormalCoefficients(constant=False)
    _,r2 = arbitrary.homological_projection(arbitrary.t.mul(arbitrary.perturbation(),-1),2)
    coefficients2 = arbitrary.coefficients(r2)
    constant = _FormalCoefficients(constant=True)
    x2 = constant.perturbation()
    g2,exact = constant.homological_projection(constant.t.mul(x2,-1),2)
    if any(exact):
        raise ArithmeticError("restricted epsilon² reconstruction did not verify")
    source4 = constant.t.mul(constant.c.bracket(g2,x2),Q(1,2))
    _,r4 = constant.homological_projection(source4,4)
    coefficients4 = constant.coefficients(r4)
    return ScreeningResiduals(
        {key:arbitrary.readable(c) for key,c in coefficients2.items()},
        arbitrary.ideal(coefficients2),
        {key:constant.readable(c) for key,c in coefficients4.items()},
        constant.ideal(coefficients4),
    )


def _leader(ctx):
    """One nontrivial Tsarev speed, not the x-flow or the input flow."""
    b,t = ctx.b,ctx.t
    v = t.mono(t.V_INDEX)
    leader = (v+9*b**3)*t.j(0,1),(v-9*b**3)*t.j(1,1)
    if any(ctx.c.bracket(ctx.c.input_flow(0),leader)):
        raise ArithmeticError("test leader does not commute at the dispersionless level")
    return leader


def _physical_poles(ctx, vector):
    t = ctx.t
    return tuple(ctx.L({k:v for k,v in row.items()
                        if k[1]<0 or k[2]<0 or k[t.log_index(0)] or k[t.log_index(1)]})
                 for row in vector)


def symmetry_residuals():
    """Project the selected symmetry directly, independently of the G projection.

    The graded recursion directly solves [X0,Y2]=[Y0,X2] then
    [X0,Y4]=[Y2,X2] (the higher input coefficients are zero). Any exact
    candidate is additionally checked for Laurent poles and logarithms.
    Both the unprojected bracket residual and the nonpolynomial pieces
    are retained as separate exact coefficients.
    """
    arb = _FormalCoefficients(constant=False)
    x2 = arb.perturbation()
    leader = _leader(arb)
    y2,r2 = arb.homological_projection(arb.c.bracket(leader,x2),3)
    coefs2 = arb.coefficients(r2)
    poles2 = arb.coefficients(_physical_poles(arb,y2))

    const = _FormalCoefficients(constant=True)
    x2 = const.perturbation()
    y2,exact2 = const.homological_projection(const.c.bracket(_leader(const),x2),3)
    if any(exact2) or any(_physical_poles(const,y2)):
        raise ArithmeticError("restricted second-order symmetry was not polynomial and exact")
    y4,r4 = const.homological_projection(const.c.bracket(y2,x2),5)
    coefs4 = const.coefficients(r4)
    poles4 = const.coefficients(_physical_poles(const,y4))

    return {
        order: {
            "residual": {k:ctx.readable(c) for k,c in coefs.items()},
            "poles": {k:ctx.readable(c) for k,c in poles.items()},
            "ideal": ctx.ideal({**{("residual",k):v for k,v in coefs.items()},
                                **{("poles",k):v for k,v in poles.items()}}),
        }
        for order,ctx,coefs,poles in ((2,arb,coefs2,poles2),(4,const,coefs4,poles4))
    }


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--symmetry",action="store_true",
                        help="also run a direct symmetry recursion independent of the transform")
    args = parser.parse_args()
    result = transform_residuals()
    for epsilon in (2,4):
        coefficients = getattr(result,f"epsilon{epsilon}_coefficients")
        ideal = getattr(result,f"epsilon{epsilon}_ideal")
        print(f"epsilon^{epsilon}: {len(coefficients)} physical residual coefficients; ideal {ideal}")
        for key, value in sorted(coefficients.items()):
            print(f"  {key}: {value}")
    if args.symmetry:
        for order,row in symmetry_residuals().items():
            print(f"direct symmetry epsilon^{order}: {row['ideal']}")
            for kind in ("residual","poles"):
                for key,value in sorted(row[kind].items()):
                    print(f"  {kind} {key}: {value}")


if __name__ == "__main__":
    main()
