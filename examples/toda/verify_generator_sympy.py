"""Cross-check the recorded Toda generators in SymPy's physical jet algebra."""
import argparse
import json
from functools import lru_cache
from pathlib import Path
from time import monotonic

import sympy as sp
from quasitriviality import JetSpace
from quasitriviality.jets import _expr


class GapJetSpace(JetSpace):
    """A local exact simplifier using d=r−s instead of repeated expansion."""

    def __init__(self):
        super().__init__(("r", "s"))
        self.delta = sp.Dummy("delta")

    @lru_cache(maxsize=3000)
    def shifted(self, expr):
        r, s = self.states
        shifted = sp.expand(expr.xreplace({r: self.delta + s}))
        if any(p.exp.is_Integer is not True or
               (p.exp < 0 and p.base.is_Symbol is not True and not p.base.is_number)
               for p in shifted.atoms(sp.Pow)):
            shifted = sp.expand(sp.cancel(shifted))
        return shifted

    def exact(self, expr):
        r, s = self.states
        return self.shifted(_expr(expr)).xreplace({self.delta: r-s})

    def dx(self, expr, times=1):
        out = _expr(expr)
        for _ in range(times):
            out = self.exact(sum(self.jet(i, s+1)*self.partial(out,i,s)
                                 for i,s in self.active_jets(out)))
        return out

    def bracket(self, left, right):
        left, right = self._characteristic(left), self._characteristic(right)
        q_cache = {(i, 0): q for i, q in enumerate(left)}
        p_cache = {(i, 0): p for i, p in enumerate(right)}
        return tuple(self.exact(self._evolve(left, p, q_cache)-self._evolve(right, q, p_cache))
                     for q,p in zip(left,right))


def verify(source_file):
    data = json.loads(Path(source_file).read_text(encoding="utf-8"))["generator"]
    space = GapJetSpace()
    r,s = space.states
    h=(r-s)/4
    V=(r+s)/2
    indexed = ((r-s,) + tuple(space.jet(i,k) for k in range(1,14) for i in range(2))
               + (space.log(0),space.log(1)))
    g = {int(k):tuple(sp.Add(*(sp.Rational(c)*sp.prod(indexed[i]**pow for i,pow in exps)
                               for exps,c in item)) for item in pair)
         for k,pair in data.items()}
    def x(k):
        divisor=sp.factorial(k+1)*2**k
        a,b=space.dx(h*h,k+1),h*space.dx(V,k+1)
        return tuple(space.exact((a+sign*b)/divisor) for sign in (1,-1))
    bracket=space.bracket
    t=monotonic()
    gx0_2=bracket(g[2],x(0))
    if not all(space.exact(a-b)==0 for a,b in zip(gx0_2,x(2))):
        raise ArithmeticError("nonzero Lie residual at order 2")
    print("verified ε² in %.2fs" %(monotonic()-t),flush=True)
    gx0_4=bracket(g[4],x(0))
    b22=bracket(g[2],gx0_2)
    if not all(space.exact(a+b/2-c)==0 for a,b,c in zip(gx0_4,b22,x(4))):
        raise ArithmeticError("nonzero Lie residual at order 4")
    print("verified ε⁴ in %.2fs" %(monotonic()-t),flush=True)
    gx0_6=bracket(g[6],x(0))
    b24=bracket(g[2],gx0_4)
    b42=bracket(g[4],gx0_2)
    b222=bracket(g[2],b22)
    if not all(space.exact(a+(b+c)/2+d/6-e)==0
               for a,b,c,d,e in zip(gx0_6,b24,b42,b222,x(6))):
        raise ArithmeticError("nonzero Lie residual at order 6")
    print("verified ε⁶ in %.2fs" %(monotonic()-t),flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("source_file")
    verify(parser.parse_args().source_file)
