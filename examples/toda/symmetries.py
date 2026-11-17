"""Reuse the saved Toda quasi-Miura generator and test three more hierarchy flows.

The two sides of each check share exact arithmetic, but not input formulas.
One side starts from a diagonal, dispersionless characteristic; the other
uses the published difference flows, independently Taylor expanded in their
original variables. Differences are compared as formal jet expressions.
"""
import argparse
import json
from fractions import Fraction as Q
from math import factorial
from pathlib import Path

from .compute import L, N, U, U_INDEX, V, V_INDEX, ZERO_KEY, bracket, d, dx, input_flow, j, log_index, mul, plus


CUTOFF = 7
FULL = CUTOFF + 1  # retain epsilon^7 inside difference quotients
DATA = Path(__file__).with_name("toda-epsilon6.json")
W = d**2/16
h = d/4


def gen_from_file(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data["metadata"]["cutoff"] != CUTOFF:
        raise ValueError("generator must be known through O(epsilon^7)")

    def monomials(terms):
        def key(exps):
            result = list(ZERO_KEY)
            for i,e in exps:
                if not 0 <= i < V_INDEX:
                    raise ValueError("unexpected variable in archived generator")
                result[i] = e
            return tuple(result)
        return L({key(exps): Q(a) for exps,a in terms})

    generator = {int(k): tuple(monomials(row) for row in rows)
                 for k,rows in data["generator"].items()}
    if generator.keys() != {2,4,6} or any(len(row)!=2 for row in generator.values()):
        raise ValueError("expected both components of G2, G4 and G6")
    return generator


def leader(flow):
    """Dispersionless speeds in the order (r,s), with CDZ time normalization."""
    if flow == "t21":
        speeds = h*(V+h), -h*(V-h)
    elif flow == "t22":
        c = V**2/2 + W
        speeds = V*W+h*c, V*W-h*c
    elif flow == "t11":
        speeds = V+h*U, V-h*U
    else:
        raise ValueError(f"unknown Toda flow {flow!r}")
    return tuple(w*j(i,1) for i,w in enumerate(speeds))


def transport(generator, characteristic):
    """The original exp(ad_G), now acting on a different leading flow."""
    b2,b4,b6 = (bracket(generator[k], characteristic) for k in (2,4,6))
    c22 = bracket(generator[2], b2)
    c24 = bracket(generator[2], b4)
    c42 = bracket(generator[4], b2)
    c222 = bracket(generator[2], c22)
    return {
        0: characteristic,
        2: b2,
        4: plus(b4, mul(c22, Q(1,2))),
        6: plus(plus(b6, mul(plus(c24,c42), Q(1,2))), mul(c222, Q(1,6))),
    }


def formal(value, precision=FULL):
    return (value,) + (L(),)*(precision-1)


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def times(c,a):
    return tuple(c*x for x in a)


def multiply(a,b):
    return tuple(sum((a[k]*b[n-k] for k in range(n+1)),L())
                 for n in range(min(len(a),len(b))))


def shift(a, halfsteps):
    """Taylor expand exp(halfsteps*epsilon*D/2) on a finite series."""
    t = Q(halfsteps,2)
    res = [L() for _ in a]
    for k, f in enumerate(a):
        for n in range(len(a)-k):
            if not f:
                break
            res[k+n] = res[k+n] + f*(t**n/factorial(n))
            if n+1 < len(a)-k:
                f = dx(f)
    return tuple(res)


def quotient_eps(a):
    if a[0]:
        raise ArithmeticError("reference difference did not vanish at epsilon^0")
    return a[1:]


def lax_product(a,b):
    """Multiply raw difference operators using (a Λ^i)(b Λ^j)=a(Λ^i b)Λ^(i+j)."""
    out = {}
    for i,left in a.items():
        for j,right in b.items():
            term = multiply(left,shift(right,2*i))
            out[i+j] = add(out[i+j],term) if i+j in out else term
    return out


def reference_lax(m):
    """CDZ (2.10) in raw u,v, evaluated in the centered jet algebra."""
    v = shift(formal(V),1)
    w = formal(W)
    op = {1:formal(L(1)), 0:v, -1:w}
    power = {0:formal(L(1))}
    for _ in range(m):
        power = lax_product(power,op)
    positive = {i:value for i,value in power.items() if i >= 0}
    first, second = lax_product(positive,op),lax_product(op,positive)
    du = times(L(Q(1,factorial(m)))/W,quotient_eps(sub(first[-1],second[-1])))
    dv = times(Q(1,factorial(m)),quotient_eps(sub(first[0],second[0])))
    return du,dv


def bernoulli(a, plus_sign):
    """B_+=z/(e^z−1); B_−=z/(1−e^(−z)) in CDZ (2.30)."""
    coef = (Q(1), Q(-1 if plus_sign else 1,2), Q(1,12), 0,
            Q(-1,720), 0, Q(1,30240), 0)
    result=[L() for _ in a]
    for k,result_at_k in enumerate(a):
        for n in range(len(a)-k):
            if not result_at_k:
                break
            if coef[n]:
                result[n+k]+=result_at_k*coef[n]
            if n+1 < len(a)-k:
                result_at_k=dx(result_at_k)
    return tuple(result)


def reference_logarithmic():
    """CDZ (2.32) in raw u,v, expanded without using the generator."""
    v = shift(formal(V),1)
    u = formal(U)
    w = formal(W)
    bu = bernoulli(u,False)
    bv = bernoulli(v,True)
    two = formal(L(2))
    # All shifts before epsilon division act to order 7.
    flux_u = multiply(v,sub(bu,two))
    du = quotient_eps(add(sub(flux_u,shift(flux_u,-2)),
                         sub(shift(bv,2),shift(bv,-2))))
    flux_v_plus = multiply(shift(w,2),sub(shift(bu,2),two))
    flux_v_minus = multiply(w,sub(shift(bu,-2),two))
    dv = add(multiply(v,tuple(dx(a) for a in v))[:CUTOFF],
             quotient_eps(sub(flux_v_plus,flux_v_minus)))
    return du,dv


def raw_from_transported(vector):
    """Push forward by U=2log(d/4), V=(r+s)/2, then undo the half shift of v."""
    padded=tuple(vector.get(k,(L(),L())) for k in range(CUTOFF))
    udot=tuple(2*(r-s)/d for r,s in padded)
    vdot=shift(tuple((r+s)/2 for r,s in padded),1)
    return udot,vdot


def check_commutativity(flow, canonical):
    """Independent check with the Toda time, after transporting coefficients."""
    for n in range(0,CUTOFF,2):
        residual = (L(),L())
        for k in range(0,n+1,2):
            residual = plus(residual,bracket(input_flow(k),canonical[n-k]))
        if any(residual):
            raise ArithmeticError(f"{flow}, order {n}: nonzero Toda commutator")


def verify(path=DATA, flows=("t21","t22","t11"), *, output=None):
    generator=gen_from_file(path)
    report={}
    results={}
    for flow in flows:
        canonical=transport(generator,leader(flow))
        actual=raw_from_transported(canonical)
        expected=reference_logarithmic() if flow=="t11" else reference_lax(
            2 if flow=="t21" else 3
        )
        terms=[]
        for n,(actual_u,actual_v) in enumerate(zip(*actual)):
            residual=(actual_u-expected[0][n],actual_v-expected[1][n])
            if any(residual):
                raise ArithmeticError(f"{flow}, order {n}: {[len(x) for x in residual]} nonzero terms")
            for a in (actual_u,actual_v):
                if any(any(exps[i]<0 for i in range(1,1+2*N))
                       or exps[log_index(0)] or exps[log_index(1)] for exps in a):
                    raise ArithmeticError(f"{flow}, order {n}: not polynomial in positive jets")
            terms.append((len(actual_u),len(actual_v)))
        check_commutativity(flow,canonical)
        report[flow]=terms
        results[flow]={
            "raw": {str(n):[a.serial(),b.serial()] for n,(a,b) in enumerate(zip(*actual))},
            "centered": {str(n):[a.serial(),b.serial()] for n,(a,b) in canonical.items()},
        }
    if output is not None:
        metadata = {
            "cutoff": CUTOFF,
            "raw": "(u_t,v_t) as Laurent monomials in the centered state and jets",
            "centered": "(r_t,s_t); odd coefficients zero",
            "state": {"0":"d=r-s", str(V_INDEX):"V=(r+s)/2", str(U_INDEX):"U=2 log(d/4)"},
            "jets": "indices 1..26 alternate r_x,s_x,r_xx,s_xx,...",
            "encoding": "each monomial is [list of [variable_index, exponent], rational coefficient]",
        }
        Path(output).write_text(json.dumps({"metadata":metadata,"flows":results},
                                           separators=(",",":"))+"\n",encoding="utf-8")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, help="save the exact raw and centered flow coefficients")
    parser.add_argument("--flow", action="append", choices=("t21","t22","t11"))
    args=parser.parse_args()
    for flow,rows in verify(args.generator,args.flow or ("t21","t22","t11"),
                            output=args.output).items():
        print(f"{flow}: all 14 raw residual coefficients zero; monomial counts {rows}",flush=True)


if __name__=="__main__":
    main()
