"""Exact Toda benchmark using the coupled hodograph recurrence.

The sparse algebra stores rational coefficients of monomials in r−s and
the jets. The transverse and own-jet integrations are the same as in
quasitriviality.coupled. Only the even, translation-invariant Toda chart
is supported by this specialized backend.
"""
import argparse
import json
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial


MAX_JET_ORDER = 13  # the O(epsilon**7) checks use orders strictly below this
N = MAX_JET_ORDER
VARIABLE_COUNT = 1 + 2*N + 2
ZERO_KEY = (0,)*VARIABLE_COUNT


def variable_index(i, n):
    if i not in (0, 1) or not 1 <= n <= MAX_JET_ORDER:
        raise IndexError(f"jet index ({i}, {n}) is outside the Toda benchmark domain")
    return 1 + 2*(n-1) + i


def log_index(i):
    if i not in (0, 1):
        raise IndexError(f"invalid logarithm component {i}")
    return 1 + 2*N + i


def power_at(key, i):
    return key[i]


def move_power(key, i, increment):
    return key[:i] + (key[i] + increment,) + key[i+1:]


def product_key(left, right):
    return tuple(a+b for a, b in zip(left, right))


def quotient_key(left, right):
    return tuple(a-b for a, b in zip(left, right))


class L(dict):
    def __init__(self, x=0):
        if isinstance(x, dict):
            super().__init__({k:Q(v) for k,v in x.items() if v})
        else:
            super().__init__({ZERO_KEY: Q(x)} if x else {})
    def __add__(a,b):
        if not isinstance(b,L): b=L(b)
        r=L(a)
        for k,v in b.items():
            c=r.get(k,0)+v
            if c:r[k]=c
            elif k in r:del r[k]
        return r
    __radd__=__add__
    def __neg__(a): return L({k:-v for k,v in a.items()})
    def __sub__(a,b): return a+(-b)
    def __rsub__(a,b): return (-a)+b
    def __mul__(a,b):
        if isinstance(b,(int,Q)):
            return L({k:v*b for k,v in a.items()}) if b else L()
        if not a or not b:return L()
        if len(a)==1:
            (k,v),=a.items()
            return L({product_key(k,j):v*w for j,w in b.items()})
        if len(b)==1:return b*a
        r={}
        for k,v in a.items():
            for j,w in b.items():
                i=product_key(k,j)
                r[i]=r.get(i,0)+v*w
        return L(r)
    __rmul__=__mul__
    def __truediv__(a,b):
        if not isinstance(b,L):return a* (1/Q(b))
        (k,v),=b.items()
        return L({quotient_key(j,k):w/v for j,w in a.items()})
    def __pow__(a,n):
        if len(a)==1:
            (k,v),=a.items()
            return L({tuple(i*n for i in k):v**n})
        if n<0:raise ValueError("nonmonomial denominator")
        r=L(1)
        while n:
            if n&1:r=r*a
            n>>=1
            if n:a=a*a
        return r
    def partial(a,i):
        return L({move_power(k,i,-1):v*n for k,v in a.items() if (n:=power_at(k,i))})
    def coef(a,i,p):
        return L({move_power(k,i,-p):v for k,v in a.items() if power_at(k,i)==p})
    def max_exponent(a,i):
        return max((power_at(k,i) for k in a),default=0)
    def depends(a,i):
        return any(power_at(k,i) for k in a)
    def degree_order(a):
        return max(
            ((i-1)//2+1 for i in range(1,1+2*N) if a.depends(i)),
            default=int(a.depends(log_index(0)) or a.depends(log_index(1))),
        )
    def substitute(a,func):
        if not a:return L()
        symbols=list(range(1+2*N+2))
        active=[i for i in symbols if a.depends(i)]
        pow_cache={}
        out=L()
        for k,c in a.items():
            m=L(c)
            for i in active:
                e=power_at(k,i)
                if e:
                    idx=i,e
                    if idx not in pow_cache:pow_cache[idx]=func(i)**e
                    m=m*pow_cache[idx]
            out+=m
        return out
    def active(a):return [i for i in range(1+2*N+2) if a.depends(i)]
    def serial(a):
        return [[[(i,e) for i in range(1+2*N+2) if (e:=power_at(k,i))], str(v)]
                for k,v in sorted(a.items()) if v]

def mono(i,n=1): return L({move_power(ZERO_KEY,i,n):1})
def j(i,n): return mono(variable_index(i,n))
def log(i): return mono(log_index(i))
d=mono(0)

def dx(a,times=1):
    for _ in range(times):
        acc=a.partial(0)*(j(0,1)-j(1,1))
        for i in range(2):
            for n in range(1,a.degree_order()+1):
                if b:=a.partial(variable_index(i,n)):acc+=b*j(i,n+1)
            if b:=a.partial(log_index(i)):acc+=b*j(i,2)/j(i,1)
        a=acc
    return a

def evol(vector,target):
    maxorder=max((a.degree_order() for a in target),default=0)
    applicable=[]
    for a in vector:
        derivatives=[a]
        for n in range(maxorder):derivatives.append(dx(derivatives[-1]))
        applicable.append(derivatives)
    out=[]
    for a in target:
        acc=a.partial(0)*(vector[0]-vector[1])
        for i in range(2):
            for n in range(1,a.degree_order()+1):
                if b:=a.partial(variable_index(i,n)):acc+=b*applicable[i][n]
            if b:=a.partial(log_index(i)):acc+=b*applicable[i][1]/j(i,1)
        out.append(acc)
    return tuple(out)

def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(a,b):return tuple(x*b for x in a)
def bracket(a,b):return plus(evol(a,b),mul(evol(b,a),-1))

def input_flow(order):
    U=dx(d*d/16,order+1)
    V=(j(0,order+1)+j(1,order+1))*d/8
    c=Q(1,(2**order)*factorial(order+1))
    return (U+V)*c,(U-V)*c

def vspeed(i):return d * (Q(1,4) if i==0 else Q(-1,4))

def dstate(expr,i):return expr.partial(0)*(1 if i==0 else -1)

def gamma(i,jj):return dstate(vspeed(i),jj)/(vspeed(jj)-vspeed(i))
@lru_cache(maxsize=None)
def value_coeff(i,k,n):
    if n==0:return gamma(i,k)
    prev=value_coeff(i,k,n-1)
    return dstate(prev,i)-gamma(k,i)*prev
@lru_cache(maxsize=None)
def jet_coeff(i,k,n,t):
    if n==0:return L()
    if t==0:return dstate(jet_coeff(i,k,n-1,0),i)-value_coeff(i,k,n-1)
    return dstate(jet_coeff(i,k,n-1,t),i)+jet_coeff(i,k,n-1,t-1)
@lru_cache(maxsize=None)
def physical_to_hodo_symbol(i,n):
    # physical p_{i,n} as an expression of inverse jets z.
    if n==1:return j(i,1)**-1
    expr = hodo_to_physical_symbol(i,n) + j(i,n)/j(i,1)**(n+1)
    return (expr.substitute(to_hodo_var)-j(i,n))/j(i,1)**(n+1)

@lru_cache(maxsize=None)
def hodo_to_physical_symbol(i,n):
    # inverse jet z_{i,n} as an expression of physical p.
    if n==1:return j(i,1)**-1
    prev=hodo_to_physical_symbol(i,n-1)
    k=1-i
    c = sum((jet_coeff(i,k,n-1,t)*hodo_to_physical_symbol(i,t+1)
             for t in range(n-1)),L())
    return (dx(prev)-j(k,1)*c)/j(i,1)

def to_hodo_var(k):
    if k==0:return d
    if k>=1+2*N:return -mono(k)
    return physical_to_hodo_symbol((k-1)%2,(k-1)//2+1)
def to_physical_var(k):
    if k==0:return d
    if k>=1+2*N:return -mono(k)
    return hodo_to_physical_symbol((k-1)%2,(k-1)//2+1)

def normalize(a):return tuple((x/j(i,1)).substitute(to_hodo_var) for i,x in enumerate(a))
def denormalize(a):return tuple(x.substitute(to_physical_var)*j(i,1) for i,x in enumerate(a))

@lru_cache(maxsize=None)
def transport_velocity(row,i,n):
    ans=(vspeed(i)-vspeed(row))*j(i,n+1)/j(i,1)
    k=1-i
    if k!=row:
        ans+=(vspeed(k)-vspeed(row))*sum(
            (jet_coeff(i,k,n,t)*j(i,t+1) for t in range(n)),L()
        )/j(k,1)
    if n==1:ans-=dstate(vspeed(i),i)
    return ans
def transport(row,a,factor):
    k=1-row
    ans=(vspeed(k)-vspeed(row))*dstate(a,k)/j(k,1)+factor*dstate(vspeed(row),k)*a/j(k,1)
    for i in range(2):
        for n in range(1,a.degree_order()+1):
            if v:=a.partial(variable_index(i,n)):ans+=transport_velocity(row,i,n)*v
        if v:=a.partial(log_index(i)):ans+=transport_velocity(row,i,1)*v/j(i,1)
    return ans
def hodo_action(a):
    return tuple(transport(i,a[i],1)-dstate(vspeed(i),1-i)*a[1-i]/j(1-i,1)
                 for i in range(2))

def primitive(a,i,n):
    ind=variable_index(i,n)
    out={}
    for e,c in a.items():
        m=power_at(e,ind)+1
        if n!=1:
            p=move_power(e,ind,1)
            out[p]=out.get(p,0)+c/m
            continue
        k=power_at(e,log_index(i))
        if m==0:
            p=move_power(move_power(e,ind,1),log_index(i),1)
            out[p]=out.get(p,0)+c/(k+1)
        else:
            for h in range(k+1):
                p=move_power(move_power(e,ind,1),log_index(i),-h)
                out[p]=out.get(p,0)+c*Q((-1)**h*factorial(k),factorial(k-h)*m**(h+1))
    return L(out)

def firstpart(row,rhs,weight):
    ind=variable_index(row,1)
    ell=log_index(row)
    slope=dstate(vspeed(row),row)
    # Keep the own z exponent and log; no other auxiliary variables.
    coeff={}
    for e,v in rhs.items():
        if power_at(e,ind)!=-weight-1:continue
        k=power_at(e,ell)
        base=move_power(move_power(e,ind,weight+1),ell,-k)
        if any(base[1:]):continue
        if weight==0:
            p=move_power(base,ell,k+1)
            coeff[p]=coeff.get(p,0)-v/(k+1)
        else:
            for h in range(k+1):
                p=move_power(move_power(base,ind,-weight),ell,k-h)
                coeff[p]=coeff.get(p,0)+v*Q(factorial(k),factorial(k-h)*weight**(h+1))
    return L(coeff)/slope

def top_transverse(row,rhs,n):
    k=1-row
    f=rhs.partial(variable_index(k,n+1))*j(k,1)/(vspeed(k)-vspeed(row))
    require(f.degree_order()<=n, f"nonlinear highest transverse jet in row {row}, order {n}")
    return primitive(f,k,n)

def solve_transport(row,rhs,weight,order,factor):
    if not rhs:return L()
    require(order>0, "no positive-order primitive")
    top=top_transverse(row,rhs,order)
    rest=rhs-transport(row,top,factor)
    require(rest.degree_order()<=order, f"uncancelled highest jet in row {row}, order {order}")
    if order==1:
        result=top+firstpart(row,rest,weight)
    else:
        ind=variable_index(row,order)
        for k in range(rest.max_exponent(ind),0,-1):
            c=rest.coef(ind,k)
            lower=solve_transport(row,c,weight+k,order-1,factor-k)
            term=lower*j(row,order)**k
            rest-=transport(row,term,factor)
            require(not rest.coef(ind,k), f"uncancelled own-jet power {k}")
            top+=term
        require(not rest.depends(ind), "uncancelled highest own jet")
        result=top+solve_transport(row,rest,weight,order-1,factor)
    require(transport(row,result,factor)==rhs, f"scalar transport residual in row {row}, order {order}")
    return result

def solve(rhs,weight):
    m=weight-1
    rest=normalize(rhs)
    ans=(L(),L())
    def subtract(parts):
        nonlocal rest,ans
        if not any(parts):return
        rest=plus(rest,mul(hodo_action(parts),-1))
        ans=plus(ans,parts)
    for n in range(max(1,*(a.degree_order() for a in rest)),0,-1):
        subtract(tuple(top_transverse(i,r,n) for i,r in enumerate(rest)))
        require(all(a.degree_order()<=n for a in rest), f"uncancelled top transverse jet at order {n}")
        if n==1:
            subtract(tuple(firstpart(i,r,m) for i,r in enumerate(rest)))
            continue
        for i in range(2):
            ind=variable_index(i,n)
            while rest[i].depends(ind):
                k=rest[i].max_exponent(ind)
                require(k>0, "unsupported negative higher-jet dependence")
                term=j(i,n)**k*solve_transport(i,rest[i].coef(ind,k),m+k,n-1,1-k)
                subtract((term,L()) if i==0 else (L(),term))
                require(not rest[i].coef(ind,k), f"uncancelled row {i} power {k}")
    require(not any(rest), f"nonzero final normalized residual at weight {weight}")
    ansphys=denormalize(ans)
    require(bracket(input_flow(0),ansphys)==rhs, f"nonzero physical bracket residual at weight {weight}")
    return ansphys

def inverses(g2,g4,g6):
    # At order 2, 4, 6 in exp(±E_G) acting on states.
    e22=evol(g2,g2)
    e24=evol(g2,g4); e42=evol(g4,g2)
    e222=evol(g2,e22)
    return {
      "reduce": (g2, plus(g4,mul(e22,Q(1,2))),
                  plus(plus(g6,mul(plus(e24,e42),Q(1,2))),mul(e222,Q(1,6)))),
      "dress": (mul(g2,-1), plus(mul(g4,-1),mul(e22,Q(1,2))),
                 plus(plus(mul(g6,-1),mul(plus(e24,e42),Q(1,2))),mul(e222,-Q(1,6))))
    }


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def x_derivative_leading_flat(component, order):
    """D**order of U0 or V0, with the undifferentiated constants excluded."""
    if order < 1:
        raise ValueError("the translation-invariant algebra stores derivatives of U0 or V0")
    if component == 0:
        return dx(2*(j(0,1)-j(1,1))/d, order-1)
    return (j(0,order)+j(1,order))/2


def centered_flat_map(dressing):
    """Convert dispersionless -> dressed Riemann coordinates to centered U,V."""
    delta = {k: (r-s)/d for k,(r,s) in zip((2,4,6),dressing)}
    U = {
        2: 2*delta[2],
        4: 2*delta[4] - delta[2]**2,
        6: 2*delta[6] - 2*delta[2]*delta[4] + Q(2,3)*delta[2]**3,
    }
    V = {k: (r+s)/2 for k,(r,s) in zip((2,4,6),dressing)}
    exponential = {
        0: d*d/16,
        2: d*d/8*delta[2],
        4: d*d/16*(2*delta[4] + delta[2]**2),
        6: d*d/8*(delta[6] + delta[2]*delta[4]),
    }
    return U, V, exponential


def raw_flat_map(U, V):
    """Undo the half-shift V=e^(-eps D/2)v; U=u is unchanged."""
    raw = {}
    for n in range(1,7):
        shifted = Q(1,2**n*factorial(n))*x_derivative_leading_flat(1,n)
        for k in range(2,n+1,2):
            m=n-k
            shifted += Q(1,2**m*factorial(m))*dx(V[k],m)
        raw[n]=(U.get(n,L()),shifted)
    return raw


def verify_direct(U, V, exponential, raw, leading):
    """Check both PDEs directly in centered and original lattice variables."""
    # E_X0 V0 = (X0^r+X0^s)/2, E_X0 U0 = 2(X0^r-X0^s)/d.
    lead_t=(2*(leading[0]-leading[1])/d,(leading[0]+leading[1])/2)
    for n in range(0,7,2):
        left = lead_t if n==0 else evol(leading,(U[n],V[n]))
        right_u = Q(1,2**n*factorial(n+1))*x_derivative_leading_flat(1,n+1)
        right_v = L()
        for k in range(0,n+1,2):
            c = Q(1,2**k*factorial(k+1))
            right_v += c*dx(exponential[n-k],k+1)
            if n-k:
                right_u += c*dx(V[n-k],k+1)
        require(left==(right_u,right_v),f"centered lattice chain-rule residual at epsilon^{n}")
    for n in range(0,7):
        left = lead_t if n==0 else evol(leading,raw[n])
        right_u = Q((-1)**n,factorial(n+1))*x_derivative_leading_flat(1,n+1)
        right_v = L()
        for k in range(n+1):
            if n-k:
                right_u += Q((-1)**k,factorial(k+1))*dx(raw[n-k][1],k+1)
            if n-k in exponential:
                right_v += Q(1,factorial(k+1))*dx(exponential[n-k],k+1)
        require(left==(right_u,right_v),f"original lattice chain-rule residual at epsilon^{n}")


def total_primitive(source):
    """Invert D_x on the positive-weight total derivatives used here.

    Recover each highest-jet coefficient from the next jet in the source.
    A final nonzero remainder is an error rather than an integration guess.
    """
    rest=source
    result=L()
    for order in range(rest.degree_order()-1,0,-1):
        forced=tuple(rest.partial(variable_index(i,order+1)) for i in range(2))
        require(all(a.degree_order()<=order for a in forced),
                f"nonlinear top derivative in total primitive at order {order}")
        part=L()
        for i, item in enumerate(forced):
            derivative=part.partial(variable_index(i,order))
            if order==1:
                derivative+=part.partial(log_index(i))/j(i,1)
            part+=primitive(item-derivative,i,order)
        rest-=dx(part)
        require(rest.degree_order()<=order,
                f"uncancelled top derivative in total primitive at order {order}")
        result+=part
    require(not rest and dx(result)==source, "total primitive has a nonzero remainder")
    return result


def recover_free_energies(U,V,leading):
    """Extract F2,F3 from the computed map, with no imported genus data.

    Apply the inverse centered Miura shifts B=(epsilon D)/(2 sinh(epsilon D/2))
    and B**2 to the flat coordinates. In tau coordinates the corrections
    are D E_0 F_g and D**2 F_g. Check both, and the tau PDE's t,t derivative.
    """
    B={0:Q(1), 2:Q(-1,24), 4:Q(7,5760), 6:Q(-31,967680)}
    B2={0:Q(1), 2:Q(-1,12), 4:Q(1,240), 6:Q(-1,6048)}
    q2=d*d/16
    results={}
    for n in (4,6):
        uhat=B2[n]*x_derivative_leading_flat(0,n)
        vhat=B[n]*x_derivative_leading_flat(1,n)
        for k in range(2,n+1,2):
            uhat+=B2[n-k]*dx(U[k],n-k)
            vhat+=B[n-k]*dx(V[k],n-k)
        F=total_primitive(total_primitive(uhat))
        derivative,=evol(leading,(F,))
        require(dx(derivative)==vhat,f"tau-coordinate V check at epsilon^{n}")
        tt,=evol(leading,(derivative,))
        exp_term=(U[4]+U[2]**2/2) if n==4 else (
            U[6]+U[2]*U[4]+U[2]**3/6
        )
        require(tt==q2*exp_term,f"tau equation check at epsilon^{n}")
        results[n]=(F,)
    return results


def verify_reducing(reducing, flow):
    """Check E_X of the reducing map equals the dispersionless flow of its outputs."""
    coeffs = dict(zip((2,4,6), reducing))
    for n in range(0,7,2):
        left = flow[n]
        for k in range(0,n,2):
            left = plus(left,evol(flow[k],coeffs[n-k]))
        right = [L(),L()]
        for k in range(0,n+1,2):
            gap = d if k == 0 else coeffs[k][0]-coeffs[k][1]
            terms = (j(0,1),j(1,1)) if k == n else tuple(dx(f) for f in coeffs[n-k])
            for i in range(2):
                right[i] += gap*terms[i]*Q(1 if i == 0 else -1,4)
        require(left==tuple(right),f"reducing chain-rule residual at epsilon^{n}")


def compute():
    """Generate both quasi-Miura directions from the centered equation alone."""
    x0,x2,x4,x6=map(input_flow,(0,2,4,6))
    g2=solve(mul(x2,-1),2)
    r4=plus(x4,mul(bracket(g2,x2),-Q(1,2)))
    g4=solve(mul(r4,-1),4)
    r6=plus(x6,mul(plus(bracket(g2,r4),bracket(g4,x2)),-Q(1,2)))
    r6=plus(r6,mul(bracket(g2,bracket(g2,x2)),-Q(1,6)))
    g6=solve(mul(r6,-1),6)
    maps=inverses(g2,g4,g6)
    maps["generator"]=(g2,g4,g6)
    verify_reducing(maps["reduce"],{0:x0,2:x2,4:x4,6:x6})
    U,V,exp_U=centered_flat_map(maps["dress"])
    raw=raw_flat_map(U,V)
    verify_direct(U,V,exp_U,raw,x0)
    result={name:{2*k+2: value for k,value in enumerate(coefficients)}
            for name,coefficients in maps.items()}
    result["centered_flat_dress"]={k:(U[k],V[k]) for k in (2,4,6)}
    result["original_flat_dress"]=raw
    result["free_energy"]=recover_free_energies(U,V,x0)
    return result


def serialize(result):
    """Encode each term as [list of [variable index, integer exponent], coefficient]."""
    result={key:{str(order):[a.serial() for a in vector] for order,vector in rows.items()}
            for key,rows in result.items()}
    result["metadata"]={
        "cutoff":7,
        "basis":{"0":"d = r - s", "1 ... 26":"r_x, s_x, r_xx, s_xx, ... in alternating order",
                 str(log_index(0)):"log|r_x|", str(log_index(1)):"log|s_x|"},
        "leading_r_s":"(r,s), d=4exp(U0/2), V0=(r+s)/2",
        "leading_flat":"(U0,V0); raw continuum V0 agrees with centered V0 at epsilon^0",
        "coordinates":"centered V=exp(-epsilon D_x/2) v, U=u",
        "free_energy":"order 4 holds F2 and order 6 holds F3, reconstructed from the tau-coordinate map",
        "verification":"homological brackets, centered and original chain rules, and both potential identities for F2/F3",
    }
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o","--output",type=Path,help="write the full exact coefficient JSON")
    args=parser.parse_args(argv)
    result=compute()
    if args.output:
        encoded = serialize(result)
        blocks = [json.dumps(key) + ": {\n" + ",\n".join(
            "  " + json.dumps(order) + ": " + json.dumps(row, separators=(",", ":"))
            for order, row in rows.items()
        ) + "\n}" for key, rows in encoded.items()]
        args.output.write_text("{\n" + ",\n".join(blocks) + "\n}\n", encoding="utf-8")
    print("Verified centered and raw Toda through O(epsilon^7).")
    for name, rows in result.items():
        print(name, {k:tuple(map(len,value)) for k,value in rows.items()})


if __name__=="__main__":
    main()
