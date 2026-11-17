"""Exact equation-only A2 (Boussinesq) quasi-Miura computation through O(epsilon**7).

Use the Dubrovin–Zhang flat normalization
    U_t = V_x,
    V_t = 3 U U_x + epsilon**2 U_xxx / 4.
The sparse algebra stores rational Laurent monomials in b=sqrt(U/3)
and physical or generalized hodograph jets. The solver integrates each
coefficient of the flow before the free energies are extracted.
"""

import argparse
import json
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial

if __package__ in (None, ""):
    # Also support python examples/a2/compute.py from a repository checkout.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from examples.toda.compute import (
    L, N, j, log, log_index, mono, move_power, mul, plus, power_at, primitive,
    require, variable_index,
)

b = mono(0)


def dx(a,times=1):
    for _ in range(times):
        gap_x = j(0,1)-j(1,1)
        acc = a.partial(0)*gap_x/(36*b**2)
        for i in range(2):
            for n in range(1,a.degree_order()+1):
                if coefficient:=a.partial(variable_index(i,n)):acc+=coefficient*j(i,n+1)
            if coefficient:=a.partial(log_index(i)):acc+=coefficient*j(i,2)/j(i,1)
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
        gap_t = vector[0]-vector[1]
        acc = a.partial(0)*gap_t/(36*b**2)
        for i in range(2):
            for n in range(1,a.degree_order()+1):
                if coefficient:=a.partial(variable_index(i,n)):acc+=coefficient*applicable[i][n]
            if coefficient:=a.partial(log_index(i)):acc+=coefficient*applicable[i][1]/j(i,1)
        out.append(acc)
    return tuple(out)

def bracket(a,b):return plus(evol(a,b),mul(evol(b,a),-1))

def input_flow(order):
    """Riemann-coordinate epsilon coefficients, obtained from Boussinesq."""
    if order == 0:
        return tuple(vspeed(i)*j(i,1) for i in range(2))
    if order == 2:
        dispersion = dx(3*b**2,3)/4
        return dispersion,dispersion
    return L(),L()


def vspeed(i):
    return b * (3 if i==0 else -3)


def dstate(expr,i):
    # r - s = 12*b**3. Expressions in this benchmark are invariant
    # under translations of both Riemann coordinates.
    return expr.partial(0)*(1 if i==0 else -1)/(36*b**2)


def self_speed_derivative(i,n):
    derivative=vspeed(i)
    for _ in range(n):
        derivative=dstate(derivative,i)
    return derivative


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
    if k == 0:
        return b
    if k in (log_index(0), log_index(1)):
        return -mono(k)
    if not 1 <= k <= 2*N:
        raise ValueError("unexpected variable in the translation-invariant jet chart")
    return physical_to_hodo_symbol((k-1)%2,(k-1)//2+1)


def to_physical_var(k):
    if k == 0:
        return b
    if k in (log_index(0), log_index(1)):
        return -mono(k)
    if not 1 <= k <= 2*N:
        raise ValueError("unexpected variable in the translation-invariant jet chart")
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
    return ans-self_speed_derivative(i,n)
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




def dressed_flat(dressing):
    """Convert the inverse Riemann transformation into flat (U,V) fields."""
    delta = {n: (r-s)/(12*b**3) for n,(r,s) in zip((2,4,6),dressing)}
    factor = 3*b**2
    U = {
        0: factor,
        2: factor*Q(2,3)*delta[2],
        4: factor*(Q(2,3)*delta[4]-Q(1,9)*delta[2]**2),
        6: factor*(Q(2,3)*delta[6]-Q(2,9)*delta[2]*delta[4]+Q(4,81)*delta[2]**3),
    }
    V = {n: (r+s)/2 for n,(r,s) in zip((2,4,6),dressing)}
    return U,V


def verify_direct(U,V,leading):
    """Check the flat PDEs by directly differentiating the dressed fields."""
    for n in (0,2,4,6):
        left_u, = evol(leading,(U[n],))
        right_u = (j(0,1)+j(1,1))/2 if n==0 else dx(V[n])
        left_v = (leading[0]+leading[1])/2 if n==0 else evol(leading,(V[n],))[0]
        potential = sum((U[k]*U[n-k] for k in range(0,n+1,2)),L())
        right_v = Q(3,2)*dx(potential)
        if n:
            right_v += dx(U[n-2],3)/4
        require((left_u,left_v)==(right_u,right_v),f"flat chain-rule residual at epsilon^{n}")


def verify_reducing(reducing):
    """Directly check that the map sends the dispersive flow to X0."""
    corrections = dict(zip((2,4,6),reducing))
    delta = {n: (r-s)/(12*b**3) for n,(r,s) in corrections.items()}
    speed_root = {
        0: b,
        2: b*delta[2]/3,
        4: b*(delta[4]/3-delta[2]**2/9),
        6: b*(delta[6]/3-Q(2,9)*delta[2]*delta[4]+Q(5,81)*delta[2]**3),
    }
    for n in (0,2,4,6):
        left = input_flow(n)
        for k in range(0,n,2):
            left = plus(left,evol(input_flow(k),corrections[n-k]))
        right = []
        for i in range(2):
            component = speed_root[n]*j(i,1)
            for k in range(0,n,2):
                component += speed_root[k]*dx(corrections[n-k][i])
            right.append(component*(3 if i==0 else -3))
        require(left==tuple(right),f"reducing chain-rule residual at epsilon^{n}")


def total_primitive(source):
    """Integrate a verified total derivative over its finite jet support."""
    rest=source
    result=L()
    for order in range(rest.degree_order()-1,0,-1):
        forced=tuple(rest.partial(variable_index(i,order+1)) for i in range(2))
        require(all(a.degree_order()<=order for a in forced),
                f"nonlinear top derivative in total primitive at order {order}")
        part=L()
        for i,item in enumerate(forced):
            derivative=part.partial(variable_index(i,order))
            if order==1:
                derivative+=part.partial(log_index(i))/j(i,1)
            part+=primitive(item-derivative,i,order)
        rest-=dx(part)
        require(rest.degree_order()<=order,
                f"uncancelled top derivative in total primitive at order {order}")
        result+=part
    require(not rest and dx(result)==source,"total primitive has a nonzero remainder")
    return result


def recover_free_energies(U,V,leading):
    """Recover F1,F2,F3 from D_x^2 F_g=U_(2g) and verify both tau identities."""
    result={}
    for n in (2,4,6):
        F=total_primitive(total_primitive(U[n]))
        derivative,=evol(leading,(F,))
        require(dx(derivative)==V[n],f"tau-coordinate V check at epsilon^{n}")
        tt,=evol(leading,(derivative,))
        potential=sum((U[k]*U[n-k] for k in range(0,n+1,2)),L())
        require(tt==Q(3,2)*potential+dx(U[n-2],2)/4,f"tau PDE check at epsilon^{n}")
        result[n]=(F,)
    return result


def compute():
    x0,x2,x4,x6=map(input_flow,(0,2,4,6))
    g2=solve(mul(x2,-1),2)
    r4=plus(x4,mul(bracket(g2,x2),-Q(1,2)))
    g4=solve(mul(r4,-1),4)
    r6=plus(x6,mul(plus(bracket(g2,r4),bracket(g4,x2)),-Q(1,2)))
    r6=plus(r6,mul(bracket(g2,bracket(g2,x2)),-Q(1,6)))
    g6=solve(mul(r6,-1),6)
    maps=inverses(g2,g4,g6)
    maps["generator"]=(g2,g4,g6)
    U,V=dressed_flat(maps["dress"])
    verify_direct(U,V,x0)
    verify_reducing(maps["reduce"])
    result={name:{2*k+2:value for k,value in enumerate(coefficients)}
            for name,coefficients in maps.items()}
    result["flat_dress"]={k:(U[k],V[k]) for k in (2,4,6)}
    result["free_energy"]=recover_free_energies(U,V,x0)
    return result


def serialize(result):
    """Encode each term as [list of [variable index, integer exponent], coefficient]."""
    result={key:{str(order):[a.serial() for a in vector] for order,vector in rows.items()}
            for key,rows in result.items()}
    result["metadata"]={
        "cutoff":7,
        "basis":{"0":"b = sqrt(U0 / 3)",f"1 ... {2*N}":"r_x, s_x, r_xx, s_xx, ... alternating",
                 str(log_index(0)):"log|r_x|",str(log_index(1)):"log|s_x|"},
        "leading_r_s":"(r,s), r-s=12*b**3, V0=(r+s)/2",
        "leading_flat":"(U0,V0), U0=3*b**2",
        "coordinates":"U_t=V_x; V_t=3*U*U_x+epsilon**2*D_x**3(U)/4",
        "free_energy":"entries at epsilon order 2g contain F_g extracted from D_x**2(F_g)",
        "verification":"homological identities, direct dressing and reducing PDE chain rules, and tau potentials",
    }
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o","--output",type=Path,help="write the full exact coefficient JSON")
    args=parser.parse_args(argv)
    result=compute()
    if args.output:
        args.output.write_text(json.dumps(serialize(result))+"\n",encoding="utf-8")
    print("Verified A2 Boussinesq through O(epsilon^7).")
    for name,rows in result.items():
        print(name,{k:tuple(map(len,value)) for k,value in rows.items()})


if __name__=="__main__":
    main()
