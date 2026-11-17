"""Check the A2 Dubrovin–Zhang jet-space loop equation over exact Laurent polynomials.

The period gradients of the deformed intersection metric satisfy its
Levi-Civita connection. The checker contracts them using that metric, so
no choices of algebraic cubic roots or square-root branches enter.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
import json
from math import comb
from pathlib import Path

from ..toda.compute import L, V_INDEX, j, log_index, mono, variable_index
from ..toda.loop_equation import add_product, read_sparse
from .compute import b, dx as physical_dx


# The Laurent state symbols denote b = sqrt(U/3) and a = V - lambda.
# Here r,s = V +/- 6 b**3; lambda is constant under D_x.
a = mono(V_INDEX)
w = 9*b**2
delta = a**2 - 36*b**6
vx = (j(0,1)+j(1,1))/2
ux = (j(0,1)-j(1,1))/(6*b)

# Covariant intersection metric B = B_num / delta. The rows of H_V and
# H_U are the coefficients in partial_i (partial_alpha p) = H_i[alpha,:]
# dot (partial_V p, partial_U p) / delta, for any flat period p.
B_num = ((-2*b**2, a), (a, -18*b**4))
H_V = ((-2*a/3, 2*b**2/3), (12*b**4, -a/3))
H_U = ((12*b**4, -a/3), (-6*a*b**2, 6*b**4))
H_x = tuple(tuple(v*vx + u*ux for v,u in zip(row_v,row_u))
            for row_v,row_u in zip(H_V,H_U))
identity = ((L(1), L()), (L(), L(1)))


def dx(expr, times=1):
    for _ in range(times):
        expr = physical_dx(expr) + vx*expr.partial(V_INDEX)
    return expr


delta_x = dx(delta)


@lru_cache(maxsize=None)
def b_jet(order):
    return dx(b_jet(order-1)) if order else b


def flat_partial(expr, alpha, k):
    """Partial derivative with respect to a flat jet, holding all others fixed."""
    if alpha == 0:
        if not k:
            return expr.partial(V_INDEX)
        result = expr.partial(variable_index(0,k)) + expr.partial(variable_index(1,k))
        if k == 1:
            result += expr.partial(log_index(0))/j(0,1) + expr.partial(log_index(1))/j(1,1)
        return result
    result = expr.partial(0)/(6*b) if not k else L()
    for n in range(max(1,k),expr.degree_order()+1):
        difference = expr.partial(variable_index(0,n)) - expr.partial(variable_index(1,n))
        if n == 1:
            difference += (expr.partial(log_index(0))/j(0,1)
                           - expr.partial(log_index(1))/j(1,1))
        if difference:
            result += 3*comb(n,k)*b_jet(n-k)*difference
    return result


def derivatives(numerator, power, order):
    """Numerators of D_x^i(numerator / delta**power) over delta**(power+i)."""
    rows = [numerator]
    for k in range(order):
        rows.append(dx(rows[-1])*delta - (power+k)*rows[-1]*delta_x)
    return rows


def matmul(left, right):
    return tuple(tuple(left[i][0]*right[0][j] + left[i][1]*right[1][j]
                       for j in range(2)) for i in range(2))


def transpose(matrix):
    return tuple(zip(*matrix))


def period_numerators(order):
    """A_k / delta**k acts on the period covector to give its k-th x-jet."""
    rows = [identity]
    for k in range(order):
        prev = rows[-1]
        shifted = matmul(prev,H_x)
        rows.append(tuple(tuple(dx(prev[i][j])*delta - k*prev[i][j]*delta_x
                                + shifted[i][j] for j in range(2)) for i in range(2)))
    return rows


def connection_numerators():
    # Periods can be chosen as roots of z^3/9 - U z + V - lambda, hence
    # partial_lambda q = -partial_V q. Both minus signs cancel in the pairing.
    gram = matmul(matmul(H_V,B_num),transpose(H_V))
    # The metric eta raises a covector by swapping its entries.
    # The Frobenius product satisfies e_1=unit and e_2^2=w e_1.
    c1, c2 = gram[1][1] + w*gram[0][0], 2*gram[1][0]
    return (c1*vx+w*c2*ux)/2, (c2*vx+c1*ux)/2


def read_free_energies(archive):
    data = json.loads(Path(archive).read_text(encoding="utf-8"))
    return tuple(read_sparse(data["free_energy"][str(2*g)][0]) for g in (1,2,3))


def coefficients(free_energies):
    """Return (power, cleared residual) for epsilon^0, epsilon^2, epsilon^4."""
    if len(free_energies) != 3:
        raise ValueError("expected F1, F2, F3")
    orders = [f.degree_order() for f in free_energies]
    order = max(max(orders), max(orders[:2])+1)
    indices = [[(alpha,k) for alpha in range(2) for k in range(n+1)] for n in orders]
    grads = [{(alpha,k): flat_partial(f,alpha,k) for alpha,k in idx}
             for f,idx in zip(free_energies,indices)]
    # In flat coordinates (E-lambda)^(-1)=(a,-2*b**2)/delta.
    resolvent = [derivatives(n,1,order) for n in (a,-2*b**2)]
    connection = [derivatives(n,3,max(orders[:2])+1) for n in connection_numerators()]
    periods = period_numerators(order)
    weighted = [matmul(period,B_num) for period in periods]

    @lru_cache(maxsize=None)
    def pair(alpha,k,beta,l):
        if (alpha,k) > (beta,l):
            return pair(beta,l,alpha,k)
        return (weighted[k][alpha][0]*periods[l][beta][0]
                + weighted[k][alpha][1]*periods[l][beta][1])

    residuals = []
    for g in range(3):
        terms = defaultdict(dict)
        for (alpha,k),grad in grads[g].items():
            add_product(terms,k+1,grad,resolvent[alpha][k])
            if k and grad:
                mixed = sum((comb(k,m)*pair(0,m-1,1-alpha,k-m+1)
                             for m in range(1,k+1)),L())
                add_product(terms,k+1,grad,mixed)
        if not g:
            # 1/2 times the integrated period pairing, after cancellation
            # of kappa_0/lambda^2, is -(a**2 + 45*b**6)/(9*delta**2).
            add_product(terms,2,(a**2+45*b**6)/9)
        else:
            for (alpha,k),grad in grads[g-1].items():
                add_product(terms,4+k,-grad,connection[alpha][k+1])
                if grad:
                    for beta,l in indices[g-1]:
                        hess = flat_partial(grad,beta,l)
                        if hess:
                            add_product(terms,k+l+3,-hess/2,pair(1-alpha,k+1,1-beta,l+1))
            for h in range(g):
                for (alpha,k),grad in grads[h].items():
                    if grad:
                        for (beta,l),other in grads[g-1-h].items():
                            if other:
                                add_product(terms,k+l+3,-grad/2,other,
                                            pair(1-alpha,k+1,1-beta,l+1))
        power = max(terms,default=0)
        cleared = sum((L(terms[p])*delta**(power-p) for p in sorted(terms)),L())
        residuals.append((power,cleared))
    return tuple(residuals)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive",nargs="?",type=Path,
                        default=Path(__file__).with_name("a2-epsilon6.json"))
    args = parser.parse_args(argv)
    residuals = coefficients(read_free_energies(args.archive))
    for g,(p,residual) in enumerate(residuals,start=1):
        print(f"genus {g}, epsilon^{2*g-2}: Delta^{p} times residual has "
              f"{len(residual)} nonzero monomials", flush=True)
    if any(residual for _,residual in residuals):
        raise ArithmeticError("the archived A2 free energies violate the DZ loop equation")


if __name__ == "__main__":
    main()
