"""Independent A₂ genus-one and genus-two free-energy comparison.

The four contractions are Eguchi–Getzler–Xiong, hep-th/0007194 §2,
equation (2.5). They use only the dispersionless primary correlators;
the equation-only reconstruction does not consume this reference.
"""

import argparse
import json
import sys
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path


def reference_genera_1_2(leading, *, algebra):
    """Independent Eguchi–Getzler–Xiong §2/(2.5), in DZ A2 flat normalization.

    Let S_ab denote genus-zero primary two-point functions, with lower indices,
    S=[[U,V],[V,3U²/2]]. The inverse of D_x(S) is the symmetric propagator
    C^{ab}. The puncture direction is D_x and the second primary direction is
    E_X0. Thus all primary three- and four-point functions in (2.5) are obtained
    without using the computed deformation or potential.
    """
    a = algebra
    U = 3*a.mono(0)**2
    # V itself need not be a Laurent generator; D_x V and E_X0 V are explicit.
    Ut = (leading[0]-leading[1])/(6*a.mono(0))
    Vt = (leading[0]+leading[1])/2
    Vx = (a.j(0, 1)+a.j(1, 1))/2
    Ux = a.dx(U)
    if Ut != Vx or Vt != 3*U*Ux:
        raise ValueError("unexpected A2 primary t^{2,0} normalization")
    delta = a.j(0, 1)*a.j(1, 1)
    P = ((Ux, Vx), (Vx, 3*U*Ux))
    if P[0][1]**2-P[0][0]*P[1][1] != delta:
        raise ArithmeticError("A2 genus-one determinant identity failed")
    C = ((-P[1][1]/delta, P[0][1]/delta),
         (P[0][1]/delta, -P[0][0]/delta))

    def evolve(expr):
        return a.evol(leading, (expr,))[0]

    @lru_cache(maxsize=None)
    def tri(*t):
        k = sum(t)
        if k == 0:
            return Ux
        if k == 1:
            return Vx
        if k == 2:
            return Vt
        return 3*U*Ut

    @lru_cache(maxsize=None)
    def quad(*t):
        # By commutation of the two primary flows, the full correlator is
        # symmetric. Put a puncture in the derivative slot whenever possible.
        return a.dx(tri(*sorted(t)[1:])) if 0 in t else evolve(tri(1, 1, 1))

    @lru_cache(maxsize=None)
    def dx_P(i, j):
        # ∂²<<O_i O_j>> = ∂(D_x S_ij).
        return a.dx(P[i][j])

    @lru_cache(maxsize=None)
    def dx_tri(*t):
        return a.dx(tri(*t))

    @lru_cache(maxsize=None)
    def dxx_tri(*t):
        return a.dx(tri(*t), 2)

    @lru_cache(maxsize=None)
    def dx_quad(*t):
        return a.dx(quad(*t))

    first = sum((a.dx(quad(i,j,k,l), 2)*C[i][j]*C[k][l]
                 for i,j,k,l in product(range(2), repeat=4)), a.L())/1152
    second = sum((dx_P(i,j)*dx_quad(k,l,m,n)*C[i][k]*C[j][l]*C[m][n]
                  for i,j,k,l,m,n in product(range(2), repeat=6)), a.L())/1152
    third = sum((dxx_tri(i,j,k)*dx_tri(l,m,n)*C[i][l]*C[j][m]*C[k][n]
                 for i,j,k,l,m,n in product(range(2), repeat=6)), a.L())/360
    fourth = sum((dx_P(i,j)*dx_tri(k,l,m)*dx_tri(n,p,q)
                  *C[i][k]*C[j][n]*C[l][p]*C[m][q]
                  for i,j,k,l,m,n,p,q in product(range(2), repeat=8)), a.L())/360
    F1=(a.log(0)+a.log(1))/24
    return F1, first-second-third+fourth


def verify_archive(path, *, algebra):
    """Compare the archived F1 and F2 against the independent exact formulas."""
    a = algebra
    n = len(next(iter(a.mono(0))))
    with Path(path).open(encoding="utf-8") as stream:
        archive = json.load(stream)["free_energy"]

    def decode(row):
        result = {}
        for powers, coefficient in row:
            exponents = [0]*n
            for i, power in powers:
                exponents[i] = power
            key = tuple(exponents)
            if key in result:
                raise ValueError("duplicate sparse monomial")
            result[key] = Q(coefficient)
        return a.L(result)

    checked = {}
    for k, F in zip((2, 4), reference_genera_1_2(a.input_flow(0), algebra=a)):
        recovered = decode(archive[str(k)][0])
        if recovered != F:
            raise ArithmeticError(f"reference mismatch at genus {k//2}")
        checked[k//2] = len(F)
    return checked


def main(argv=None):
    if __package__ in (None, ""):
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.a2 import compute

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path,
                        default=Path(__file__).with_name("a2-epsilon6.json"))
    args = parser.parse_args(argv)
    for genus, count in verify_archive(args.archive, algebra=compute).items():
        print(f"A2 F{genus} agrees with the independent reference: {count} terms")


if __name__ == "__main__":
    main()
