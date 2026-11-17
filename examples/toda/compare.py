"""Compare a computed centered Toda inverse with the published CP1 tau gauge."""
from pathlib import Path
import json
from .reference import d, p, q, h, vder, uder, D, E, F1, F2, B, B2, R, can, reference_residuals

def unpack(rows):
    symbols = (d,) + tuple(x for pair in zip(p, q) for x in pair)
    from sympy import prod
    return sum(R(c) * prod(symbols[i]**e for i, e in powers) for powers,c in rows)

def normal_coefficients(dress):
    inv = [tuple(unpack(row) for row in dress[str(k)]) for k in (2,4,6)]
    gap2,gap4,gap6=(can((r-s)/d) for r,s in inv)
    V=tuple(can((r+s)/2) for r,s in inv)
    U=(can(2*gap2),can(2*gap4-gap2**2),can(2*gap6-2*gap2*gap4+2*gap2**3/3))
    # t^1=B V, t^2=B² U, where U's zeroth order is 2 log(d/4).
    w1=tuple(can(V[g-1]+B[g]*vder(2*g)+sum(B[l]*D(V[g-l-1],2*l) for l in range(1,g)))
             for g in (1,2))
    w2=tuple(can(U[g-1]+B2[g]*uder(2*g)+sum(B2[l]*D(U[g-l-1],2*l) for l in range(1,g)))
             for g in (1,2))
    return V,U,w1,w2

def compare(data):
    V,U,w1,w2=normal_coefficients(data["dress"])
    residuals=[]
    for g,(F,a,b) in enumerate(zip((F1,F2), w1, w2),1):
        difference=(can(a-D(E(F))),can(b-D(F,2)))
        print(f"genus {g}: (t^1_(2g) - D E F_g, t^2_(2g) - D² F_g) =", difference)
        residuals.append(difference)
    if residuals != [(0,0),(0,0)] or reference_residuals() != (0,0):
        raise ArithmeticError("literature map did not agree")
    f2_difference = can(unpack(data["free_energy"]["4"][0])-F2)
    if f2_difference != 0:
        raise ArithmeticError("reconstructed F2 differs from (3.10.116)")
    print("Reconstructed F2 − DZ (3.10.116) = 0")
    print("Reconstructed F3:", len(data["free_energy"]["6"][0]), "monomials")
    print("Inverse map V,U number of Laurent monomials:",[[len(a.as_ordered_terms()) for a in row] for row in zip(V,U)])
    return residuals

def verify_flow_direct(dress):
    """Check both centered Toda equations in flat fields (V,W=e^U)."""
    from .reference import A
    inv=[tuple(unpack(row) for row in dress[str(k)]) for k in (2,4,6)]
    V=[None]+[can((r+s)/2) for r,s in inv]
    gap=[d]+[r-s for r,s in inv]
    W=[can(sum(gap[i]*gap[g-i]/16 for i in range(g+1))) for g in range(4)]
    WD=[[D(W[k],2*g+1) for g in range(4-k)] for k in range(4)]
    VD=[[vder(2*g+1) if k==0 else D(V[k],2*g+1) for g in range(4-k)] for k in range(4)]
    for g in range(4):
        Vt = D(W[0]) if g==0 else E(V[g])
        Wt = E(W[g])
        a = can(Vt-sum(A[k]*WD[g-k][k] for k in range(g+1)))
        b = can(Wt-sum(W[i]*A[k]*VD[g-i-k][k] for i in range(g+1) for k in range(g+1-i)))
        if a != 0 or b != 0:
            raise ArithmeticError((g,a,b))
        print(f"centered Toda direct ε^{2*g} residual = (0,0)",flush=True)

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--direct", action="store_true",
                        help="also redo the direct flow verification in SymPy")
    args = parser.parse_args()
    data = json.loads(args.result.read_text(encoding="utf-8"))
    compare(data)
    if args.direct:
        verify_flow_direct(data["dress"])

if __name__ == "__main__":
    main()
