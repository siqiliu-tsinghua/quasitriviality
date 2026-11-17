"""Exact necessary lift certificates for A2 third-derivative perturbations.

Equation: U_t=V_x, V_t=3 U U_x+ε²(A U_xxx+B U_x U_xx+C U_x³).
The calculation is over a real connected domain U>0, at formal differential
weights 3 and 5, with no ε⁴ term in the input flow.

It uses the single explicit dispersionless commuting leader obtained from
the t^{1,1} principal flow. Polynomial lift coefficients can be taken
independent of V at each of these two steps: differentiating the lift
equation in V then invokes the positive-weight homogeneous centralizer
vanishing of SemiHamiltonian/Jet/SymmetryVanishing.lean. The same vanishing
justifies the dilation projection of the ε⁴ equation.
"""

from functools import cache

import sympy as sp


U = sp.Symbol("U", positive=True)
V = sp.Symbol("V")
A, B, C = (sp.Function(name)(U) for name in ("A", "B", "C"))
p = sp.symbols("p1:9")
q = sp.symbols("q1:9")
p1, p2, p3, p4 = p[:4]
q1, q2, q3, q4 = q[:4]
X0 = (q1, 3*U*p1)
Y0 = (V*p1+U*q1, V*q1+3*U**2*p1)


@cache
def total_derivative(expr):
    """Apply D_x to the finite-state and finite-jet expression."""
    if expr.has(p[-1], q[-1]):
        raise ValueError("extend the certificate's jet symbols before differentiating")
    ans = expr.diff(U)*p1 + expr.diff(V)*q1
    for jets in (p, q):
        ans += sum(expr.diff(jets[k])*jets[k+1]
                   for k in range(len(jets)-1) if expr.has(jets[k]))
    return sp.expand(ans)


def dx(expr, n):
    for _ in range(n):
        expr = total_derivative(expr)
    return expr


def evolve(vector, target):
    result = []
    for value in target:
        ans = value.diff(U)*vector[0] + value.diff(V)*vector[1]
        for i, jets in enumerate((p, q)):
            ans += sum(value.diff(jet)*dx(vector[i], k+1)
                       for k, jet in enumerate(jets) if value.has(jet))
        result.append(sp.expand(ans))
    return tuple(result)


def bracket(left, right):
    """Package convention: [left,right] = E_left right - E_right left."""
    return tuple(sp.expand(a-b) for a, b in zip(evolve(left, right), evolve(right, left)))


@cache
def jet_monomials(weight):
    """All polynomial jet monomials of the specified differential weight."""
    if weight == 0:
        return (sp.S.One,)
    return tuple(sorted({j*m for k in range(1, weight+1)
                         for j in (p[k-1], q[k-1])
                         for m in jet_monomials(weight-k)}, key=sp.default_sort_key))


def second_order_witness(a, b):
    """The explicit surviving ε² symmetry correction, with A=a, B=b, C=0."""
    return (b*(p1*q2+p2*q1)/3 + 2*a*q3/3,
            (2*U*b+6*a)*p1*p2 + 3*U*a*p3 + b*p1**3 + b*q1*q2/3)


def second_order_conditions():
    """Derive and check that the ε² lift requires A'=B'=0 and C=0."""
    if bracket(X0, Y0) != (0, 0):
        raise ArithmeticError("the selected leading field must commute")
    X2 = (sp.S.Zero, A*p3+B*p1*p2+C*p1**3)
    rhs = bracket(Y0, X2)
    if any(term.has(V) for term in rhs):
        raise ArithmeticError("unexpected V-dependent ε² source")
    basis = jet_monomials(3)
    u = tuple(tuple(sp.Function(f"z{i}_{k}")(U) for k in range(len(basis)))
              for i in range(2))
    value = u[0]+u[1]
    derivative = tuple(z.diff(U) for z in value)
    unknown = value+derivative
    lift = tuple(sum(z*m for z, m in zip(row, basis)) for row in u)
    equations = tuple(c for lhs, r in zip(bracket(X0, lift), rhs)
                      for c in sp.Poly(lhs-r, *p[:4], *q[:4]).coeffs())
    linear = sp.linsolve(equations, unknown)
    solved = dict(zip(unknown, next(iter(linear))))
    free = [key for key, val in solved.items() if key == val]
    # These four free matrix coordinates are derivatives of coefficients
    # already determined algebraically; impose this exact jet consistency.
    if any(key not in derivative or solved[key.expr].has(*free) for key in free):
        raise ArithmeticError("cannot determine the remaining coefficient derivatives")
    forced = {key: solved[key.expr].diff(U) for key in free}
    solved = {key: sp.expand(val.xreplace(forced)) for key, val in solved.items()}
    consistency = tuple(sp.factor(solved[d.expr].diff(U)-solved[d]) for d in derivative)
    if not any(eq == 2*A.diff(U)/3 for eq in consistency):
        raise ArithmeticError("missing explicit A' condition")
    if not any(eq == 2*B.diff(U)/3 for eq in consistency):
        raise ArithmeticError("missing explicit B' condition")
    constant_tests = {d: 0 for eq in consistency for d in eq.atoms(sp.Derivative)
                      if d.expr in (A, B)}
    residuals = [eq.xreplace(constant_tests) for eq in consistency]
    if sp.linsolve(residuals, (C, C.diff(U))) != {(0, 0)}:
        raise ArithmeticError("C was not forced to vanish")
    a, b = sp.symbols("a b")
    Y2 = second_order_witness(a, b)
    if bracket(X0, Y2) != bracket(Y0, (sp.S.Zero, a*p3+b*p1*p2)):
        raise ArithmeticError("explicit surviving ε² correction did not verify")
    return A.diff(U), B.diff(U), C


def _state_weight(monomial):
    """Scaling weight, not differential weight: U,jets=2; V,jets=3."""
    exponents = monomial.as_powers_dict()
    return 2*sum(exponents[j] for j in p) + 3*sum(exponents[j] for j in q)


def fourth_order_certificate():
    """Return the exact obstruction -2/3 for the normalized B² branch.

    Under the weighted dilation ad_Z for Z=(2U,3V), the a², ab, b² parts
    of [Y2,X2] have distinct eigenweights 0,2,4. A putative lift projects
    to one solving the b² branch, of dilation eigenweight 3, since X0 has
    weight 1 and its positive differential-weight centralizer vanishes.
    Each coefficient in the V-independent homogeneous lift is therefore a
    constant times a determined power of U. After substituting this form,
    evaluating the coefficient equations at U=1 is exact: each equation has
    an overall single power of U. The linear combination below cancels
    every unspecified constant, using only 15 jet coefficient equations.
    """
    basis = jet_monomials(5)
    x = (sp.symbols(f"z0_0:{len(basis)}"), sp.symbols(f"z1_0:{len(basis)}"))
    lift = tuple(sum(c*m*U**((sp.Integer(3)+w-_state_weight(m))/2)
                     for c, m in zip(row, basis)) for w, row in zip((2, 3), x))
    a, b = sp.symbols("a b")
    rhs_full = bracket(second_order_witness(a, b), (sp.S.Zero, a*p3+b*p1*p2))
    dilation = (2*U, 3*V)
    if bracket(dilation, X0) != X0:
        raise ArithmeticError("the dispersionless field must have dilation weight one")
    rhs = None
    for powers, eigen in (((2, 0), 0), ((1, 1), 2), ((0, 2), 4)):
        piece = tuple(sp.Poly(f, a, b).coeff_monomial(powers) for f in rhs_full)
        if bracket(dilation, piece) != tuple(eigen*x for x in piece):
            raise ArithmeticError("unexpected dilation grade of ε⁴ source")
        if powers == (0, 2):
            rhs = piece
    rows = tuple(sp.Poly(sp.expand(lhs-r).subs(U, 1), *p[:6], *q[:6])
                 for lhs, r in zip(bracket(X0, lift), rhs))
    certificate = (
        (-1, 0, p1**3*p3), (-6, 0, p1*q1*q4), (9, 0, p1*q2*q3),
        (sp.Rational(2, 3), 0, p2**3), (-3, 0, p2*q1*q3),
        (-10, 0, p2*q2**2), (9, 0, p3*q1*q2), (-10, 0, p4*q1**2),
        (-24, 0, q1**3*q3), (-3, 0, q3**2), (-1, 1, p1**3*q3),
        (1, 1, p1*p2*q3), (-sp.Rational(4, 3), 1, p2**2*q2),
        (-8, 1, p3*q1**3), (8, 1, q2**3),
    )
    obstruction = sp.expand(sum(c*rows[i].coeff_monomial(m) for c, i, m in certificate))
    if obstruction != -sp.Rational(2, 3):
        raise ArithmeticError("the exact ε⁴ residual identity failed")
    return obstruction


def verify_certificate():
    """Verify the ε² differential and ε⁴ algebraic obstructions."""
    return second_order_conditions(), fourth_order_certificate()


if __name__ == "__main__":
    conditions, obstruction = verify_certificate()
    print("ε²: " + ", ".join(f"{x} = 0" for x in conditions))
    print(f"ε⁴: {obstruction} * B**2 = 0")
