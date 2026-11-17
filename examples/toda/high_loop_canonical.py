"""Independent exact CP1 loop check with a Laurent spectral variable.

Set y**2=(s-lambda)/(r-lambda) and hold lambda fixed under all jet
derivatives. Compare every Laurent coefficient; no sampling is used.
"""
import argparse
from math import comb
from pathlib import Path
from time import perf_counter

from .high_loop_equation import read_archive


class CanonicalLoop:
    def __init__(self, chart):
        self.chart = chart
        self.d = chart.d
        self.y = chart.V  # reuse this sparse index for the spectral variable
        self.dx_gap = chart.j(0, 1)-chart.j(1, 1)
        self.period = (1-self.y**2)/(self.d*self.y)
        self.w = self.d**2/16

    def add(self, terms):
        result = {}
        for term in terms:
            for key, value in term.items():
                result[key] = result.get(key, 0)+value
        return self.chart.L(result)

    def spectral_variation(self, r, s):
        # Delta y at fixed lambda.
        return self.period*(s-self.y**2*r)/2

    def apply(self, expr, r, s):
        """Contract the canonical jet partials with two sequences of variations."""
        c = self.chart
        gap = r[0]-s[0]
        result = [
            expr.partial(0)*gap,
            expr.partial(c.U_INDEX)*2*gap/self.d,
            expr.partial(c.V_INDEX)*self.spectral_variation(r[0], s[0]),
        ]
        max_order = expr.degree_order()
        for i, row in enumerate((r, s)):
            for n in range(1, max_order+1):
                if part := expr.partial(c.variable_index(i, n)):
                    result.append(part*row[n])
            if part := expr.partial(c.log_index(i)):
                result.append(part*row[1]/c.j(i, 1))
        return self.add(result)

    def dx(self, expr):
        n = expr.degree_order()
        return self.apply(expr, *(tuple(self.chart.j(i, k+1) for k in range(n+1))
                                  for i in (0, 1)))

    def prolong(self, expr, n):
        rows = [expr]
        for _ in range(n):
            rows.append(self.dx(rows[-1]))
        return rows

    def residuals(self, densities):
        """Yield (genus, exact Laurent residual), with no denominators to clear."""
        if not densities:
            raise ValueError("at least F1 is required")
        n = max(f.degree_order() for f in densities)
        previous = max((f.degree_order() for f in densities[:-1]), default=0)
        period = self.prolong(self.period, max(n-1, 1))
        j = self.chart.j

        # Resolvent J and the pushforward A of the differentiated flat period
        # vector D_x (a, -2)/sqrt(Delta).
        j0 = (self.period*self.y, self.period/self.y)
        a0 = (self.period*self.d*j0[0]*j(0, 1)/2,
              -self.period*self.d*j0[1]*j(1, 1)/2)
        a = tuple(self.prolong(x, max(n, previous+1, 1)) for x in a0)
        linear = tuple(
            tuple(self.add([resolvent[k]] + [
                comb(k,m)*period[m-1]*component[k-m] for m in range(1,k+1)])
                  for k in range(n+1))
            for resolvent, component in zip(
                (self.prolong(x, n) for x in j0), a, strict=True))

        if len(densities) > 1:
            # Hessian = E_A^2 - E_(E_A A - B), where B is the second
            # variation of canonical coordinates in a frozen flat direction.
            half_curvature = self.d*period[1]**2/2
            hessian_subtraction = tuple(
                self.prolong(self.apply(x, *a)-sign*half_curvature, previous)
                for x, sign in zip(a0, (1, -1), strict=True))

            # The negative of the RHS connection in canonical coordinates.
            c_r = self.w*j0[0]**3*j0[1]*j(0, 1)
            c_s = self.w*j0[0]*j0[1]**3*j(1, 1)
            chain = self.dx_gap*(c_r-c_s)/(2*self.d)
            connection = tuple(self.prolong(x, previous) for x in
                               (self.dx(c_r)-chain, self.dx(c_s)+chain))

        contractions = []
        for g, density in enumerate(densities):
            lhs = self.apply(density, *linear)
            if g == 0:
                yield 1, lhs+self.w*self.period**4
                continue
            prev = densities[g-1]
            contractions.append(self.apply(prev, *a))
            hess = self.apply(contractions[-1], *a)-self.apply(prev, *hessian_subtraction)
            product = self.add(
                contractions[k]*contractions[g-1-k]*(1 if 2*k == g-1 else 2)
                for k in range((g+1)//2))
            yield g+1, self.add((lhs, (hess+product)/4, self.apply(prev, *connection)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path,
                        default=Path(__file__).with_name("toda-free-energies-g1-g5.json"))
    parser.add_argument("--genus", type=int, help="check only through this genus")
    args = parser.parse_args(argv)
    chart, densities = read_archive(args.archive)
    if args.genus is not None:
        if not 1 <= args.genus <= len(densities):
            parser.error("--genus must be available in the archive")
        densities = densities[:args.genus]
    start = perf_counter()
    for g, residue in CanonicalLoop(chart).residuals(densities):
        print(f"genus {g}, epsilon^{2*g-2}: {len(residue)} nonzero Laurent monomials "
              f"({perf_counter()-start:.2f}s elapsed)", flush=True)
        if residue:
            raise ArithmeticError(f"the archived Toda F{g} violates the DZ loop equation")


if __name__ == "__main__":
    main()
