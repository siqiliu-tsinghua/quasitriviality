"""Check the full maximal-jet slices of the stored Toda dressing map."""
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
import sys


def verify_top_jets(archive):
    bound = archive["basis"]["jet_order"]
    actual = {}
    for g in range(1, archive["cutoff"]//2 + 1):
        order = 3*g
        high = (2*order-1, 2*order)
        result = {}
        for sign, row in zip((1, -1), archive["coefficients"][str(2*g)]["dress"], strict=True):
            for exponents, coefficient in row:
                key = dict(exponents)
                if any(high[1] < i <= 2*bound for i in key):
                    raise ArithmeticError(f"dressing order exceeds {order} at genus {g}")
                if high[0] not in key and high[1] not in key:
                    continue
                # U = 2 log(d/4). Products of lower-genus corrections
                # do not contribute at jet order 3g.
                key[0] = key.get(0, 0) - 1
                monomial = tuple(sorted((i,e) for i, e in key.items() if e))
                result[monomial] = result.get(monomial, 0) + 2*sign*Fraction(coefficient)
        weight = Fraction(1, 2**(g-1)*24**g*factorial(g))
        expected = {
            tuple(([(0, g-1)] if g>1 else []) + [(i+1, -g), (h, 1)]):
            weight*(-1 if i==1 and g%2==0 else 1)
            for i,h in enumerate(high)
        }
        actual[g] = {k: v for k, v in result.items() if v}
        if actual[g] != expected:
            raise ArithmeticError(f"unexpected top-jet slice at genus {g}")
    return actual


if __name__ == "__main__":
    verified = verify_top_jets(json.loads(Path(sys.argv[1]).read_text(encoding="utf-8")))
    for g, row in verified.items():
        print(f"genus {g}: {row}")
