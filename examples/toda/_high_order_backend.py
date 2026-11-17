"""Private exact sparse chart for the high-order Toda benchmark.

It reuses the base benchmark's transport and differential operations. Larger
jet capacity and packed arithmetic are local: lower-order examples continue
to use their original tuple indices and cached jet symbols.
"""
from functools import lru_cache
from math import lcm
from types import FunctionType, SimpleNamespace

from . import compute


def _copy(function, scope):
    return FunctionType(function.__code__, scope, function.__name__,
                        function.__defaults__, function.__closure__)


def make_chart(max_jet_order=20):
    scope = vars(compute).copy()
    for name, item in vars(compute).items():
        bare = getattr(item, "__wrapped__", item)
        if isinstance(bare, FunctionType) and bare.__globals__ is vars(compute):
            local = _copy(bare, scope)
            scope[name] = (lru_cache(**item.cache_parameters())(local)
                           if bare is not item else local)

    # Signed exponents packed in balanced base 2**32. Every exponent in
    # this bounded jet calculation is many orders of magnitude smaller.
    bits = 32
    half = 1 << (bits - 1)
    mask = 2*half - 1
    count = 5 + 2*max_jet_order
    shifts = tuple(bits*i for i in range(count))
    strides = tuple(1 << s for s in shifts)
    bias = half*sum(strides)
    Q = compute.Q

    def power_at(key, i):
        return (((key + bias) >> shifts[i]) & mask) - half

    def move_power(key, i, increment):
        return key + strides[i]*increment

    def exponents(key):
        return tuple(power_at(key, i) for i in range(count))

    def add(a, b):
        if not isinstance(b, L):
            b = L(b)
        out = L()
        out.update(a)
        _accumulate(out, b)
        return out

    def mul(a, b):
        out = L()
        if isinstance(b, (int, Q)):
            if b:
                out.update((k, v*b) for k, v in a.items())
            return out
        if not a or not b:
            return out
        if len(a) > len(b):
            a, b = b, a
        if len(a) == 1:
            (k, v), = a.items()
            out.update((k+j, w if v == 1 else v*w) for j, w in b.items())
            return out
        ad = lcm(*(v.denominator for v in a.values()))
        bd = lcm(*(v.denominator for v in b.values()))
        ai = [(k, v.numerator*(ad//v.denominator)) for k, v in a.items()]
        bi = [(k, v.numerator*(bd//v.denominator)) for k, v in b.items()]
        acc = {}
        for k, v in ai:
            for j, w in bi:
                i = k+j
                acc[i] = acc.get(i, 0) + v*w
        out.update((k, Q(v, ad*bd)) for k, v in acc.items() if v)
        return out

    def power(a, n):
        if len(a) == 1:
            (k, v), = a.items()
            return L({k*n: v**n})
        return base_power(a, n)

    def substitute(a, func):
        if not a:
            return L()
        mono, multi = [], []
        for i in a.active():
            mapped = func(i)
            if len(mapped) == 1:
                (key, coefficient), = mapped.items()
                mono.append((i, key, coefficient, {}))
            else:
                multi.append((i, mapped, {0: L(1), 1: mapped}))
        # Group the sparse higher-order jet powers before expanding the
        # lower-order factors. Monomial substitutions need no expansion.
        multi.reverse()
        rows = []
        for key, value in a.items():
            outkey = 0
            for i, base, coefficient, cache in mono:
                exponent = power_at(key, i)
                if exponent:
                    outkey += base*exponent
                    if coefficient != 1:
                        if exponent not in cache:
                            cache[exponent] = coefficient**exponent
                        value *= cache[exponent]
            rows.append((tuple(power_at(key, i) for i, _, _ in multi), outkey, value))

        def expand(group, depth):
            out = L()
            if depth == len(multi):
                _accumulate(out, ((key, value) for _, key, value in group))
                return out
            _, base, powers = multi[depth]
            groups = {}
            for row in group:
                groups.setdefault(row[0][depth], []).append(row)
            for exponent, block in groups.items():
                sub = expand(block, depth+1)
                if exponent:
                    if exponent not in powers:
                        powers[exponent] = base**exponent
                    sub = sub*powers[exponent]
                _accumulate(out, sub)
            return out

        return expand(rows, 0)

    def serial(a):
        return [[[[i, e] for i, e in enumerate(key) if e], str(v)]
                for key, v in sorted((exponents(k), v) for k, v in a.items())]

    def firstpart(row, rhs, weight):
        # Same own-first-jet primitive as compute.firstpart. Only its
        # state-only key test differs (packed instead of sliced tuple).
        ind = scope["variable_index"](row, 1)
        ell = scope["log_index"](row)
        slope = scope["dstate"](scope["vspeed"](row), row)
        coefs = {}
        factorial = compute.factorial
        for key, value in rhs.items():
            if power_at(key, ind) != -weight-1:
                continue
            k = power_at(key, ell)
            base = move_power(move_power(key, ind, weight+1), ell, -k)
            if abs(base) >= half:
                continue
            if weight == 0:
                p = move_power(base, ell, k+1)
                coefs[p] = coefs.get(p, 0) - value/(k+1)
            else:
                for h in range(k+1):
                    p = move_power(move_power(base, ind, -weight), ell, k-h)
                    coefs[p] = coefs.get(p, 0) + value*Q(
                        factorial(k), factorial(k-h)*weight**(h+1))
        return L(coefs)/slope

    def _accumulate(out, incoming):
        if isinstance(incoming, dict):
            incoming = incoming.items()
        for key, value in incoming:
            if value := out.get(key, 0) + value:
                out[key] = value
            else:
                out.pop(key, None)

    methods = {name: _copy(item, scope) for name, item in vars(compute.L).items()
               if isinstance(item, FunctionType)}
    base_power = methods["__pow__"]
    methods.update(__add__=add, __radd__=add, __mul__=mul, __rmul__=mul,
                   __pow__=power, substitute=substitute, serial=serial)
    L = type("PackedTodaLaurent", (compute.L,), methods)
    scope.update(L=L, MAX_JET_ORDER=max_jet_order, N=max_jet_order,
                 V_INDEX=count-2, U_INDEX=count-1, VARIABLE_COUNT=count,
                 ZERO_KEY=0, power_at=power_at, move_power=move_power,
                 product_key=int.__add__, quotient_key=int.__sub__, firstpart=firstpart)
    scope["d"], scope["V"], scope["U"] = (
        scope["mono"](i) for i in (0, count-2, count-1))
    return SimpleNamespace(**scope)
