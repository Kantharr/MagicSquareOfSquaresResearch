#!/usr/bin/env python3
"""
Route 1: valuations of leg products, and centres with one or two prime factors.

T(E) = { pq : p > q > 0, p^2 + q^2 = E^2 }.  For a prime l = 1 mod 4 with l^a || E, write
p + iq = gamma in Z[i] (norm E^2), l = pi * conj(pi), x = v_pi(gamma). Then
    v_l(pq) = 2 min(x, 2a - x)        if min(x, 2a - x) < a,
    v_l(pq) >= 2a                     if x = a.
For E = l^k the case x = k gives q = 0, so the k elements of T(E) have the distinct
valuations 0, 2, ..., 2k-2 (Theorem A follows from the ultrametric inequality).

    python valuations.py lemma N        # check the lemma against brute force for E <= N
    python valuations.py theoremA N     # no {s, r, r+s}, {s, r+s, r+2s}, {s, r, r+2s}, {r, r+s, r+2s}
                                        #   inside T(l^k) for prime powers l^k <= N
    python valuations.py patterns AMAX  # two primes E = p^a q^b: valuation patterns of
                                        #   (s, r, r+s, r+2s) surviving the ultrametric test
"""

import itertools
import sys
from math import gcd, isqrt

sys.setrecursionlimit(10000)


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def vp(n, p):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def T_brute(E):
    out = []
    for q in range(1, isqrt(E * E // 2) + 1):
        p2 = E * E - q * q
        p = isqrt(p2)
        if p * p == p2 and p > q:
            out.append((p, q))
    return out


def gauss_pi(l):
    """pi = a + bi with a^2 + b^2 = l (l = 1 mod 4)."""
    for a in range(1, isqrt(l) + 1):
        b2 = l - a * a
        b = isqrt(b2)
        if b * b == b2:
            return (a, b)
    raise ValueError(l)


def gdiv_count(g, pi):
    """v_pi(g) for Gaussian integers given as pairs."""
    (x, y), (a, b) = g, pi
    n = a * a + b * b
    k = 0
    while True:
        # g / pi = g * conj(pi) / n
        re, im = x * a + y * b, y * a - x * b
        if re % n or im % n:
            return k
        x, y = re // n, im // n
        k += 1


def check_lemma(N):
    bad = checked = 0
    for E in range(2, N + 1):
        f = factor(E)
        if any(l % 4 == 3 for l in f) or 2 in f:
            continue
        for p, q in T_brute(E):
            for l, a in f.items():
                x = gdiv_count((p, q), gauss_pi(l))
                m = min(x, 2 * a - x)
                v = vp(p * q, l)
                ok = (v == 2 * m) if m < a else (v >= 2 * a)
                checked += 1
                bad += not ok
    print(f'lemma: {checked} (pair, prime) checks for odd E <= {N} with only 1 mod 4 primes, {bad} failures')
    return bad == 0


SHAPES = {'{s,r,r+s}': lambda s, r: (s, r, r + s), '{s,r+s,r+2s}': lambda s, r: (s, r + s, r + 2 * s),
          '{s,r,r+2s}': lambda s, r: (s, r, r + 2 * s), '{r,r+s,r+2s}': lambda s, r: (r, r + s, r + 2 * s)}


def check_theoremA(N):
    found = tested = 0
    for l in range(5, N + 1, 4):
        if len(factor(l)) != 1 or factor(l).get(l) != 1:
            continue
        k, E = 1, l
        while E <= N:
            T = set(p * q for p, q in T_brute(E)) if E <= 20000 else None
            if T is None:
                break
            vals = sorted(vp(t, l) for t in T)
            assert vals == list(range(0, 2 * k, 2)), (E, vals)
            for s in T:
                for r in T:
                    if r == s:
                        continue
                    for nm, f in SHAPES.items():
                        if all(x in T for x in f(s, r)):
                            found += 1
            tested += 1
            k += 1
            E *= l
    print(f'theorem A: {tested} prime powers l^k <= {min(N, 20000)} tested, valuations 0,2,...,2k-2 confirmed, '
          f'{found} three-pair subsets found')
    return found == 0


def ultra_ok(vals):
    """vals: (v_alpha, v_beta, v_gamma) for alpha + beta = gamma; the two smallest must agree."""
    a = sorted(vals)
    return a[0] == a[1]


def patterns(amax):
    """Two primes: each element has class (xp, xq), xp in {0..a-1, 'A'}, value 2*xp or a free
    value >= 2a for 'A'; multiplicity 2 if both classes ordinary, 1 if one is 'A', 0 if both.
    Relations among (s, r, S = r+s, U = r+2s):  s + r = S,  s + S = U,  2S = r + U  (2 is a unit)."""
    survivors = {}
    for a in range(1, amax + 1):
        for b in range(1, amax + 1):
            if (2 * a + 1) * (2 * b + 1) < 9:
                continue
            cls_p = list(range(a)) + ['A']
            cls_q = list(range(b)) + ['A']
            classes = [(x, y) for x in cls_p for y in cls_q if not (x == 'A' and y == 'A')]
            mult = {c: (1 if 'A' in c else 2) for c in classes}
            # free values for 'A': any value >= 2a; only order matters, so try 2a + {0..4}
            def vals(c, k, base, free):
                return [2 * c] if c != 'A' else [2 * base + j for j in range(free)]
            n_ok = 0
            examples = []
            for roles in itertools.product(classes, repeat=4):          # (s, r, S, U)
                cnt = {}
                for c in roles:
                    cnt[c] = cnt.get(c, 0) + 1
                if any(cnt[c] > mult[c] for c in cnt):
                    continue

                def prime_ok(idx, base):
                    choices = [vals(c[idx], None, base, 5) for c in roles]
                    for vs in itertools.product(*choices):
                        s, r, S, U = vs
                        if ultra_ok((s, r, S)) and ultra_ok((s, S, U)) and ultra_ok((S, r, U)):
                            return True
                    return False
                if prime_ok(0, a) and prime_ok(1, b):
                    n_ok += 1
                    if len(examples) < 3:
                        examples.append(roles)
            survivors[(a, b)] = (n_ok, examples)
            print(f'(a,b)=({a},{b}): {n_ok} role-class assignments survive; e.g. {examples[:2]}')
    return survivors


def main(argv):
    if argv[0] == 'lemma':
        return 0 if check_lemma(int(argv[1])) else 1
    if argv[0] == 'theoremA':
        return 0 if check_theoremA(int(argv[1])) else 1
    if argv[0] == 'patterns':
        patterns(int(argv[1]))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
