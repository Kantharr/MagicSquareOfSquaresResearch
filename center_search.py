#!/usr/bin/env python3
"""
Search for 3x3 magic squares of squares by centre E (centre cell E^2).

    python center_search.py N            # every centre E <= N
    python center_search.py N --quiet    # summary only

Reduction (research_notes.md, section 2). In any 3x3 magic square the magic
sum is 3e and the four lines through the centre give a + i = b + h = c + g =
d + f = 2e. With centre E^2 each such pair is a 3-term arithmetic progression
of squares X^2, E^2, Y^2, i.e. X^2 + Y^2 = 2E^2. Writing X = p - q, Y = p + q
gives p^2 + q^2 = E^2, and the common difference is E^2 - X^2 = 2pq. So let

    T(E) = { p*q : p > q > 0, p^2 + q^2 = E^2 }   (legs of right triangles with hypotenuse E)

A magic square of 9 distinct squares with centre E^2 exists exactly when
T(E) contains s, r, r+s, r+2s for some r, s > 0 with r != s. The grid is then

    E^2 + 2s        E^2 + 2r      E^2 - 2(r+s)
    E^2 - 2(r+2s)   E^2           E^2 + 2(r+2s)
    E^2 + 2(r+s)    E^2 - 2r      E^2 - 2s

If only 3 of the 4 numbers lie in T(E), the grid is a magic square with 7
square entries: the centre and three complete pairs. (Bremner's 1999 example,
centre 425^2, is a different 7-square shape: two complete pairs plus one
square from each of the other two pairs, so this search does not find it.)
The search reports every primitive 4-of-4 hit (a solution) and every
primitive 3-of-4 hit (a 7-square near miss of the three-pair shape).

Pure Python, no dependencies. Memory and time grow like N log N.
"""

import sys
from math import gcd, isqrt


def triangle_products(N):
    """T[E] for every E <= N, from Euclid's parametrisation plus scaling."""
    T = [None] * (N + 1)
    m = 2
    while m * m + 1 <= N:
        for n in range(1, m):
            if (m - n) % 2 == 0 or gcd(m, n) != 1:
                continue
            h = m * m + n * n
            if h > N:
                break
            prim = (m * m - n * n) * (2 * m * n)
            k = 1
            while k * h <= N:
                E = k * h
                if T[E] is None:
                    T[E] = set()
                T[E].add(k * k * prim)
                k += 1
        m += 1
    return T


def candidates(ts):
    """All (r, s), r != s, with at least 3 of s, r, r+s, r+2s in ts."""
    out = {}
    vals = sorted(ts)
    for x in vals:
        for y in vals:
            if x == y:
                continue
            # x, y take two of the four roles; recover (r, s) from each choice
            trial = [(y, x), (y - x, x), (y - 2 * x, x), (x, y - x), (2 * x - y, y - x)]
            if (y - x) % 2 == 0:
                trial.append((x, (y - x) // 2))
            for r, s in trial:
                if r <= 0 or s <= 0 or r == s or (r, s) in out:
                    continue
                four = (s, r, r + s, r + 2 * s)
                present = sum(t in ts for t in four)
                if present >= 3:
                    out[(r, s)] = (present, tuple(t for t in four if t not in ts))
    return out


def grid(E, r, s):
    e = E * E
    return [e + 2 * s, e + 2 * r, e - 2 * (r + s),
            e - 2 * (r + 2 * s), e, e + 2 * (r + 2 * s),
            e + 2 * (r + s), e - 2 * r, e - 2 * s]


def is_square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def primitive(g):
    d = 0
    for x in g:
        d = gcd(d, x)
    # a scaled copy of a smaller square has some k^2 > 1 dividing every entry
    k = 2
    while k * k <= d:
        if d % (k * k) == 0:
            return False
        k += 1
    return True


def fmt(x):
    return f"{isqrt(x)}^2" if is_square(x) else str(x)


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    N = int(argv[0])
    quiet = "--quiet" in argv
    T = triangle_products(N)
    full, near = [], []
    for E in range(1, N + 1):
        ts = T[E]
        if ts is None or len(ts) < 2:
            continue
        for (r, s), (present, missing) in sorted(candidates(ts).items()):
            g = grid(E, r, s)
            if not primitive(g):
                continue
            rec = (E, r, s, g, missing)
            (full if present == 4 else near).append(rec)
            if not quiet:
                tag = "SOLUTION" if present == 4 else "7-square"
                pos = "" if min(g) > 0 else "  (has entry <= 0)"
                print(f"{tag}  E={E}  r={r}  s={s}  missing t={list(missing)}{pos}")
                for row in range(3):
                    print("    " + "  ".join(f"{fmt(x):>14}" for x in g[3 * row:3 * row + 3]))
    with_T = sum(1 for E in range(1, N + 1) if T[E] and len(T[E]) >= 2)
    print(f"\nE <= {N}: {with_T} centres with |T(E)| >= 2; "
          f"primitive solutions: {len(full)}; primitive 7-square near misses: {len(near)} "
          f"({sum(1 for x in near if min(x[3]) > 0)} with all entries positive)")
    return 0 if full else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
