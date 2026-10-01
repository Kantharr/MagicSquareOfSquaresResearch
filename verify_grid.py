#!/usr/bin/env python3
"""
Check a 3x3 grid: the eight line sums, which entries are perfect squares, and
whether the entries are distinct.

    python verify_grid.py                      # runs the built-in examples
    python verify_grid.py n1 n2 ... n9         # row by row; "x^2" or "x**2" is accepted

Built-in examples:
  bremner  Bremner (1999): magic, 7 of 9 entries square, centre 425^2.
  parker   The "Parker square": 9 squares, 6 distinct, one diagonal fails.

Exit code 0 if the grid is a magic square of 9 distinct squares, 1 otherwise.
"""

import sys
from math import isqrt

LINES = {
    "row 1": (0, 1, 2), "row 2": (3, 4, 5), "row 3": (6, 7, 8),
    "col 1": (0, 3, 6), "col 2": (1, 4, 7), "col 3": (2, 5, 8),
    "diag": (0, 4, 8), "anti": (2, 4, 6),
}

EXAMPLES = {
    "bremner": [373**2, 289**2, 565**2, 360721, 425**2, 23**2, 205**2, 527**2, 222121],
    "parker": [29**2, 1, 47**2, 41**2, 37**2, 1, 23**2, 41**2, 29**2],
}


def is_square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def parse(tok):
    for sep in ("**", "^"):
        if sep in tok:
            b, e = tok.split(sep)
            return int(b) ** int(e)
    return int(tok)


def report(name, g):
    sums = {k: sum(g[i] for i in ix) for k, ix in LINES.items()}
    common = max(set(sums.values()), key=list(sums.values()).count)
    squares = [is_square(x) for x in g]
    distinct = len(set(g)) == 9
    magic = len(set(sums.values())) == 1

    print(f"== {name}")
    for r in range(3):
        cells = []
        for x in g[3 * r:3 * r + 3]:
            cells.append(f"{isqrt(x)}^2" if is_square(x) else str(x))
        print("   " + "  ".join(f"{c:>12}" for c in cells))
    bad = [k for k, v in sums.items() if v != common]
    print(f"   line sums   : {common} on {8 - len(bad)}/8 lines"
          + (f"; off: {', '.join(f'{k}={sums[k]}' for k in bad)}" if bad else ""))
    print(f"   squares     : {sum(squares)}/9")
    print(f"   distinct    : {len(set(g))}/9")
    if magic:
        e = g[4]
        print(f"   centre check: S = 3*centre is {sums['row 1'] == 3 * e}"
              f" (S a square: {is_square(sums['row 1'])})")
    ok = magic and all(squares) and distinct
    print(f"   SOLUTION    : {ok}")
    return ok


def main(argv):
    if len(argv) == 9:
        return 0 if report("input", [parse(t) for t in argv]) else 1
    if argv:
        print(__doc__)
        return 2
    for name, g in EXAMPLES.items():
        report(name, g)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
