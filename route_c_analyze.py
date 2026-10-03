#!/usr/bin/env python3
"""
Summarize a Route C run directory (FIB lines from route_c_run.py).

    python route_c_analyze.py OUT [OUT2 ...]

Writes OUT/summary.txt and OUT/hits.txt:
  - fibres searched, by status (ok / gens-incomplete / rank-unproven / timeout / error)
  - rank distribution, torsion sizes, total points tested, sanity failures (must be 0)
  - every 7-or-more-square grid, deduplicated up to rotation/reflection (the grids are
    already divided by their largest common square), each re-verified independently:
    8 equal line sums, which cells are squares, all nine distinct.
"""

import ast
import glob
import os
import sys
from collections import Counter
from math import isqrt

LINES = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]


def d4(g):
    m = [g[0:3], g[3:6], g[6:9]]
    out = []
    for _ in range(4):
        out.append(sum(m, []))
        out.append(sum([r[::-1] for r in m], []))
        m = [[m[2 - j][i] for j in range(3)] for i in range(3)]
    return out


def canon(g):
    return min(tuple(x) for x in d4(list(g)))


def issq(x):
    return x >= 0 and isqrt(x) ** 2 == x


def verify(g):
    sums = {sum(g[i] for i in l) for l in LINES}
    return len(sums) == 1, sum(issq(x) for x in g), len(set(g)) == 9


def parse(path):
    for line in open(path):
        if not line.startswith('FIB '):
            continue
        parts = line.rstrip('\n').split(' | ')
        n, d = map(int, parts[0].split()[1:3])
        status = parts[1]
        rk = ast.literal_eval(parts[2])
        ngens, ntors, npts, nbad = (int(x) for x in parts[3:7])
        hits = ast.literal_eval(parts[7]) if len(parts) > 7 else []
        yield n, d, status, rk, ngens, ntors, npts, nbad, hits


def main(dirs):
    seen = {}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, 'fibres_w*.txt'))):
            for rec in parse(f):
                seen[(rec[0], rec[1])] = rec          # last record per lambda wins (re-runs)
    recs = list(seen.values())
    st = Counter(r[2] for r in recs)
    rk = Counter(tuple(r[3]) for r in recs)
    tors = Counter(r[5] for r in recs)
    pts = sum(r[6] for r in recs)
    bad = sum(r[7] for r in recs)
    dmax = max((r[1] for r in recs), default=0)
    hits = {}
    for r in recs:
        for h in r[8]:
            n, d, p, q, nsq, flags, grid = h
            hits.setdefault(canon(grid), []).append((n, d, p, q, nsq, flags))
    out = dirs[0]
    with open(os.path.join(out, 'summary.txt'), 'w') as fh:
        fh.write(f'Route C summary: {len(recs)} fibres (lambda = n/d, 1 <= n < d <= {dmax}), dirs {dirs}\n')
        fh.write(f'status: {dict(st)}\n')
        fh.write(f'rank [lower, upper]: {dict(sorted(rk.items()))}\n')
        fh.write(f'torsion sizes: {dict(tors)}\n')
        fh.write(f'points tested: {pts}\n')
        fh.write(f'sanity failures (c33 not square): {bad}\n')
        fh.write(f'distinct grids with >= 7 squares (up to D4): {len(hits)}\n')
        for g, where in hits.items():
            magic, nsq, distinct = verify(list(g))
            fh.write(f'  grid {list(g)}: magic {magic}, squares {nsq}, distinct {distinct}; '
                     f'found at {len(where)} (lambda, p, q), e.g. {where[:4]}\n')
    with open(os.path.join(out, 'hits.txt'), 'w') as fh:
        for g, where in hits.items():
            fh.write(f'{list(g)}\t{where}\n')
    print(open(os.path.join(out, 'summary.txt')).read())
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
