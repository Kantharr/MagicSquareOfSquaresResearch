# Route C stage 3: find the missing Mordell-Weil generators of incomplete fibres with mwrank.
#
#   sage route_c_mwrank.sage IN.jsonl OUT.jsonl [--limits 12,15] [--procs 10] [--tlim 120]
#
# IN: EXP lines from rc_export() in route_c_fibre.gp (the "EXP " prefix is optional): the exact
# minimal model of the fibre as PARI computed it, PARI's rank bounds and known generators.
# For each fibre, mwrank's two-descent is run on that same model (never changed) with
# second_limit = each value in --limits until it returns rank-many generators, then saturated.
# OUT: one JSON line per fibre:
#   {"n","d","rank","bound","gens":[[x,y],...],"limit","secs","status"}
# status: "complete" (len(gens) == rank == bound), "partial", "timeout" or "error ...".
# Resumable: fibres already in OUT are skipped. Run inside sagemath/sagemath (Docker).

import json, os, sys, time
from multiprocessing import Pool

def parse_args(argv):
    a = {'limits': [12, 15], 'procs': 10, 'tlim': 120}
    pos = []
    i = 0
    while i < len(argv):
        if argv[i] == '--limits':
            a['limits'] = [int(x) for x in argv[i + 1].split(',')]; i += 2
        elif argv[i] == '--procs':
            a['procs'] = int(argv[i + 1]); i += 2
        elif argv[i] == '--tlim':
            a['tlim'] = int(argv[i + 1]); i += 2
        else:
            pos.append(argv[i]); i += 1
    a['inp'], a['out'] = pos[0], pos[1]
    return a

ARGS = None

def work(rec):
    n, d = rec['n'], rec['d']
    res = {'n': n, 'd': d}
    if 'error' in rec:
        res['status'] = 'error export: ' + rec['error']
        return res
    t0 = time.time()
    try:
        alarm(ARGS['tlim'])
        E = EllipticCurve(QQ, [QQ(x) for x in rec['a']])
        best = None
        for L in ARGS['limits']:
            C = E.mwrank_curve()
            C.two_descent(verbose=False, second_limit=L)
            r, b = C.rank(), C.rank_bound()
            g = C.gens()
            best = (L, r, b, g)
            if len(g) >= max(r, rec['rk'][0]) and r == b:
                break
        cancel_alarm()
        L, r, b, g = best
        pts = []
        for P in g:
            P = E(P)                      # points on the same (minimal) model
            pts.append([str(P[0]), str(P[1])])
        res.update(rank=int(r), bound=int(b), gens=pts, limit=int(L),
                   status='complete' if (len(pts) == r == b) else 'partial')
    except AlarmInterrupt:
        res['status'] = 'timeout'
    except Exception as e:
        cancel_alarm()
        res['status'] = 'error ' + type(e).__name__ + ': ' + str(e)[:200]
    res['secs'] = float('%.1f' % (time.time() - t0))   # plain float: Sage preparses literals
    return res

def main():
    global ARGS
    ARGS = parse_args(sys.argv[1:])
    done = set()
    if os.path.exists(ARGS['out']):
        for line in open(ARGS['out']):
            try:
                j = json.loads(line); done.add((j['n'], j['d']))
            except Exception:
                pass
    todo = []
    for line in open(ARGS['inp']):
        line = line.strip()
        if line.startswith('EXP '):
            line = line[4:]
        if not line.startswith('{'):
            continue
        rec = json.loads(line)
        if (rec['n'], rec['d']) not in done:
            todo.append(rec)
    print('mwrank: %d fibres to do, %d already done, limits %s, %d procs, tlim %ds'
          % (len(todo), len(done), ARGS['limits'], ARGS['procs'], ARGS['tlim']), flush=True)
    with open(ARGS['out'], 'a') as fh, Pool(ARGS['procs']) as pool:
        for k, res in enumerate(pool.imap_unordered(work, todo), 1):
            fh.write(json.dumps(res) + '\n'); fh.flush()
            if k % 50 == 0 or k == len(todo):
                print('mwrank: %d/%d' % (k, len(todo)), flush=True)

main()
