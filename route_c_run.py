#!/usr/bin/env python3
"""
Route C runner: search the fibres of configuration VII for 7-square magic squares.

For every lambda = n/d with 1 <= n < d, gcd(n, d) = 1 and DMIN <= d <= DMAX
(the symmetries lambda -> -lambda, 1/lambda, -1/lambda give the same grids up to
rotation, see route_c_model.py notes), run rc_line(n, d, HMAX) from
route_c_fibre.gp and append its FIB line to OUT/fibres_w<k>.txt.

    python route_c_run.py OUT DMIN DMAX HMAX [--workers 10] [--batch 100]
                          [--deadline "2026-10-04 22:00"] [--wait-gate]

Resumable: lambdas that already have a FIB line in OUT are skipped.
CPU gates (the project rule): before every batch the script asks
~/.claude/hooks/cpu-gate.ps1 whether a heavy job may run (weekday 8-17 Central,
Perfect Cuboid jobs on the CPU). If not, it pauses and re-checks every 10 minutes
(--wait-gate), or stops. Workers stop taking batches after --deadline (Central).
gp is started on a batch file inside OUT, so its command line names this project
and the gate does not mistake our workers for cuboid jobs.
"""

import argparse
import datetime as dt
import glob
import json
import math
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
GP = os.environ.get('PARI_GP_PATH', r'C:\Users\Owner\AppData\Local\Programs\PariGP\gp.exe')
GATE = os.path.expanduser(r'~\.claude\hooks\cpu-gate.ps1')
CENTRAL = dt.timezone(dt.timedelta(hours=-5))   # CDT (valid until 2026-11-01); the gate uses the real zone


def gate_ok():
    """True if the CPU gate allows a heavy job now. Returns (ok, reason)."""
    if not os.path.exists(GATE):
        return True, 'no gate script'
    payload = json.dumps({'tool_name': 'Bash', 'tool_input': {'command': 'docker run route-c-gate-probe'}})
    r = subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', GATE],
                       input=payload, capture_output=True, text=True, timeout=120)
    out = r.stdout.strip()
    if '"deny"' in out:
        m = re.search(r'"permissionDecisionReason":"([^"]*)"', out)
        return False, m.group(1) if m else out
    return True, ''


def done_set(out):
    done = set()
    for f in glob.glob(os.path.join(out, 'fibres_w*.txt')):
        with open(f) as fh:
            for line in fh:
                if line.startswith('FIB '):
                    n, d = line.split()[1:3]
                    done.add((int(n), int(d)))
    return done


def run_batch(out, k, batch, hmax):
    src = os.path.join(out, f'batch_w{k}.gp').replace('\\', '/')
    with open(src, 'w', newline='\n') as f:
        f.write('default(parisize, 200000000); default(parisizemax, 2000000000);\n')
        f.write(f'read("{HERE}/route_c_fibre.gp");\n')
        for n, d in batch:
            f.write(f'rc_line({n}, {d}, {hmax});\n')
        f.write('quit\n')
    r = subprocess.run([GP, '-q', src], capture_output=True, text=True, cwd=HERE)
    lines = [l for l in r.stdout.splitlines() if l.startswith('FIB ')]
    with open(os.path.join(out, f'fibres_w{k}.txt'), 'a') as fh:
        for l in lines:
            fh.write(l + '\n')
    got = {(int(l.split()[1]), int(l.split()[2])) for l in lines}
    missing = [x for x in batch if x not in got]
    if missing:
        with open(os.path.join(out, 'errors.txt'), 'a') as fh:
            fh.write(f'{dt.datetime.now().isoformat()} worker {k}: {len(missing)} lambdas without output, '
                     f'first {missing[:3]}; stderr: {r.stderr[-500:]!r}\n')
    return len(lines), sum(1 for l in lines if not l.rstrip().endswith('| []'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('dmin', type=int)
    ap.add_argument('dmax', type=int)
    ap.add_argument('hmax', type=float)
    ap.add_argument('--workers', type=int, default=10)
    ap.add_argument('--batch', type=int, default=100)
    ap.add_argument('--deadline', default='2026-10-04 22:00')
    ap.add_argument('--wait-gate', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    deadline = dt.datetime.strptime(a.deadline, '%Y-%m-%d %H:%M').replace(tzinfo=CENTRAL)

    todo = [(n, d) for d in range(max(2, a.dmin), a.dmax + 1) for n in range(1, d) if math.gcd(n, d) == 1]
    done = done_set(a.out)
    todo = [x for x in todo if x not in done]
    # interleave so that every batch mixes small and large d (balanced batch times)
    batches = [todo[i::max(1, math.ceil(len(todo) / a.batch))] for i in range(max(1, math.ceil(len(todo) / a.batch)))]
    log = open(os.path.join(a.out, 'progress.txt'), 'a', buffering=1)
    log.write(f'{dt.datetime.now().isoformat()} start: d in [{a.dmin},{a.dmax}], hmax {a.hmax}, '
              f'{len(todo)} lambdas to do ({len(done)} already done), {len(batches)} batches, '
              f'{a.workers} workers, deadline {a.deadline} Central\n')
    lock = threading.Lock()
    state = {'next': 0, 'fibres': 0, 'hits': 0, 'stop': False}

    def take():
        with lock:
            if state['stop'] or state['next'] >= len(batches):
                return None
            i = state['next']
            state['next'] += 1
            return i

    def gate_wait():
        while True:
            if dt.datetime.now(CENTRAL) >= deadline:
                return False
            ok, why = gate_ok()
            if ok:
                return True
            log.write(f'{dt.datetime.now().isoformat()} gate closed: {why}\n')
            if not a.wait_gate:
                return False
            time.sleep(600)

    def worker(k):
        while True:
            if not gate_wait():
                with lock:
                    state['stop'] = True
                return
            i = take()
            if i is None:
                return
            t0 = time.time()
            nf, nh = run_batch(a.out, k, batches[i], a.hmax)
            with lock:
                state['fibres'] += nf
                state['hits'] += nh
                log.write(f'{dt.datetime.now().isoformat()} w{k} batch {i + 1}/{len(batches)}: {nf} fibres, '
                          f'{nh} with hits, {time.time() - t0:.0f}s; total {state["fibres"]}/{len(todo)}, '
                          f'hit lines {state["hits"]}\n')

    with ThreadPoolExecutor(a.workers) as ex:
        list(ex.map(worker, range(a.workers)))
    log.write(f'{dt.datetime.now().isoformat()} finished: {state["fibres"]} fibres this session, '
              f'stopped early: {state["stop"]}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
