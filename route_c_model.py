#!/usr/bin/env python3
"""
Route C model: Bremner's configuration VII as an elliptic fibration, checked
symbolically and written out for the PARI/GP fibre worker.

Configuration VII (Bremner 2001, Acta Arith. 99, section 2, eqs (11)-(13)).
Magic square in Bremner's centre-c form

    [ a+c     -a-b+c   b+c   ]
    [ -a+b+c   c       a-b+c ]
    [ -b+c     a+b+c   -a+c  ]

with the six cells  +-b+c, +-(a-b)+c, -a+c, c  square. Put lambda = n/d and
parametrize the first quadric of (11) by (p : q):

    alpha = n(p^2-q^2) + 2dpq,   delta = d(p^2-q^2) - 2npq,   beta = n, gamma = d
    m = alpha*beta + gamma*delta,  nn = alpha*gamma - beta*delta,
    r = alpha*gamma + beta*delta,  s = alpha*beta - gamma*delta,
    c = m^2 + nn^2 (= r^2 + s^2),  b = 2*m*nn,  a = b + 2*r*s.

Then five cells are squares identically, cell (3,3) = -a+c is the quartic
Q12(p,q), and the fibre over lambda is the genus-one curve  v^2 = Q12(p,q).
The three remaining cells (1,1), (1,2), (3,2) are Bremner's quartics (13).
A seventh square = one of them square at a point of the fibre.

The fibre has the rational point (p:q) = (1:0), v = d^2 + 2nd - n^2 (never 0
over Q). Writing t = q/p, v^2 = A t^4 + B t^3 + C t^2 + D t + q0^2 and the
Connell/Mordell transformation gives a Weierstrass model; both directions of
the map are checked here.

Usage:
    python route_c_model.py            # checks, then writes route_c_cells.gp
    python route_c_model.py --check    # checks only
"""

import os
import sys

import sympy as sp

n, d, p, q, t, v, X, Y = sp.symbols('n d p q t v X Y')

HERE = os.path.dirname(os.path.abspath(__file__))


def build_cells():
    alpha = n * (p**2 - q**2) + 2 * d * p * q
    delta = d * (p**2 - q**2) - 2 * n * p * q
    beta, gamma = n, d
    m = alpha * beta + gamma * delta
    nn = alpha * gamma - beta * delta
    r = alpha * gamma + beta * delta
    s = alpha * beta - gamma * delta
    c = sp.expand(m**2 + nn**2)
    assert sp.expand(c - (r**2 + s**2)) == 0
    b = sp.expand(2 * m * nn)
    a = sp.expand(b + 2 * r * s)
    cells = [a + c, -a - b + c, b + c, -a + b + c, c, a - b + c, -b + c, a + b + c, -a + c]
    return [sp.expand(x) for x in cells], (a, b, c)


def check_magic(g):
    lines = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]
    sums = [sp.expand(g[i] + g[j] + g[k]) for i, j, k in lines]
    assert all(sp.expand(x - sums[0]) == 0 for x in sums), "not magic"
    assert sp.expand(sums[0] - 3 * g[4]) == 0


def check_squares(g):
    # cells 2,3,4,5,6 (0-based) are identically squares; cell 8 is Q12
    roots = {}
    for i in (2, 3, 4, 5, 6):
        f = sp.factor_list(g[i])
        assert all(e % 2 == 0 for _, e in f[1]) and f[0] > 0 and sp.sqrt(f[0]).is_rational, f"cell {i} not a square"
        roots[i] = sp.sqrt(f[0]) * sp.Mul(*[b**(e // 2) for b, e in f[1]])
        assert sp.expand(roots[i]**2 - g[i]) == 0
    return roots


def bremner_quartics():
    lam = n / d
    Q12 = ((1 + 2*lam - lam**2)**2 * p**4 - 32*lam**2 * p**3*q + 2*(1 - 12*lam + 2*lam**2 + 12*lam**3 + lam**4) * p**2*q**2
           + 32*lam**2 * p*q**3 + (1 + 2*lam - lam**2)**2 * q**4)
    Q13a = ((1 - 2*lam - lam**2)**2 * p**4 + 32*lam**2 * p**3*q + 2*(1 + 12*lam + 2*lam**2 - 12*lam**3 + lam**4) * p**2*q**2
            - 32*lam**2 * p*q**3 + (1 - 2*lam - lam**2)**2 * q**4)
    Q13b = ((1 + 2*lam - lam**2)**2 * p**4 - 4*(1 + 10*lam**2 + lam**4) * p**3*q + 2*(1 - 12*lam + 2*lam**2 + 12*lam**3 + lam**4) * p**2*q**2
            + 4*(1 + 10*lam**2 + lam**4) * p*q**3 + (1 + 2*lam - lam**2)**2 * q**4)
    Q13c = ((1 - 2*lam - lam**2)**2 * p**4 + 4*(1 + 10*lam**2 + lam**4) * p**3*q + 2*(1 + 12*lam + 2*lam**2 - 12*lam**3 + lam**4) * p**2*q**2
            - 4*(1 + 10*lam**2 + lam**4) * p*q**3 + (1 - 2*lam - lam**2)**2 * q**4)
    return [sp.expand(x * d**4) for x in (Q12, Q13a, Q13b, Q13c)]


def connell(quart):
    """v^2 = quart(t) with quart(0) = q0^2: Weierstrass coefficients and both maps."""
    P = sp.Poly(quart, t)
    co = [P.coeff_monomial(t**k) for k in range(5)]
    e0, D, C, B, A = co
    q0 = d**2 + 2*n*d - n**2
    assert sp.expand(q0**2 - e0) == 0, "constant term is not (d^2+2nd-n^2)^2"
    a1 = D / q0
    a2 = C - D**2 / (4 * q0**2)
    a3 = 2 * q0 * B
    a4 = -4 * q0**2 * A
    a6 = a2 * a4
    Xf = (2 * q0 * (v + q0) + D * t) / t**2
    Yf = (4 * q0**2 * (v + q0) + 2 * q0 * (D * t + C * t**2) - D**2 * t**2 / (2 * q0)) / t**3
    tf = (2 * q0 * (X + C) - D**2 / (2 * q0)) / Y
    vf = -q0 + tf * (tf * X - D) / (2 * q0)
    return dict(q0=q0, A=A, B=B, C=C, D=D, a=[a1, a2, a3, a4, a6], Xf=Xf, Yf=Yf, tf=tf, vf=vf)


def check_connell(cn, quart, nv, dv):
    """Exact check at a numeric lambda: forward map lands on E, inverse undoes it."""
    sub = {n: nv, d: dv}
    a1, a2, a3, a4, a6 = [x.subs(sub) for x in cn["a"]]  # no nsimplify: it can corrupt exact rationals
    Q = sp.expand(quart.subs(sub))
    # a generic point on the quartic over Q(t, v) with v^2 = Q: check E-equation modulo v^2 - Q
    Xs, Ys = cn['Xf'].subs(sub), cn['Yf'].subs(sub)
    lhs = sp.together(Ys**2 + a1 * Xs * Ys + a3 * Ys - (Xs**3 + a2 * Xs**2 + a4 * Xs + a6))
    num = sp.Poly(sp.expand(sp.numer(lhs)), v)
    red = sp.rem(num, sp.Poly(v**2 - Q, v))
    assert sp.expand(red.as_expr()) == 0, "forward map not on E"
    tb = sp.simplify(cn['tf'].subs(sub).subs({X: Xs, Y: Ys}))
    tb = sp.together(tb)
    numt = sp.Poly(sp.expand(sp.numer(tb) - t * sp.denom(tb)), v)
    assert sp.expand(sp.rem(numt, sp.Poly(v**2 - Q, v)).as_expr()) == 0, "inverse map wrong"


def main(argv):
    g, (a, b, c) = build_cells()
    check_magic(g)
    roots = check_squares(g)
    Qs = bremner_quartics()
    Q12, Q13a, Q13b, Q13c = Qs
    # identify our homogenized cells with Bremner's quartics (up to the d^4 scaling)
    ratio = sp.cancel(g[8] / Q12)
    assert ratio.free_symbols <= {n, d} and not ratio.has(p, q), "cell (3,3) != Q12 up to a constant"
    for cell, Qb, nm in ((g[0], Q13a, '(1,1)'), (g[1], Q13b, '(1,2)'), (g[7], Q13c, '(3,2)')):
        rr = sp.cancel(cell / Qb)
        assert rr == ratio, f"cell {nm} != Bremner quartic times {ratio}"
    print(f"cells: magic OK; cells (1,3),(2,1),(2,2),(2,3),(3,1) squares OK; "
          f"(3,3),(1,1),(1,2),(3,2) = Bremner Q12,Q13a,Q13b,Q13c times {ratio}")
    # Bremner's square at lambda = 13, (p, q) = (9, 2)
    vals = [int(x.subs({n: 13, d: 1, p: 9, q: 2})) for x in g]
    from math import gcd, isqrt
    from functools import reduce
    G = reduce(gcd, [abs(x) for x in vals])
    prim = [x // G for x in vals]
    brem = sorted([373**2, 289**2, 565**2, 360721, 425**2, 23**2, 205**2, 527**2, 222121])
    assert sorted(prim) == brem, "lambda=13 does not give Bremner's square"
    print(f"lambda = 13, (p,q) = (9,2): Bremner's square times {G} OK")
    # Weierstrass model of the fibre, via t = q/p
    quart = sp.expand(g[8].subs({p: 1, q: t}))
    cn = connell(quart)
    for nv, dv in ((13, 1), (3, 7), (-5, 2)):
        check_connell(cn, quart, nv, dv)
    print("Connell model: forward and inverse maps checked at lambda = 13, 3/7, -5/2")
    if '--check' in argv:
        return 0
    write_gp(g, cn)
    return 0


def gp_poly(expr):
    return str(sp.expand(expr)).replace('**', '^')


def write_gp(g, cn):
    out = os.path.join(HERE, 'route_c_cells.gp')
    names = ['c11', 'c12', 'c13', 'c21', 'c22', 'c23', 'c31', 'c32', 'c33']
    with open(out, 'w', newline='\n') as f:
        f.write('/* Generated by route_c_model.py -- do not edit.\n'
                '   Configuration VII cells as forms in (n, d, p, q), lambda = n/d.\n'
                '   Squares identically: c13, c21, c22, c23, c31.  Fibre: v^2 = c33.\n'
                '   Test cells (7th square): c11, c12, c32.  Order: row by row. */\n')
        f.write('rc_cells(n, d, p, q) = [' + ', '.join(gp_poly(x) for x in g) + '];\n')
        quart = sp.expand(g[8].subs({p: 1, q: t}))
        P = sp.Poly(quart, t)
        co = [sp.expand(P.coeff_monomial(t**k)) for k in range(5)]
        f.write('/* c33(1, t) = e0 + D t + C t^2 + B t^3 + A t^4, e0 = q0^2 */\n')
        f.write('rc_quartic(n, d) = [' + ', '.join(gp_poly(x) for x in co) + '];\n')
        f.write('rc_q0(n, d) = ' + gp_poly(d**2 + 2*n*d - n**2) + ';\n')
    print('wrote', out)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
