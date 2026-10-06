#!/usr/bin/env python3
"""
Route 2: symbolic facts about the configuration-VII fibres E_lambda, lambda = u/w.

Checks (all exact, in sympy):
  1. The Connell model of the fibre v^2 = c33(1, t) (route_c_model.py) has 2-division
     cubic with roots 2*kappa^2 and two further quartics; shifting by 2*kappa^2 gives
         E_lambda:  y^2 = x (x + 4 D A)(x + 4 D B),
     D = u^2+w^2, A = u^2-4uw+w^2, B = u^2+4uw+w^2, kappa = w^2+2uw-u^2.
     So E_lambda is the quadratic twist by D of E': y^2 = x(x+A)(x+B).
  2. Delta = 2^22 u^2 w^2 D^6 A^2 B^2,  c4 = 2^8 D^2 (4D^2 - 3AB),  c6 = 2^12 D^4 (9AB - 8D^2),
     and AB = D^2 - 16 u^2 w^2.
  3. The point x = -12 D^2 satisfies x (x+4DA)(x+4DB) = -3 (16 D^2 (u^2-w^2))^2, i.e.
     P3 = (-12 D^2, 16 sqrt(-3) D^2 (u^2-w^2)) is a section of E_lambda over Q(sqrt(-3))(lambda).
  4. For u, w of opposite parity, c6 / 2^12 = D^4 (9AB - 8D^2) is = 1 mod 8.

    python rootno.py
"""

import sympy as sp

import route_c_model as M

n, d, t, X = M.n, M.d, M.t, sp.symbols('X')
u, w = n, d                     # lambda = u/w; route_c_model calls them n, d
D = u**2 + w**2
A = u**2 - 4*u*w + w**2
B = u**2 + 4*u*w + w**2


def main():
    g, _ = M.build_cells()
    cn = M.connell(sp.expand(g[8].subs({M.p: 1, M.q: t})))
    a1, a2, a3, a4, a6 = [sp.together(x) for x in cn['a']]
    b2, b4, b6 = sp.together(a1**2 + 4*a2), sp.together(2*a4 + a1*a3), sp.together(a3**2 + 4*a6)
    cub = sp.expand(sp.together(4*X**3 + b2*X**2 + 2*b4*X + b6))
    roots = sorted(sp.roots(sp.Poly(cub, X)).keys(), key=str)
    kappa = w**2 + 2*u*w - u**2
    assert any(sp.expand(r - 2*kappa**2) == 0 for r in roots)
    others = [sp.expand(r - 2*kappa**2) for r in roots if sp.expand(r - 2*kappa**2) != 0]
    assert sorted(map(str, others)) == sorted(map(str, [sp.expand(-4*D*A), sp.expand(-4*D*B)])), others
    print("1. E_lambda ~ y^2 = x(x + 4DA)(x + 4DB)  (2-torsion roots of the Connell model) OK")

    # invariants of y^2 = x^3 + a2 x^2 + a4 x
    A2, A4 = 4*D*(A + B), 16*D**2*A*B
    c4 = sp.expand(16*A2**2 - 48*A4)
    c6 = sp.expand(-64*A2**3 + 288*A2*A4)
    Delta = sp.expand(16*A4**2*(A2**2 - 4*A4))
    assert sp.expand(A*B - (D**2 - 16*u**2*w**2)) == 0
    assert sp.expand(Delta - 2**22*u**2*w**2*D**6*A**2*B**2) == 0
    assert sp.expand(c4 - 2**8*D**2*(4*D**2 - 3*A*B)) == 0
    assert sp.expand(c6 - 2**12*D**4*(9*A*B - 8*D**2)) == 0
    print("2. Delta, c4, c6 factorizations OK")

    x0 = -12*D**2
    assert sp.expand(x0*(x0 + 4*D*A)*(x0 + 4*D*B) + 3*(16*D**2*(u**2 - w**2))**2) == 0
    assert sp.expand(x0 + 4*D*A + 8*D*(u + w)**2) == 0 and sp.expand(x0 + 4*D*B + 8*D*(u - w)**2) == 0
    print("3. P3 = (-12 D^2, 16 sqrt(-3) D^2 (u^2 - w^2)) lies on E_lambda OK")

    # 4. residue of c6/2^12 mod 8 for opposite parity: check all (u, w) mod 8
    f = sp.lambdify((u, w), D**4*(9*A*B - 8*D**2))
    bad = [(a, b) for a in range(8) for b in range(8) if (a + b) % 2 == 1 and f(a, b) % 8 != 1]
    assert not bad, bad
    print("4. c6/2^12 = 1 mod 8 whenever u + w is odd OK")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
