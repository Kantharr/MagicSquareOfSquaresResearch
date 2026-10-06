/* Route 2: root numbers of the configuration-VII fibres E_lambda, lambda = u/w.

   E_lambda ~ y^2 = x (x + 4 D A)(x + 4 D B),   D = u^2+w^2, A = u^2-4uw+w^2, B = u^2+4uw+w^2,
   i.e. the quadratic twist by D of E': y^2 = x(x+A)(x+B)  (rootno.py checks this symbolically).
   For odd p the primes dividing u, w, A, B, D are pairwise disjoint (gcd(u,w) = 1), and
       p | uw : w_p = -(-1/p)      p | AB : w_p = -(2/p)      p | D : w_p = +1.
   (At p | A or p | B the reduction is y^2 = x^2 (x + C) with C = +-8uw, and the twist by
    D = -+4uw mod p multiplies the split test by (D/p).)
   rn_check(u, w) compares every local factor with PARI (ellrootno(E, p)) and returns
   [global root number, w_2, number of local mismatches].                                    */

read("route_c_fibre.gp");

rn_model(u, w) =
{
  my(D = u^2 + w^2, A = u^2 - 4*u*w + w^2, B = u^2 + 4*u*w + w^2);
  ellinit([0, 4*D*(A + B), 0, 16*D^2*A*B, 0]);
}

/* predicted local root number at an odd prime p */
rn_local(u, w, p) =
{
  my(D = u^2 + w^2, A = u^2 - 4*u*w + w^2, B = u^2 + 4*u*w + w^2);
  if((u*w) % p == 0, return(-kronecker(-1, p)));
  if(A % p == 0, return(-kronecker(2, p)));
  if(B % p == 0, return(-kronecker(2, p)));
  if(D % p == 0, return(1));
  1;
}

/* odd part of the predicted global root number: prod over odd bad primes */
rn_odd(u, w) =
{
  my(D = u^2 + w^2, A = u^2 - 4*u*w + w^2, B = u^2 + 4*u*w + w^2, P, s = 1);
  P = Set(concat([factor(abs(x))[, 1]~ | x <- [u, w, A, B, D], x != 0]));
  foreach(P, p, if(p > 2, s *= rn_local(u, w, p)));
  s;
}

rn_check(u, w) =
{
  my(E = rn_model(u, w), F = ellinit(rc_curve(u, w)[1]), bad = 0);
  if(!ellisisom(E, F), error("rn_model not isomorphic to rc_curve at ", u, "/", w));
  my(N = ellglobalred(E)[1], P = factor(N)[, 1]);
  foreach(P, p, if(p > 2 && ellrootno(E, p) != rn_local(u, w, p), bad++));
  [ellrootno(E), ellrootno(E, 2), bad];
}

/* closed formula (w_2 = -1 for every coprime u, w, checked numerically by rn_check):
   W(u, w) = (-1)^( #{odd p | uw : p = 1 mod 4} + #{p | AB : p = +-1 mod 8} ),  primes counted once */
rn_W(u, w) =
{
  my(A = u^2 - 4*u*w + w^2, B = u^2 + 4*u*w + w^2, e = 0);
  foreach(Set(concat(factor(abs(u))[, 1]~, factor(abs(w))[, 1]~)), p, if(p % 4 == 1, e++));
  foreach(Set(concat(factor(abs(A))[, 1]~, factor(abs(B))[, 1]~)), p, if(p % 8 == 1 || p % 8 == 7, e++));
  (-1)^e;
}
