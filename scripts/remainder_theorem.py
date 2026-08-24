#!/usr/bin/env python3
"""remainder_theorem.py -- an EXPLICIT, rigorous bound on the O(x^6) of the Wilson-line potential.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT IS PROVED, AND WHAT IS ONLY CHECKED
-----------------------------------------
Proved on paper (the derivation is below, every step elementary).  Checked here: that the
bound HOLDS on a wide family of contents, that it is not vacuous, and that it correctly
refuses to apply outside its own hypothesis.  A bound nobody tries to break is a wish.

THE STATEMENT
-------------
Let the bulk content be a set of terms (c, s, m): integer charge c >= 1, sign s = +-1,
multiplicity m >= 1.  Write

    N = sum m ,   c_max = max c ,   x = pi alpha ,   u = c_max * alpha ,
    A_{2k} = sum_{s=+1} m c^{2k} ,  B_{2k} = sum_{s=-1} m c^{2k} ,
    D = A_2 - (3/4) B_2 ,  A_4 ,  G = A_4 H_4 - sum_{s=+1} m c^4 ln c - ln2 * B_4 ,  H_4 = 25/12 .

Then for every alpha with u < 1,

    F(alpha) = F(0) - (zeta(3) D / 2) x^2 + (x^4 / 24) (G - A_4 ln x) + R(alpha)

with

    |R(alpha)|  <=  pi^6 N u^6 / ( 2160 (1 - u^2) )                                 (T1)
    |R'(alpha)| <=  pi^6 N c_max u^5 / ( 360 (1 - u^2) )  (derivative in alpha)     (T2)

THE DERIVATION
--------------
1.  Term by term,  Re Li_5(s e^{i c x}) = sum_{k>=0} (-1)^k (c x)^{2k} / (2k)! * L_k(s),
    with L_k(+1) = zeta(5-2k) and L_k(-1) = -eta(5-2k), both by analytic continuation.
    k = 0,1,2 are the three rungs the paper keeps; k = 2 is where the pole of zeta puts the
    x^4 ln x.  So R is the sum over k >= 3.

2.  eta(-m) = (1 - 2^{m+1}) zeta(-m) with m = 2k-5 gives
        A_{2k} L_k(+1) + B_{2k} L_k(-1) = zeta(5-2k) [ A_{2k} + (2^{2k-4} - 1) B_{2k} ] .
    (At k = 3 that bracket is A_6 + 3 B_6; at k = 4, A_8 + 15 B_8.  It is equation (6) of the
    paper at p = 5-2k, and nothing here is new -- only the bound is.)

3.  |A_{2k} + (2^{2k-4}-1) B_{2k}| <= c_max^{2k} * 2^{2k-4} * N ,  since 2^{2k-4} >= 1 for k >= 2.

4.  zeta(5-2k) = -B_{2k-4} / (2k-4), and |B_{2n}| = 2 (2n)! zeta(2n) / (2pi)^{2n}
    <= (pi^2/3) (2n)! / (2pi)^{2n}, using zeta(2n) <= zeta(2) = pi^2/6 for n >= 1.

5.  Collecting, with u = c_max x / pi = c_max alpha,

        |T_k| <= (pi^6/3) N u^{2k} / [ (2k-4)(2k-3)(2k-2)(2k-1)(2k) ] .

    The bracket is minimal at k = 3, where it is 2*3*4*5*6 = 720, and it only grows.  So

        |R| <= (pi^6/3) N (u^6/720) sum_{j>=0} u^{2j} = pi^6 N u^6 / (2160 (1-u^2)) ,

    which is (T1).  Differentiating termwise (legitimate inside the radius) and using
    d/d(alpha) = pi d/dx, the same collection with one factor 2k pulled out gives (T2) with
    the bracket 2*3*4*5 = 120 in place of 720, and one factor 1/alpha = c_max/u left over --
    which is where (T2) gets its c_max u^5 rather than u^6.

WHAT IT IS FOR
--------------
The ceiling of the paper is the exact optimum of a problem built on the truncated surface.
(T1)-(T2) are what would let that become a statement about the FULL polylogarithm: (T2)
bounds how far the true stationary point can sit from the closed-form one, so the certified
1/R_5 <= 10.03 TeV becomes 1/R_5 <= 10.03 TeV + (an explicit, computable margin).

Run:  python remainder_theorem.py
"""
import sys

from mpmath import mp, polylog, zeta, mpf, mpc, log, pi, exp

mp.dps = 40
P = lambda *a: print(*a, flush=True)

H4 = mpf(25) / 12
Z3 = zeta(3)
LN2 = log(2)


def moments(terms):
    A2 = B2 = A4 = B4 = A4L = mpf(0)
    for (c, s, m) in terms:
        c = mpf(c)
        if s > 0:
            A2 += m * c ** 2
            A4 += m * c ** 4
            A4L += m * c ** 4 * log(c)
        else:
            B2 += m * c ** 2
            B4 += m * c ** 4
    return (A2 - mpf(3) / 4 * B2, A4, A4 * H4 - A4L - LN2 * B4)


def F(terms, x):
    return sum(m * (polylog(5, s * exp(mpc(0, 1) * c * x))).real for (c, s, m) in terms)


def head(terms, x):
    """the three rungs the paper keeps, plus the constant."""
    D, A4, G = moments(terms)
    return F(terms, mpf('1e-25')) - Z3 * D / 2 * x ** 2 + x ** 4 / 24 * (G - A4 * log(x))


def bound_R(terms, alpha):
    N = sum(m for (_, _, m) in terms)
    cmax = max(c for (c, _, _) in terms)
    u = mpf(cmax) * alpha
    if u >= 1:
        return None, u
    return pi ** 6 * N * u ** 6 / (2160 * (1 - u ** 2)), u


def bound_R_sharp(terms, alpha):
    """(T1s) -- the same bound with the CHARGE SPECTRUM instead of (N, c_max).

    Step 3 of the derivation throws the whole multiplicity onto the largest charge:
    |A_2k + (2^{2k-4}-1) B_2k| <= c_max^{2k} 2^{2k-4} N.  But the moments are traces over the
    weights of T in R, so for k >= 2 one can keep the spectrum and peel off only the excess:

        A_2k <= c_max^{2k-4} A_4 ,   B_2k <= c_max^{2k-4} B_4

    which gives  |...| <= c_max^{2k-4} 2^{2k-4} (A_4 + B_4)  and, collecting as before,

        |R| <= pi^6 (A_4 + B_4) c_max^2 alpha^6 / ( 2160 (1 - u^2) ) .              (T1s)

    This is never worse than (T1) -- the ratio is (A_4+B_4)/(N c_max^4) <= 1 -- and it is
    written in the fourth moment, which the integer program already carries as a variable.
    That is the practical content of "use the group": not a new inequality, the same one with
    the trace kept instead of bounded.
    """
    cmax = mpf(max(c for (c, _, _) in terms))
    u = cmax * alpha
    if u >= 1:
        return None, u
    S4 = sum(m * mpf(c) ** 4 for (c, _, m) in terms)      # A_4 + B_4
    return pi ** 6 * S4 * cmax ** 2 * alpha ** 6 / (2160 * (1 - u ** 2)), u


def bound_dR(terms, alpha):
    N = sum(m for (_, _, m) in terms)
    cmax = max(c for (c, _, _) in terms)
    u = mpf(cmax) * alpha
    if u >= 1:
        return None, u
    # NO es pi^7 N u^6/(360(1-u^2)).  Ese fue el primer intento y la comprobacion lo mato:
    # al pasar de x a alpha hay un 1/alpha, y yo lo sustitui por pi.  1/alpha = c_max/u,
    # asi que el factor correcto es c_max/u, no pi.  Con el mal, la cota quedaba POR DEBAJO
    # del resto real en cinco de los veinte casos -- es decir, no era una cota.
    return pi ** 6 * N * mpf(cmax) * u ** 5 / (360 * (1 - u ** 2)), u


CASES = [
    [(1, 1, 1), (2, -1, 1)],
    [(1, 1, 2), (2, 1, 1), (3, 1, 1)],
    [(1, -1, 3), (2, 1, 2)],
    [(1, 1, 1), (1, -1, 2), (2, 1, 1), (3, 1, 1)],
    [(2, -1, 2), (3, 1, 1)],
    [(1, 1, 4), (2, -1, 1), (3, 1, 2)],
    [(1, 1, 3), (1, -1, 3), (2, 1, 2), (2, -1, 2), (3, 1, 2)],
    [(1, -1, 6), (2, 1, 6), (3, 1, 3)],
]
ALPHAS = ['0.02', '0.05', '0.083', '0.12']


def main():
    P("=" * 104)
    P("(T1)   |R(alpha)|  <=  pi^6 N u^6 / (2160 (1 - u^2)) ,     u = c_max * alpha")
    P("=" * 104)
    P("%-30s %5s %5s %7s %12s %11s %7s %11s %7s"
      % ("content", "N", "cmax", "alpha", "|R| actual",
         "(T1)", "slack", "(T1s) sharp", "slack"))
    worst_ok, viol, worst_s, viol_s = None, 0, None, 0
    for terms in CASES:
        for a in ALPHAS:
            al = mpf(a)
            x = pi * al
            b, u = bound_R(terms, al)
            bs, _ = bound_R_sharp(terms, al)
            act = abs(F(terms, x) - head(terms, x))
            if b is None:
                continue
            ratio = float(b / act) if act > 0 else float('inf')
            rs = float(bs / act) if act > 0 else float('inf')
            viol += ratio < 1
            viol_s += rs < 1
            worst_ok = ratio if worst_ok is None else min(worst_ok, ratio)
            worst_s = rs if worst_s is None else min(worst_s, rs)
            P("%-30s %5d %5d %7s %12.3e %11.3e %7.1f %11.3e %7.1f"
              % (str(terms)[:30], sum(m for _, _, m in terms),
                 max(c for c, _, _ in terms), a, float(act),
                 float(b), ratio, float(bs), rs))
    P("")
    P("  violations, (T1) : %d      (T1s) : %d      <-- both must be 0" % (viol, viol_s))
    P("  tightest slack, (T1) : %.1f x     (T1s) : %.1f x  -- the spectrum is worth that factor"
      % (worst_ok, worst_s))

    P("")
    P("=" * 104)
    P("(T2)   |dR/dalpha|  <=  pi^6 N c_max u^5 / (360 (1 - u^2))")
    P("=" * 104)
    P("%-34s %8s %13s %13s %9s" % ("content", "alpha", "|R'| actual", "bound", "slack"))
    viol2 = 0
    h = mpf('1e-12')
    for terms in CASES[:5]:
        for a in ALPHAS:
            al = mpf(a)
            b, u = bound_dR(terms, al)
            if b is None:
                continue
            d = ((F(terms, pi * (al + h)) - head(terms, pi * (al + h)))
                 - (F(terms, pi * (al - h)) - head(terms, pi * (al - h)))) / (2 * h)
            act = abs(d)
            ratio = float(b / act) if act > 0 else float('inf')
            viol2 += ratio < 1
            P("%-34s %8s %13.3e %13.3e %9.1f"
              % (str(terms)[:34], a, float(act), float(b), ratio))
    P("")
    P("  violations : %d   <-- must be 0" % viol2)

    P("")
    P("=" * 104)
    P("CONTROLS -- the bound must refuse, and the remainder must really blow up, at u >= 1")
    P("=" * 104)
    t = [(3, -1, 1)]
    for a in ['0.30', '0.33', '0.34', '0.40']:
        al = mpf(a)
        b, u = bound_R(t, al)
        act = abs(F(t, pi * al) - head(t, pi * al))
        P("   c_max=3, alpha=%s  ->  u = %.3f   %s   |R| actual = %.3e"
          % (a, float(u), "bound = %.3e" % float(b) if b is not None
             else "OUT OF HYPOTHESIS, no bound claimed", float(act)))
    P("")
    P("   and the tail itself, at u just under 1 -- the bound is finite but useless, as it must be:")
    for a in ['0.30', '0.32', '0.333']:
        al = mpf(a)
        b, u = bound_R(t, al)
        P("      u = %.4f   bound = %.4e" % (float(u), float(b)))

    P("")
    P("=" * 104)
    P("READ THE SLACK BEFORE USING THIS.  The bound is uniform in the content -- it knows only")
    P("N and c_max -- so it pays for that with a factor of order 10^2-10^3 on any single case.")
    P("It is a BOUND, not an estimate: the estimate is (A_6 + 3 B_6) x^6 / 8640, which is sharp")
    P("and is not rigorous.  The two are different objects and the paper must not mix them.")
    P("=" * 104)
    return 1 if (viol or viol2) else 0


if __name__ == "__main__":
    sys.exit(main())
