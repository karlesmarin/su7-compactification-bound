#!/usr/bin/env python3
"""crossed_equations.py -- what falls out when the paper's sections are crossed against each other.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Carles read Part VII hunting for the apparently secondary sentences, the limit cases and the
consequences that only appear when equations from different sections are multiplied together.
This file measures the four of his findings that are measurable, and it is built so each one can
FAIL.  Two of them are new results, one is a wording defect in the paper, and one is a decoy test
that is meant to come out negative if the idea is wrong.

  1  THE LADDER PAST THE POLE.  p = 5-2k does not stop at p=1: it continues, and for p<0 every
     zeta is RATIONAL (Bernoulli).  So exactly ONE logarithm exists in the whole expansion and the
     transcendental content of the potential is {zeta(5), zeta(3), ln x, ln 2, ln 3} and nothing
     else.  Already derived in REMAINDER_BOUND.md section 1 step 2; reproduced here because
     section 2 needs the same machinery and because the paper's own text has not caught up.

  2  A4 = 0 IS AN ANALYTIC ISLAND, AND IT COSTS THE QUANTUM.  A4 is the coefficient of x^4 ln x.
     Theorem 2 says 8D = 2 A4 + 3 (mod 6), so A4 = 0 forces 8D = 3 (mod 6): 8D in {3, 9, 15, ...}
     and NEVER 1.  Cancelling the branch point forbids the cheapest curvature, which is exactly
     the one the ceiling of section 8 is attained on.  Maximal hierarchy and maximal analyticity
     are incompatible, and the congruence is what makes them so.

  3  THE KK MASS SITS ON AN ARITHMETIC COMB.  Crossing (II) with Theorem 2: at fixed 8D = k the
     admissible A4 move in steps of three, so M_KK^2 is quantised with spacing
     8 pi^2 m_W^2 / (zeta(3) k), independent of the content.

  4  THE c = 2 PARITY FLIP HAS A LOG-FREE SIGNATURE, AND IT IS THE ONLY CHARGE THAT DOES.
     Moving one unit of multiplicity at charge two from periodic to antiperiodic shifts
     (D, A4, G) by (-7, -16, -100/3) exactly -- the logarithms cancel, because ln 2 is BOTH the
     log of the charge and the universal log of the antiperiodic tower.  At c = 1 and c = 3 the
     shift carries ln 2 and ln(3/2) and is irrational.  So if Part VI's anchor residual were a
     misread c=2 parity, the discrepancy would be a RATIONAL multiple of that vector -- which is
     a test, run below, with two decoys.

  5  'ATTAINED' NEEDS A WITNESS.  The paper says the ceiling is 'attained at (A4, 8D) = (215, 1)'.
     The dual certifies the BOUND; 'attained' is a further claim and needs an integer content.
     Measured here.

Run:  python crossed_equations.py
"""
import itertools
import math
import pathlib
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_LO, MH_HI = 125.0, 127.0
LN2, LN3 = math.log(2), math.log(3)

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
GV = [moments([(r, e, p, 1)])["G"] - _g["G"] for r, e, p in REPS]
A0, K0, G0 = round(_g["A4"]), round(8 * _g["D"]), _g["G"]


def mom_of(vec):
    """(A4, 8D, G) of the content with multiplicity vector vec, gauge included."""
    return (A0 + sum(v * AV[j] for j, v in enumerate(vec)),
            K0 + sum(v * KV[j] for j, v in enumerate(vec)),
            G0 + sum(v * GV[j] for j, v in enumerate(vec)))


# ============================================================== 1
P("=" * 100)
P("1 -- THE LADDER PAST THE POLE: after p = 1 the transcendence stops")
P("=" * 100)
P("  Re Li_5(s e^{icx}) = SUM_k (-1)^k (cx)^{2k}/(2k)! L_k(s),  L_k(+1) = zeta(5-2k),")
P("  L_k(-1) = -eta(5-2k).  The paper's three rungs are k = 0, 1, 2 (p = 5, 3, 1).  For k >= 3")
P("  the argument p = 5-2k is a NEGATIVE ODD integer, where zeta is rational by Bernoulli.")
P("")
P("  %-4s %-5s %-16s %-22s %s" % ("k", "p", "zeta(p)", "eta(p)/zeta(p)", "bracket weight 2^{2k-4}-1"))
for k in range(0, 8):
    p = 5 - 2 * k
    if p == 1:
        P("  %-4d %-5d %-16s %-22s %s" % (k, p, "POLE", "-> the ln x", "-- the branch point"))
        continue
    zp = Fr(0)
    try:
        import mpmath
        z = mpmath.zeta(p)
        e = (1 - mpmath.mpf(2) ** (1 - p)) * z
        rat = e / z if z != 0 else mpmath.mpf("nan")
        zs = mpmath.nstr(z, 10)
        rs = mpmath.nstr(rat, 10)
    except ImportError:
        zs = rs = "?"
    w = 2 ** (2 * k - 4) - 1 if k >= 3 else None
    P("  %-4d %-5d %-16s %-22s %s" % (k, p, zs, rs, "" if w is None else str(w)))
P("")
P("  the bracket sequence for k = 3, 4, 5, 6 is  %s  =  2^{2k-4} - 1"
  % ", ".join(str(2 ** (2 * k - 4) - 1) for k in (3, 4, 5, 6)))
P("  so the rung is  zeta(5-2k) [ A_{2k} + (2^{2k-4}-1) B_{2k} ],  which is the paper's eq. (6)")
P("  at p = 5-2k -- the SAME ladder, continued past its own pole.")
P("")
P("  CONSEQUENCE, and it is the one that matters for the text: there is EXACTLY ONE logarithm in")
P("  the whole expansion, at order x^4, and everything beyond it is rational in the moments.  The")
P("  regular part is therefore a CONVERGENT series, not an asymptotic one, inside the radius set")
P("  by the nearest singularity of Li_5(s e^{i c pi alpha}) -- s e^{ic pi alpha} = 1, i.e.")
P("     periodic     alpha = 2n/c      nearest non-zero  2/c_max")
P("     antiperiodic alpha = (2n+1)/c  nearest           1/c_max")
P("  and with c_max = 3 that is alpha = 1/3.  Measured, with a control outside the radius:")
P("")


def term(c, s, alpha, k):
    """|the k-th term| of Re Li_5(s e^{i c pi alpha}).  Convergence is decided by whether the
    TERMS shrink -- comparing partial sums is not a test, it is a tautology."""
    import mpmath
    mpmath.mp.dps = 60
    x = mpmath.pi * alpha
    p = 5 - 2 * k
    z = mpmath.zeta(p)
    L = z if s > 0 else -(1 - mpmath.mpf(2) ** (1 - p)) * z
    # BOTH towers are expanded in the SAME variable c x -- the parity lives in the WEIGHT
    # (zeta for periodic, -eta for antiperiodic), not in the argument.  Mixing the eta weights
    # with the shifted argument c x + pi is expanding twice and it reads as a divergence.
    return abs((c * x) ** (2 * k) / mpmath.factorial(2 * k) * L)


try:
    import mpmath
    P("  the terms behave like (c x / 2pi)^{2k} for the periodic tower and (c x / pi)^{2k} for the")
    P("  antiperiodic one -- the eta weights carry an extra 2^{2k}, which halves the radius.  So")
    P("  the two radii are alpha < 2/c and alpha < 1/c, exactly the two singularity families above.")
    P("  The terms have to shrink, and past the radius they must not: that is the control.")
    P("")
    P("  %-33s %9s %13s %13s %13s %11s"
      % ("case", "c*alpha", "|term k=6|", "|term k=30|", "ratio at k=200", "verdict"))
    for lbl, c, s, al in (("c=3 periodic  alpha=0.12", 3, +1, 0.12),
                          ("c=3 antiper.  alpha=0.12", 3, -1, 0.12),
                          ("c=3 antiper.  alpha=0.25", 3, -1, 0.25),
                          ("c=3 antiper.  alpha=0.32", 3, -1, 0.32),
                          ("CONTROL antiper. alpha=0.40 > 1/3", 3, -1, 0.40),
                          ("CONTROL period.  alpha=0.40 < 2/3", 3, +1, 0.40),
                          ("CONTROL period.  alpha=0.70 > 2/3", 3, +1, 0.70),
                          ("CONTROL alpha=1.50", 3, -1, 1.50)):
        t6, t30 = term(c, s, al, 6), term(c, s, al, 30)
        # the verdict is the ASYMPTOTIC ratio, not the first few terms: at alpha just past the
        # radius the terms still fall for a long while before the geometric factor takes over,
        # so reading k = 6..30 alone would call a divergent series convergent.
        r = float(term(c, s, al, 201) / term(c, s, al, 200))
        verd = "converges" if r < 1 else "DIVERGES"
        P("  %-33s %9.3f %13.3e %13.3e %13.4f %11s"
          % (lbl, c * al, float(t6), float(t30), r, verd))
    P("")
    P("  the two controls past alpha = 1/3 are what make this a measurement and not a hope.")
except ImportError:
    P("  (mpmath missing -- skipped)")

# ============================================================== 2
P("")
P("=" * 100)
P("2 -- A4 = 0 IS AN ANALYTIC ISLAND, AND THE CONGRUENCE CHARGES IT THE QUANTUM")
P("=" * 100)
P("  Theorem 2:  8D = 2 A4 + 3  (mod 6).   At A4 = 0 that is  8D = 3 (mod 6),")
P("  so 8D in {3, 9, 15, 21, ...} and 8D = 1 is FORBIDDEN.")
P("")
NMAX = 6
found = []
for n in range(1, NMAX + 1):
    for vec in itertools.combinations_with_replacement(range(8), n):
        v = [0] * 8
        for j in vec:
            v[j] += 1
        A4, k8, G = mom_of(v)
        if A4 == 0:
            found.append((tuple(v), A4, k8, G))
seen, uniq = set(), []
for v, A4, k8, G in found:
    if v not in seen:
        seen.add(v)
        uniq.append((v, A4, k8, G))
P("  contents of at most %d multiplets with A4 = 0 exactly : %d   (the paper says sixteen)"
  % (NMAX, len(uniq)))
pos = [u for u in uniq if u[2] > 0]
P("  of those, with D > 0 : %d   (the paper says five)" % len(pos))
P("")
bad = [u for u in uniq if u[2] % 6 != 3]
P("  CONTROL -- every one of them must satisfy 8D = 3 (mod 6): violations = %d" % len(bad))
P("  CONTROL -- and none of them may have 8D = 1: %d have it" % sum(1 for u in uniq if u[2] == 1))
P("")
P("  and on this branch (I) degenerates -- ln x = (G - 3mu)/A4 - 3/4 has A4 in the denominator --")
P("  so realisability is not an inequality but the EQUALITY  G = 3 mu.  The Higgs mass is then")
P("  fixed by the content with no window left, and so is alpha_min:")
P("")
P("     mu = G/3 ,  x^2 = 12 z3 D / (6 mu) = 6 z3 D / G ,  m_h = K pi^2 sqrt(mu) ,  1/R5 = 2 m_W/alpha")
P("")
P("  %-46s %5s %12s %10s %10s %11s"
  % ("content, D > 0 only", "8D", "G", "m_h GeV", "alpha_min", "1/R5 GeV"))
best = None
for v, A4, k8, G in sorted(pos, key=lambda u: u[2]):
    nm = " + ".join("%dx%s" % (v[j], NAMES[j]) for j in range(8) if v[j])
    if G <= 0:
        P("  %-46s %5d %12.4f %10s %10s %11s"
          % (nm, k8, G, "--", "--", "G <= 0"))
        continue
    mu = G / 3.0
    mh = KK * math.pi ** 2 * math.sqrt(mu)
    x2 = 6 * Z3 * (k8 / 8.0) / G
    al = math.sqrt(x2) / math.pi
    invR = 2 * MW / al
    P("  %-46s %5d %12.4f %10.1f %10.5f %11.1f %s"
      % (nm, k8, G, mh, al, invR, "  <- m_h in window" if MH_LO <= mh <= MH_HI else ""))
    if best is None or invR > best[0]:
        best = (invR, k8, mh, nm)
P("")
P("  AND THE BRANCH IS EMPTY, for a reason that is one line of algebra:")
P("     A4 = 0  =>  G = A4 H4 - A4L - ln2 B4 = -(A4L + ln2 B4)")
P("  so G > 0 needs A4L + ln2 B4 < 0, i.e. the WHOLE fourth-moment log budget negative.  Measured")
P("  over the %d contents with A4 = 0 : %d have G > 0."
  % (len(uniq), sum(1 for u in uniq if u[3] > 0)))
if best is None:
    P("")
    P("  => the closed form admits NO Higgs mass on the analytic branch at all.  It is not that")
    P("     cancelling the branch point costs hierarchy -- within the closed form it costs the")
    P("     stationary point.  The congruence 8D = 3 (mod 6) is then a second, independent")
    P("     obstruction sitting behind the first.")
else:
    P("")
    P("  best on the analytic branch : %.0f GeV = %.2f TeV   at 8D = %d"
      % (best[0], best[0] / 1000, best[1]))
    P("  against 10034 GeV = 10.03 TeV on the full lattice, which sits at 8D = 1.")
P("")
P("  => CANCELLING THE NON-ANALYTICITY AND MAXIMISING THE HIERARCHY ARE INCOMPATIBLE, and the")
P("     obstruction is a congruence, not a dynamical effect.  The smallest 8D actually realised")
P("     with A4 = 0 and D > 0 is %s, not 3." % (min(u[2] for u in pos) if pos else "--"))

# ============================================================== 3
P("")
P("=" * 100)
P("3 -- THE ARITHMETIC COMB IN M_KK^2")
P("=" * 100)
P("  (II) with D = k/8 gives   M^2 = (8 pi^2 m_W^2 / 3 z3) (6 mu + A4) / k ,   M = 1/R5 = 2 m_W/alpha")
P("  and Theorem 2 fixes A4 modulo 3 once k is fixed, so the admissible A4 step by 3:")
P("")
P("     Delta M^2 = 8 pi^2 m_W^2 / (z3 k)      -- independent of the content and of m_h")
P("")
C_COMB = 8 * math.pi ** 2 * MW ** 2 / Z3
P("  The exactly-spaced quantity is M^2, NOT M.  A spacing in mass is dM = dM^2 / 2M and is")
P("  therefore meaningless without saying at which M -- and each rung has its OWN ceiling, so")
P("  quoting them all at 10 TeV is quoting most of them where they cannot reach.")
P("")
P("  %-5s %14s %14s %12s %14s" % ("8D=k", "dM^2 (GeV^2)", "dM^2 (TeV^2)", "its ceiling", "dM there"))
CEIL = {1: 10013.0, 3: 6244.0, 5: 5121.0, 7: 4555.0, 11: 3967.0}
for k in (1, 3, 5, 7, 11):
    d2 = C_COMB / k
    M = CEIL[k]
    P("  %-5d %14.0f %14.4f %10.0f GeV %11.1f GeV" % (k, d2, d2 / 1e6, M, d2 / (2 * M)))
P("")
P("  CONTROL -- rebuild the certificate row from the comb.  At k = 1 the ceiling row is")
P("  A4 = 215, m_h = 127.0, 1/R5 = 10034 GeV.  From M^2 = (8pi^2 mW^2/3z3)(6mu+A4)/k :")
mu127 = MU(127.0)
M2 = (8 * math.pi ** 2 * MW ** 2 / (3 * Z3)) * (6 * mu127 + 215) / 1
P("     rebuilt 1/R5 = %.1f GeV   against 10034 GeV in outputs/ceiling_ilp.txt   (%s)"
  % (math.sqrt(M2), "ok" if abs(math.sqrt(M2) - 10034) < 5 else "MISMATCH"))
P("")
P("  the comb is NECESSARY, not sufficient: every real content lands on a tooth, not every tooth")
P("  holds a content.  And its ABSOLUTE position carries the anchor band and the g4 choice; what")
P("  does not is the SPACING, which is arithmetic.")

# ============================================================== 4
P("")
P("=" * 100)
P("4 -- THE c = 2 PARITY FLIP IS THE ONLY LOG-FREE ONE")
P("=" * 100)
P("  move one unit of multiplicity at charge c from the periodic to the antiperiodic tower:")
P("")
P("  %-4s %14s %14s %28s %s" % ("c", "Delta D", "Delta A4", "Delta G", "log-free?"))
H4v = float(H4)
for c in (1, 2, 3):
    dD = -c ** 2 - 0.75 * c ** 2
    dA4 = -c ** 4
    dG = -(c ** 4 * H4v - c ** 4 * math.log(c)) - LN2 * c ** 4
    # symbolic: dG = -c^4 H4 + c^4 ln c - c^4 ln 2 = -c^4 (H4 - ln(c/2))
    logfree = abs(math.log(c / 2.0)) < 1e-15
    P("  %-4d %14.4f %14d %28.6f %s"
      % (c, dD, dA4, dG, "YES  = -%s" % Fr(c ** 4 * 25, 12) if logfree else
         "no   carries ln(%d/2)" % c))
P("")
P("  the exact c = 2 signature is  (Delta D, Delta A4, Delta G) = (-7, -16, -100/3) per unit,")
P("  and it is log-free because ln c = ln 2 IS the universal log of the antiperiodic tower.")
P("  Equivalently  J = G - H4 A4 = -A4L - ln2 B4  is blind to the parity of a c = 2 state:")
for s in (+1, -1):
    A4c = 16.0 if s > 0 else 0.0
    A4L = 16.0 * LN2 if s > 0 else 0.0
    B4 = 0.0 if s > 0 else 16.0
    G = A4c * H4v - A4L - LN2 * B4
    P("     c=2 %-13s  A4 = %5.1f   G = %10.6f   J = G - H4*A4 = %12.6f"
      % ("periodic" if s > 0 else "antiperiodic", A4c, G, G - H4v * A4c))

P("")
P("  THE TEST.  If Part VI's anchor residual were a misread c = 2 parity, then for each published")
P("  row there is a t with  (D, A4, G)_theirs = (D, A4, G)_ours + t (-7, -16, -100/3).  Solve t")
P("  from (II) using THEIR alpha and THEIR m_h, then check (I), which t did not see.  Two")
P("  conditions, one unknown -- so it can fail, and two decoys are run beside it.")
P("")
VECS = [("c=2 flip  (-7, -16, -100/3)", -7.0, -16.0, -100.0 / 3.0),
        ("DECOY c=1 flip", -1.75, -1.0, -(25.0 / 12.0) - LN2),
        ("DECOY c=3 flip", -15.75, -81.0, -81.0 * (25.0 / 12.0) + 81.0 * (LN3 - LN2))]
P("  %-6s %-28s %10s %14s %14s" % ("row", "shift vector", "t", "(I) residual", "t integer?"))
for label, content, a_them, mh_them, invR in T1:
    mo = moments(content)
    Do, A4o, Go = mo["D"], mo["A4"], mo["G"]
    x = math.pi * a_them
    mu = MU(mh_them)
    for vname, vD, vA, vG in VECS:
        den = 12 * Z3 * vD - x * x * vA
        if abs(den) < 1e-12:
            continue
        t = (x * x * (6 * mu + A4o) - 12 * Z3 * Do) / den
        A4t, Gt = A4o + t * vA, Go + t * vG
        res = float("nan") if abs(A4t) < 1e-9 else math.log(x) - ((Gt - 3 * mu) / A4t - 0.75)
        near = abs(t - round(t)) < 0.02
        P("  %-6s %-28s %10.4f %14.4e %14s"
          % (label, vname, t, res, "YES t=%d" % round(t) if near else "no"))
    P("")
P("  read it as a falsification: a hit needs BOTH an integer t AND an (I) residual at the level")
P("  of the closed form's own accuracy.  Anything else is the vector not being the explanation.")

# ============================================================== 5
P("")
P("=" * 100)
P("5 -- 'ATTAINED AT (A4, 8D) = (215, 1)' -- IS THERE AN INTEGER WITNESS?")
P("=" * 100)
P("  the two-variable dual certifies the BOUND for arbitrary content.  'Attained' is a further")
P("  claim: it needs a multiplicity vector with A4 = 215, 8D = 1 and G <= G*(215, 1).")
P("")
mu_hi = MU(MH_HI)
x_star = math.sqrt(12 * Z3 * (1 / 8.0) / (6 * mu_hi + 215))
G_star = 215 * (math.log(x_star) + 0.75) + 3 * mu_hi
P("     G*(215, 1) = %.4f      x* = %.6f      1/R5 = %.1f GeV"
  % (G_star, x_star, 2 * MW / (x_star / math.pi)))
P("")
NW = 12
hits, bestk1 = [], None
for n in range(1, NW + 1):
    for vec in itertools.combinations_with_replacement(range(8), n):
        v = [0] * 8
        for j in vec:
            v[j] += 1
        A4, k8, G = mom_of(v)
        if k8 != 1:
            continue
        if A4 == 215 and G <= G_star:
            hits.append((tuple(v), A4, G))
        if bestk1 is None or A4 > bestk1[1]:
            mu_ = mu_hi
            x_ = math.sqrt(12 * Z3 * (1 / 8.0) / (6 * mu_ + A4))
            Gs_ = A4 * (math.log(x_) + 0.75) + 3 * mu_
            if G <= Gs_:
                bestk1 = (tuple(v), A4, G, 2 * MW / (x_ / math.pi))
P("  search over every content of at most %d multiplets at 8D = 1:" % NW)
P("     witnesses with A4 = 215 and G <= G* : %d" % len(hits))
if bestk1:
    nm = " + ".join("%dx%s" % (bestk1[0][j], NAMES[j]) for j in range(8) if bestk1[0][j])
    P("     largest A4 actually realised at 8D = 1 with G <= G* : A4 = %d  ->  1/R5 = %.0f GeV"
      % (bestk1[1], bestk1[3]))
    P("     that is %.2f TeV against the %.2f TeV the dual certifies -- the gap between a"
      % (bestk1[3] / 1000.0, 2 * MW / (x_star / math.pi) / 1000.0))
    P("     relaxation and its lattice.")
    P("        %s" % nm)
P("")
if not hits:
    P("  => NO INTEGER WITNESS at (215, 1) up to N = %d.  The dual bound stands -- it is a bound" % NW)
    P("     for arbitrary content and nothing here touches it -- but the word 'attained' claims")
    P("     more than the certificate delivers.  The paper should say the RELAXATION attains it")
    P("     at that vertex, and state the best realised value beside it.  Three places carry the")
    P("     word: the abstract, the anchor-trend paragraph, and the boxed result of the ceiling.")
else:
    P("  => witness found; 'attained' stands as written.")
P("")
P("  NOTE, and it is what keeps this honest: absence up to N = %d is not absence.  What is" % NW)
P("  measured is that the enumeration does not reach the vertex, not that nothing does.")
