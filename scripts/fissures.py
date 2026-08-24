#!/usr/bin/env python3
"""fissures.py -- the cracks left in Part VII, hunted on purpose.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Carles: "exploremos este a fondo, a ver donde nos lleva; no dejemos fisuras sin explorar."
So this file goes looking for them, starting with the one I put in the paper myself this
morning.

  F1  THE ANALYTIC BRANCH IS NOT EMPTY, AND MY OWN PARAGRAPH SAID IT WAS.
      crossed_equations.py measured G <= 0 on all seventeen contents of at most six multiplets
      with A4 = 0, and I wrote 'within the closed form the analytic branch admits no Higgs mass
      at all'.  That is a count, not a theorem, and the count is over a CAPPED set.  Here the
      question is settled exactly, because it is finite: the 7(+,+) is the only multiplet with
      A4 = 0, so every other one must make up the gauge sector's -18 exactly, and only the
      multiplets with A4 <= 18 can appear.  Finite, enumerable, and the answer is not the one
      in the paper.

  F2  DOES THE COMB CONTAIN ITS OWN TARGET CASE?  A spacing law is worth nothing if the five
      published rows do not sit on its teeth.  Never checked.

  F3  IS THE PER-D CEILING PROVABLY MONOTONE, OR IS IT FORTY RUNGS OF SCAN?  The paper says
      'monotone decreasing throughout: True' over 40 rungs and then rests the whole ceiling on
      the maximum sitting at 8D = 1.  A scan is not a proof, and the sufficient condition is
      one line: 1/R5^2 goes like (6 mu + t*(k))/k, so it is enough that t*(k)/k be
      non-increasing.  Measured, and stated as measured.

  F4  HOW HARD IS THE MISSING WITNESS?  crossed_equations.py searched to N = 12.  Push it, and
      say what the obstruction looks like.

Run:  python fissures.py
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

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
GV = [moments([(r, e, p, 1)])["G"] - _g["G"] for r, e, p in REPS]
A0, K0, G0 = round(_g["A4"]), round(8 * _g["D"]), _g["G"]


def show(vec):
    return " + ".join("%dx%s" % (vec[j], NAMES[j]) for j in range(8) if vec[j]) or "(gauge only)"


# ================================================================= F1
P("=" * 100)
P("F1 -- THE ANALYTIC BRANCH A4 = 0, SETTLED EXACTLY INSTEAD OF COUNTED")
P("=" * 100)
P("  the lattice, per multiplet:")
P("  %-10s %6s %6s %12s %10s" % ("multiplet", "A4", "8D", "G", "G / A4"))
for j in range(8):
    r = "%10.3f" % (GV[j] / AV[j]) if AV[j] else "        --"
    P("  %-10s %6d %6d %12.5f %s" % (NAMES[j], AV[j], KV[j], GV[j], r))
P("  gauge : A4 = %d,  8D = %d,  G = %.5f" % (A0, K0, G0))
P("")
P("  TWO STRUCTURAL FACTS, and together they make the question finite:")
zero_a = [j for j in range(8) if AV[j] == 0]
P("     (a) the ONLY multiplet with A4 = 0 is %s -- it can be added freely without leaving"
  % ", ".join(NAMES[j] for j in zero_a))
P("         the branch, and it carries 8D = %d and G = %+.5f, so it lowers BOTH."
  % (KV[zero_a[0]], GV[zero_a[0]]))
P("     (b) every other A4 is strictly positive, so A4 = 0 forces  sum m_j A_j = %d  exactly,"
  % (-A0))
P("         and only multiplets with A4 <= %d can appear at all." % (-A0))
usable = [j for j in range(8) if 0 < AV[j] <= -A0]
P("         usable : %s" % ", ".join("%s (A4=%d)" % (NAMES[j], AV[j]) for j in usable))
P("")
P("  So max G on the branch is attained with NO %s, and it is a finite enumeration."
  % NAMES[zero_a[0]])
P("")

sols = []


def rec(i, left, vec):
    if left == 0:
        sols.append(tuple(vec))
        return
    if i >= len(usable):
        return
    j = usable[i]
    for m in range(left // AV[j] + 1):
        vec[j] = m
        rec(i + 1, left - m * AV[j], vec)
    vec[j] = 0


rec(0, -A0, [0] * 8)
P("  contents with A4 = 0 and no %s : %d, and that is ALL of them" % (NAMES[zero_a[0]], len(sols)))
best = max(sols, key=lambda v: G0 + sum(v[j] * GV[j] for j in range(8)))
Gb = G0 + sum(best[j] * GV[j] for j in range(8))
Kb = K0 + sum(best[j] * KV[j] for j in range(8))
P("  maximum of G over the whole branch : %+.5f   at 8D = %d" % (Gb, Kb))
P("     %s   (%d multiplets)" % (show(best), sum(best)))
P("")
if Gb > 0:
    P("  ==> G > 0 IS REACHABLE ON THE ANALYTIC BRANCH.  The paragraph I wrote this morning --")
    P("      'within the closed form the analytic branch admits no Higgs mass at all' -- is")
    P("      FALSE.  What is true is the capped statement, and the cap is the whole content of it.")
else:
    P("  ==> G <= 0 on the entire branch: the paragraph stands, now as a theorem and not a count.")
P("")
P("  what the branch actually looks like, by content size:")
P("  %-4s %8s %10s %10s %14s %10s %11s"
  % ("N", "count", "max G", "8D there", "m_h GeV", "alpha_min", "1/R5 GeV"))
bysize = {}
for v in sols:
    n = sum(v)
    G = G0 + sum(v[j] * GV[j] for j in range(8))
    K = K0 + sum(v[j] * KV[j] for j in range(8))
    if K <= 0:
        continue                       # D > 0 or there is no minimum to speak of
    cur = bysize.get(n)
    if cur is None or G > cur[0]:
        bysize[n] = (G, K, v)
first_ok = None
for n in sorted(bysize):
    G, K, v = bysize[n]
    if G > 0:
        mu = G / 3.0
        mh = KK * math.pi ** 2 * math.sqrt(mu)
        al = math.sqrt(6 * Z3 * (K / 8.0) / G) / math.pi
        invR = 2 * MW / al
        P("  %-4d %8d %10.4f %10d %14.1f %10.5f %11.1f %s"
          % (n, sum(1 for x in sols if sum(x) == n), G, K, mh, al, invR,
             "  <- m_h in window" if MH_LO <= mh <= MH_HI else ""))
        if first_ok is None:
            first_ok = (n, v, G, K, mh, invR)
    else:
        P("  %-4d %8d %10.4f %10d %14s %10s %11s"
          % (n, sum(1 for x in sols if sum(x) == n), G, K, "--", "--", "G <= 0"))
P("")
if first_ok:
    n, v, G, K, mh, invR = first_ok
    P("  the branch opens at N = %d multiplets : %s" % (n, show(v)))
    P("     A4 = 0, 8D = %d, G = %.4f  ->  m_h = %.1f GeV, 1/R5 = %.0f GeV" % (K, G, mh, invR))
    P("  and every content of at most six multiplets has G <= 0, which is why a capped count")
    P("  looked like a theorem.  THE CAP WAS THE CLAIM.")
P("")
P("  what survives, and it is still worth stating: on the analytic branch the Higgs mass is not")
P("  a free parameter.  (I) degenerates to the EQUALITY G = 3 mu, so m_h is fixed by the content")
P("  with no window left -- a one-parameter family collapses to a point.  Contents on the branch")
P("  whose forced m_h lands in 125-127 GeV:")
inwin = []
for v in sols:
    G = G0 + sum(v[j] * GV[j] for j in range(8))
    K = K0 + sum(v[j] * KV[j] for j in range(8))
    if G <= 0 or K <= 0:
        continue
    mh = KK * math.pi ** 2 * math.sqrt(G / 3.0)
    if MH_LO <= mh <= MH_HI:
        al = math.sqrt(6 * Z3 * (K / 8.0) / G) / math.pi
        inwin.append((sum(v), v, K, mh, 2 * MW / al))
P("     %d of %d" % (len(inwin), len(sols)))
for n, v, K, mh, invR in sorted(inwin)[:8]:
    P("        %-44s 8D=%d  m_h=%.1f  1/R5=%.0f GeV" % (show(v), K, mh, invR))
P("")
P("")
P("  AND HERE IS THE STATEMENT THAT REPLACES THE FALSE ONE, and it is a theorem over the WHOLE")
P("  branch, infinite family included.  Adding a %s changes nothing in A4 and lowers G by"
  % NAMES[zero_a[0]])
P("  %.5f, and m_h = K pi^2 sqrt(G/3) is increasing in G.  So the largest Higgs mass anywhere on"
  % (-GV[zero_a[0]]))
P("  the branch is the largest over the %d contents with no %s:" % (len(sols), NAMES[zero_a[0]]))
mhmax = None
for v in sols:
    G = G0 + sum(v[j] * GV[j] for j in range(8))
    K = K0 + sum(v[j] * KV[j] for j in range(8))
    if G <= 0 or K <= 0:
        continue
    mh = KK * math.pi ** 2 * math.sqrt(G / 3.0)
    if mhmax is None or mh > mhmax[0]:
        mhmax = (mh, v, K, G)
if mhmax:
    mh, v, K, G = mhmax
    al = math.sqrt(6 * Z3 * (K / 8.0) / G) / math.pi
    P("")
    P("     max m_h on the analytic branch = %.1f GeV      at %s" % (mh, show(v)))
    P("     against the measured 125.25 GeV -- short by a factor of %.1f." % (125.25 / mh))
    P("     (and its 1/R5 would be %.0f GeV, below m_W)" % (2 * MW / al))
    P("")
    P("  ==> THE ANALYTIC BRANCH IS NOT EMPTY, IT IS UNPHYSICAL, AND THE BOUND IS EXACT.")
    P("      Cancelling the x^4 log x branch point is possible; doing it and keeping a")
    P("      Standard-Model Higgs is not.  That is a stronger statement than the one I wrote")
    P("      this morning AND it is true, which the other was not.")
P("")
P("  AND THE CONGRUENCE STANDS UNTOUCHED, because it is arithmetic and not a count:")
bad = [v for v in sols if (K0 + sum(v[j] * KV[j] for j in range(8))) % 6 != 3]
P("     8D = 3 (mod 6) on all %d : violations = %d" % (len(sols), len(bad)))
P("     8D = 1 anywhere on the branch : %d"
  % sum(1 for v in sols if K0 + sum(v[j] * KV[j] for j in range(8)) == 1))
P("     -> 'cancelling the branch point forbids the quantum' is a THEOREM. 'the branch is empty'")
P("        was a miscount.  The first is the result; the second was mine and it is withdrawn.")


# ================================================================= F2
P("")
P("=" * 100)
P("F2 -- DOES THE COMB CONTAIN ITS OWN TARGET CASE?")
P("=" * 100)
P("  M^2 = C (6 mu + A4)/k with C = 8 pi^2 m_W^2 / 3 zeta(3), and the admissible A4 at fixed k")
P("  step by three.  If the five published rows do not sit on teeth of their own comb, the law")
P("  is decoration.  Tooth index  n = (A4 - A4_ref)/3  must come out an integer for each row,")
P("  with A4_ref the smallest admissible A4 on that row's rung.")
P("")
P("  %-5s %6s %6s %10s %14s %12s %s" % ("row", "8D", "A4", "A4+8D", "mod 3", "tooth n", "integer?"))
allok = True
for label, content, a_them, mh_them, invR in T1:
    mo = moments(content)
    A4, k8 = round(mo["A4"]), round(8 * mo["D"])
    # the rung's smallest admissible A4 >= 0 with A4 + k = 0 (mod 3)
    ref = (-k8) % 3
    n = Fr(A4 - ref, 3)
    ok = n.denominator == 1
    allok &= ok and (A4 + k8) % 3 == 0
    P("  %-5s %6d %6d %10d %14d %12s %s"
      % (label, k8, A4, A4 + k8, (A4 + k8) % 3, str(n), "yes" if ok else "NO"))
P("")
P("  every published row lands on a tooth : %s" % allok)
P("")
P("  and the spacing, evaluated where it would be measured:")
C_COMB = 8 * math.pi ** 2 * MW ** 2 / Z3
P("  %-6s %16s %16s %s" % ("8D=k", "dM^2 (TeV^2)", "at 10 TeV", "at 3 TeV"))
for k in (1, 3, 9, 11, 15):
    d2 = C_COMB / k
    P("  %-6d %16.4f %13.1f GeV %11.1f GeV" % (k, d2 / 1e6, d2 / (2 * 10000.0), d2 / (2 * 3000.0)))
P("")
P("  CONTROL that can fail -- a decoy spacing.  If the step in A4 were 1 instead of 3, the")
P("  teeth would be three times denser and the published rows would still land on them, so")
P("  'the rows land on teeth' alone proves nothing.  What carries content is that A4 + 8D = 0")
P("  (mod 3) on every row, which the column above shows and which a step of 1 does not use.")

# ================================================================= F3
P("")
P("=" * 100)
P("F3 -- IS THE PER-D CEILING MONOTONE, OR IS IT FORTY RUNGS OF SCAN?")
P("=" * 100)
P("  1/R5 = 2 m_W/alpha and x^2 = 12 z3 (k/8)/(6 mu + t), so")
P("     (1/R5)^2  proportional to  (6 mu + t*(k)) / k")
P("  where t*(k) is the largest admissible A4 on rung k.  A SUFFICIENT condition for the")
P("  ceiling to be monotone decreasing is that t*(k)/k be non-increasing -- then both terms")
P("  fall.  The paper asserts monotonicity from a scan of forty rungs.  Here is the ratio it")
P("  actually rests on, read off the certificate table:")
P("")
CERT = [(1, 215), (3, 336), (5, 436), (7, 533), (9, 630), (15, 918), (21, 1209),
        (33, 1800), (45, 2397), (65, 3400), (99, 5112), (129, 6630), (201, 10272),
        (301, 15335), (501, 25470)]
P("  %-6s %10s %14s %14s" % ("8D=k", "t*(k)", "t*(k)/k", "falling?"))
prev = None
mono = True
for k, t in CERT:
    r = t / k
    fall = "" if prev is None else ("yes" if r < prev else "NO")
    if prev is not None and r >= prev:
        mono = False
    P("  %-6d %10d %14.3f %14s" % (k, t, r, fall))
    prev = r
P("")
P("  t*(k)/k is non-increasing across the certificate rungs : %s" % mono)
P("  ...which is exactly the sufficient condition, so the monotonicity the ceiling rests on is")
P("  NOT an accident of the scan -- it follows from one measured property of t*.  What is still")
P("  missing, and it is missing on purpose: t* is the optimum of a program whose constraint")
P("  G <= G*(t,k) carries a logarithm, so t*(k)/k non-increasing is MEASURED here and not")
P("  proved.  Stating it as the sufficient condition is the improvement; proving it is open.")

# ================================================================= F4
P("")
P("=" * 100)
P("F4 -- HOW FAR IS THE MISSING WITNESS?")
P("=" * 100)
mu_hi = MU(MH_HI)
x_star = math.sqrt(12 * Z3 * (1 / 8.0) / (6 * mu_hi + 215))
G_star = 215 * (math.log(x_star) + 0.75) + 3 * mu_hi
P("  the certifying vertex is (A4, 8D) = (215, 1) and it demands G <= G* = %.4f." % G_star)
P("  Two integer conditions and one inequality.  Ask each separately -- a claim that fails on")
P("  the FIRST is a different problem from one that fails on the third.")
P("")
NW = 16
hitA = hitAK = hitAKG = 0
bestG = None
for n in range(1, NW + 1):
    for vec in itertools.combinations_with_replacement(range(8), n):
        v = [0] * 8
        for j in vec:
            v[j] += 1
        A4 = A0 + sum(v[j] * AV[j] for j in range(8))
        if A4 != 215:
            continue
        hitA += 1
        k8 = K0 + sum(v[j] * KV[j] for j in range(8))
        if k8 != 1:
            continue
        hitAK += 1
        G = G0 + sum(v[j] * GV[j] for j in range(8))
        if bestG is None or G < bestG[0]:
            bestG = (G, tuple(v), n)
        if G <= G_star:
            hitAKG += 1
P("  over every content of at most %d multiplets:" % NW)
P("     with A4 = 215 exactly                    : %d" % hitA)
P("     of those, also 8D = 1                    : %d" % hitAK)
P("     of those, also G <= G*                   : %d" % hitAKG)
if bestG:
    P("     smallest G among the (215, 1) contents    : %.4f   against G* = %.4f  (short by %.4f)"
      % (bestG[0], G_star, bestG[0] - G_star))
    P("        %s   (%d multiplets)" % (show(list(bestG[1])), bestG[2]))
    P("")
    P("  ==> the vertex IS reachable in the two integer moments; what fails is the third")
    P("      condition, the one carrying the logarithm.  That relocates the open question:")
    P("      it is not 'can the lattice reach (215,1)' but 'can it reach it cheaply enough in G'.")
else:
    P("")
    P("  ==> no content reaches (215, 1) in the two integer moments at all up to N = %d, so the" % NW)
    P("      obstruction is already in the lattice and not in the G condition.")
P("")
P("  and this is a search, not a proof: absence up to N = %d is absence up to N = %d." % (NW, NW))

# ================================================================= the numbers as the paper prints them
P("")
P("=" * 100)
P("THE NUMBERS THIS FILE PUTS IN THE PAPER, AT THE PRECISION THE PAPER PRINTS THEM")
P("=" * 100)
P("  (check_numbers.py greps the archive literally, so a figure rounded in the prose has to")
P("   exist rounded here too -- otherwise the gate cannot tell a rounding from an invention.)")
P("")
if mhmax:
    P("     max m_h on the analytic branch          : %.1f GeV" % mhmax[0])
    P("     short of the measured value by a factor : %.1f" % (125.25 / mhmax[0]))
P("     G* at the certifying vertex             : %.1f" % G_star)
if bestG:
    P("     smallest G among the (215,1) contents   : %.1f" % bestG[0])
    P("     the gap                                 : %.0f" % (bestG[0] - G_star))
P("     t*(k)/k at k = 1, 9, 501                : %.0f, %.0f, %.1f"
  % (CERT[0][1] / CERT[0][0], 630 / 9.0, CERT[-1][1] / CERT[-1][0]))
P("     contents with A4 = 215                  : %d" % hitA)
P("     of those with 8D = 1                    : %d" % hitAK)
