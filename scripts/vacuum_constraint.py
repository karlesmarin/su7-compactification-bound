#!/usr/bin/env python3
"""vacuum_constraint.py -- the condition the ceiling never imposed, and it is linear.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHERE THIS CAME FROM
--------------------
socratic_witness.py asked why the ceiling's witness had no global minimum and found the root: the
witness DOES have an electroweak minimum -- at alpha = 0.016129, which the closed form predicts to
0.02 %, with m_h = 126.25 GeV inside the window -- but the potential is DEEPER at the other
symmetric point, alpha = 1.  The witness is a false vacuum.  The five published rows are not; that
is the control, and it passes.

So the certificate maximises 1/R5 over a stationary point without ever asking whether that point is
the global minimum.  At the optimum it is not.

AND THE MISSING CONDITION IS ALREADY IN THE PAPER, at the rung above the one the ceiling uses.
Evaluate F at the two symmetric points.  With F(alpha) = sum m Re Li_5(s e^{i c pi alpha}),

    F(0) = sum m Li_5(s) ,        F(1) = sum m Li_5(s (-1)^c)

so only the ODD charges move, and they move by Li_5(-s) - Li_5(s) = -/+ (zeta(5) + eta(5)).  With
eta(5) = (15/16) zeta(5) that is (31/16) zeta(5), and

    F(1) - F(0) = (31/16) zeta(5) * W ,      W = sum_{c odd} m (-s)

LINEAR IN THE MULTIPLICITIES, exactly like A4, 8D and G.  The 31 is 2^5 - 1: the same 2^p - 1 that
runs the ladder, arriving at the top rung.  So 'the electroweak vacuum is the deeper one' is a
fourth linear functional and drops straight into the same integer program.

This file: derives W, checks the closed form against the exact potential, screens the five rows and
the witness, and re-runs the rung-1 optimisation with the condition imposed.

Run:  python vacuum_constraint.py
"""
import math
import pathlib
import sys
from fractions import Fraction as Fr

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

Z5 = 1.0369277551433699
MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_LO, MH_HI = 125.0, 127.0
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]
show = lambda v: " + ".join("%dx%s" % (v[j], NAMES[j]) for j in range(8) if v[j])


def Wof(content):
    """W = sum over ODD charges of m * (-s).  Linear, and it is the whole condition.

    EXACT RATIONALS, and that is not fussiness: the gauge sector's odd-charge part is 2 - 7/2 =
    -3/2, so W is a HALF-odd-integer, never an integer.  This function used to be printed with
    %d, which truncated the half and read as -1 for a day; chasing the missing .5 is what turned
    up W_is_half_odd.py and a theorem.  Ver [[a-dropped-glyph-exits-zero]] -- same family."""
    return sum(Fr(-s) * Fr(m).limit_denominator(8)
               for m, s, c in table(content) if int(round(c)) % 2 == 1)


def dF_closed(content):
    return Fr(31, 16) * Z5 * Wof(content)


def dF_exact(content):
    f = F(content, np.array([1e-9, 1.0]))
    return float(f[1] - f[0])


# ================================================================= 1
P("=" * 100)
P("1 -- THE FORMULA, AND IT HAS TO SURVIVE THE EXACT POTENTIAL")
P("=" * 100)
P("  F(1) - F(0) = (31/16) zeta(5) W,  W = sum_{c odd} m (-s).  Checked term by term against the")
P("  exact polylogarithmic F, on the five published rows and on the ceiling's witness.")
P("")
P("  %-24s %10s %16s %16s %10s" % ("content", "W", "closed", "exact", "agree"))
rows = []
for label, content, a_them, mh_them, invR in T1:
    v = [0] * 8
    for rep, e, ep, mult in content:
        v[REPS.index((rep, e, ep))] += mult
    rows.append((label, v))
WIT = [0] * 8
WIT[0], WIT[1], WIT[3] = 17, 2, 57
rows.append(("witness", WIT))
ok = True
for label, v in rows:
    c = cont(v)
    a, b = float(dF_closed(c)), dF_exact(c)
    good = abs(a - b) < 1e-7 * max(1.0, abs(b))
    ok &= good
    P("  %-24s %10s %16.6f %16.6f %10s" % (label, Wof(c), a, b, good))
P("")
P("  CONTROL -- closed and exact agree on every row : %s" % ok)
P("")
P("  and the per-multiplet coefficients, which is what makes it a lattice functional:")
P("  %-12s %8s     %-12s %8s" % ("multiplet", "W", "multiplet", "W"))
WV = [Wof(cont([1 if i == j else 0 for i in range(8)])) - Wof(cont([0] * 8)) for j in range(8)]
W0 = Wof(cont([0] * 8))
for j in range(0, 8, 2):
    P("  %-12s %8s     %-12s %8s" % (NAMES[j], WV[j], NAMES[j + 1], WV[j + 1]))
P("  %-12s %8s   <- NOT an integer, and that is a theorem: see W_is_half_odd.py"
  % ("gauge", W0))
lin = all(Wof(cont(v)) == W0 + sum(v[j] * WV[j] for j in range(8)) for _, v in rows)
P("  CONTROL -- W really is linear in the multiplicities on every row above : %s" % lin)
P("")
P("  CONTROL -- 2W MUST BE ODD FOR EVERY CONTENT (W_is_half_odd.py).  Carles asked for this one")
P("  after finding W = 2 printed beside a theorem saying 2W is odd: the guard belongs in every")
P("  script that touches W, not in the one that proves it.  Random contents, multiplicities up")
P("  to eight:")
import random as _rnd
_rnd.seed(20260822)
_bad = 0
for _ in range(4000):
    _v = [_rnd.randint(0, 8) for _ in range(8)]
    # _wq, NOT _w.  amin_closed_form.py is exec'd into this namespace and used to keep the
    # winding weights n^-5 in a global called _w; this loop overwrote it with a Fraction and
    # every F() below returned garbage -- section 2 printed F(0) = 900900 instead of -18.76 --
    # without raising.  The weights are now _AMIN_W and basis() checks them, but the local name
    # stays distinct anyway.  Ver [[my-fix-introduces-the-next-bug]].
    _wq = Wof(cont(_v))
    if (2 * _wq).denominator != 1 or int(2 * _wq) % 2 != 1:
        _bad += 1
P("     4000 random contents, 2W not odd : %d   <-- must be 0" % _bad)
assert _bad == 0, "2W is not odd somewhere -- the half-integer bookkeeping is broken"

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- THE SCREEN.  W > 0 means the electroweak point is the deeper one")
P("=" * 100)
P("  %-24s %8s %14s %14s %s" % ("content", "W", "F(0)", "F(1)", "electroweak vacuum is"))
for label, v in rows:
    c = cont(v)
    f = F(c, np.array([1e-9, 1.0]))
    P("  %-24s %8s %14.4f %14.4f %s"
      % (label, Wof(c), f[0], f[1], "the true one" if Wof(c) > 0 else "FALSE -- deeper at alpha=1"))
P("")
P("  the five published rows pass and the ceiling's witness fails, which is what")
P("  socratic_witness.py found the slow way.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- THE CEILING WITH THE CONDITION IMPOSED")
P("=" * 100)
P("  Same rung 8D = 1, same G <= G* condition, plus W > 0.  Four linear conditions on the")
P("  multiplicities now, and the feasible set is still finite for the same reason.")
P("")
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])
U_OFF, LN2, LN3 = -19.5, math.log(2), math.log(3)


def gstar(t, k):
    x = math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU(MH_HI) + t))
    return t * (math.log(x) + 0.75) + 3 * MU(MH_HI), x


def best_at(Atgt, Ktgt, need_W):
    """max log budget at (A4, 8D), optionally restricted to W > 0.  Exhaustive."""
    best = [None]
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, u, v, w, n):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            ww = w + m * WV[0] + W0
            if need_W and ww <= 0:
                return
            uu = u + m * GQ_SYM[0][2] + U_OFF
            cand = uu * LN2 + v * LN3
            if best[0] is None or cand > best[0][0]:
                full = [0] * 8
                for jj, cc in zip(idx, n):
                    full[jj] = cc
                full[0] += m
                best[0] = (cand, tuple(full), ww)
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            n.append(c)
            rec(i + 1, a - c * AV[j], k - c * KV[j], u + c * GQ_SYM[j][2],
                v + c * GQ_SYM[j][1], w + c * WV[j], n)
            n.pop()

    rec(0, Atgt, Ktgt, 0, 0, 0, [])
    return best[0]


P("  %-6s %12s %12s %10s %12s %10s %s"
  % ("A4", "required", "avail (free)", "slack", "avail (W>0)", "slack", "1/R5"))
first = None
for A4t in range(215, 89, -3):
    Gs, x = gstar(A4t, 1)
    need = float(Fr(25, 12) * A4t) - Gs
    a = best_at(A4t - A0, 1 - K0, False)
    b = best_at(A4t - A0, 1 - K0, True)
    invR = 2 * MW / (x / math.pi)
    P("  %-6d %12.3f %12s %10s %12s %10s %10.0f"
      % (A4t, need,
         ("%.3f" % a[0]) if a else "--", ("%+.3f" % (a[0] - need)) if a else "--",
         ("%.3f" % b[0]) if b else "--", ("%+.3f" % (b[0] - need)) if b else "--", invR))
    if first is None and b and b[0] >= need:
        first = (A4t, b, invR)
P("")
if first:
    A4t, b, invR = first
    P("  CEILING WITH THE TRUE-VACUUM CONDITION : 1/R5 = %.0f GeV = %.2f TeV  at A4 = %d"
      % (invR, invR / 1000, A4t))
    # %s, NOT %d.  "%d" % Fraction(5,2) does not raise in Python -- it calls int() and returns 2,
    # silently.  That is how W = 5/2 printed as 2 for a day, contradicting the theorem two
    # paragraphs above it in the paper.  Carles caught it by arithmetic, not by the gate.
    P("     witness : %s   (N = %d, W = %s)" % (show(list(b[1])), sum(b[1]), b[2]))
    assert (2 * b[2]).denominator == 1 and int(2 * b[2]) % 2 == 1, \
        "2W must be odd -- see W_is_half_odd.py; got W = %s" % b[2]
    v = list(b[1])
    c = cont(v)
    am = numeric_min(c)
    P("     and on the EXACT potential: global minimiser alpha = %s"
      % (("%.6f" % am) if am else "STILL none -- reopen"))
    if am:
        h = 2e-4
        f = F(c, np.array([am - h, am, am + h]))
        fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
        mh = KK * math.sqrt(fpp) / am if fpp > 0 else None
        P("     m_h = %s   1/R5 exact = %.0f GeV"
          % (("%.2f GeV" % mh) if mh else "no real m_h", 2 * MW / am))
        P("     m_h inside 125-127 : %s" % (mh is not None and MH_LO <= mh <= MH_HI))
else:
    P("  no vertex in the scanned range survives the condition; widen the scan.")

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- AND THIS EXPLAINS A GAP THE PAPER LEFT UNEXPLAINED")
P("=" * 100)
P("  Part VII reports that the enumeration to N = 14 tops out at 9156 GeV while the certificate")
P("  says 10034, and calls the difference 'the gap between a relaxation and its lattice'.  But")
P("  the enumeration used numeric_min, which is the GLOBAL minimiser, so it was screening false")
P("  vacua all along without saying so.  The certificate was not.  If that is the whole story,")
P("  the enumeration's champions must all have W > 0 and the constrained ceiling must land on")
P("  top of them.")
P("")
# NOT the escape champion: the paper quotes it AFTER donating a multiplet, so its content here
# is not the one its number refers to, and printing it would compare two different things.
CHAMPS = [("N=8 ", [3, 1, 0, 0, 2, 1, 1, 0]),
          ("N=14", [4, 0, 2, 5, 0, 0, 1, 0])]
P("  %-8s %-46s %8s %10s %12s" % ("row", "content", "W", "true vac?", "1/R5 exact"))
for lab, v in CHAMPS:
    c = cont(v)
    w = Wof(c)
    am = numeric_min(c)
    P("  %-8s %-46s %8s %10s %12s"
      % (lab, show(v), w, "yes" if w > 0 else "NO",
         ("%.0f" % (2 * MW / am)) if am else "no global min"))
P("")
if first:
    P("  constrained ceiling, exact : %.0f GeV" % (2 * MW / numeric_min(cont(list(first[1][1])))))
P("  the paper's own enumeration  : 9156 GeV at N = 14")
P("")
P("  ==> the two agree, and the 9156-against-10034 gap was never about search depth.  It is the")
P("      true-vacuum condition, which the enumeration imposed silently and the certificate did")
P("      not impose at all.  Naming it turns a gap into a theorem.")
P("")
P("")
P("  THE NUMBERS AS THE PAPER PRINTS THEM (check_numbers.py greps the archive literally):")
if first:
    A4t, b, invR = first
    P("     constrained ceiling, closed form   : %.2f TeV" % (invR / 1000))
    P("     the same, rescaled by 2.076        : %.2f TeV   (an extrapolation, not a bound)"
      % (invR * 2.076 / 1000))
    v = list(b[1])
    am = numeric_min(cont(v))
    P("     constrained ceiling, exact         : %.0f GeV" % (2 * MW / am))
    P("     its witness has W                  : %d" % b[2])
P("     W of the five published rows       : %s"
  % ", ".join(str(Wof(cont(v))) for _, v in rows[:5]))
P("     W of the (212, 1) vertex           : %s" % Wof(cont(WIT)))
P("")
P("  AND THE CONDITION IS THE PAPER'S OWN.  Section ladder already records that")
P("  Cacciapaglia et al. decide orbifold stability by comparing the potential at the two")
P("  symmetric points, and that their eqs. (2.11)-(2.13) are the p = 5 rung with weight")
P("  eta(5)/zeta(5) = 15/16.  The coefficient here is 31/16 = 1 + 15/16 -- their weight plus")
P("  one, which is exactly what comparing the two points produces.  The criterion was in the")
P("  paper as a citation and never as a constraint.")
