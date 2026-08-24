#!/usr/bin/env python3
"""necessity_of_W.py -- W > 0 is NECESSARY exactly where it has to be, so 9.22 TeV is a theorem.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE GAP THIS CLOSES
-------------------
The paper bounds 1/R5 by 9.22 TeV over contents satisfying W > 0, and says honestly that W > 0 has
only been shown SUFFICIENT for the electroweak point to be the true vacuum.  Until it is also
necessary, 9.22 is a bound on a screened family and not on all of them.

Carles's route, and it does not require proving necessity in general: a counterexample only matters
if it BEATS 9.22 TeV.  Such a content has

    alpha_EW < 2 m_W / 9220 GeV = 0.01744 ,   i.e.  x = pi alpha < 0.05479 .

Now put the two facts against each other.

  THE GAP AT alpha = 1 IS QUANTISED AND CANNOT BE SMALL.  Since 2W is an odd integer,
  |F(1) - F(0)| = (31/16) zeta(5) |W| >= (31/32) zeta(5) = 1.0045.  A content with W < 0 therefore
  has F(1) <= F(0) - 1.0045 -- an order-one drop.

  THE ELECTROWEAK WELL IS MICROSCOPIC.  Eliminating G with stationarity and then A_4 with (II),
  the depth collapses to a two-term expression:

      F(x*) - F(0) = - zeta(3) D x*^2 / 8  -  mu x*^4 / 16                     (derived below)

  which at x* < 0.055 is of order 1e-4 for D = 1/8.

So if the well is shallower than the quantum gap, W < 0 forces F(1) < F(x*): the electroweak point
is not the minimum, and the content is not a counterexample after all.  Hence every content that
could beat 9.22 TeV has W > 0 -- and the integer program over W > 0 already says none does.

This file derives the depth formula, checks it against the exact potential, and computes how far
the argument reaches.

Run:  python necessity_of_W.py
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

QUANTUM = Fr(31, 32) * Z5

P("=" * 100)
P("1 -- THE DEPTH OF THE ELECTROWEAK WELL, IN CLOSED FORM")
P("=" * 100)
P("  F(x) - F(0) = -zeta(3) D x^2/2 + x^4 (G - A_4 ln x)/24 + R.")
P("  Stationarity (eq. amin) gives  G - A_4 ln x = [24 zeta(3) D / x^2 + A_4]/4, so the quartic")
P("  term becomes  zeta(3) D x^2/4 + A_4 x^4/96  and the depth collapses to")
P("")
P("      F(x*) - F(0) = - zeta(3) D x*^2/4 + A_4 x*^4/96 .")
P("")
P("  Then (II), which holds once m_h is pinned, gives 6 mu + A_4 = 12 zeta(3) D / x*^2, hence")
P("  A_4 = 12 zeta(3) D/x*^2 - 6 mu, and substituting:")
P("")
P("      F(x*) - F(0) = - zeta(3) D x*^2 / 8  -  mu x*^4 / 16 .")
P("")
P("  Two terms, both negative, and the second is negligible.  Checked against the exact")
P("  potential on the five published rows and on the constrained witness:")
P("")


def depth_formula(D, mu, x):
    return -Z3 * D * x ** 2 / 8 - mu * x ** 4 / 16


P("  %-40s %13s %13s %10s" % ("content", "formula", "exact", "err %"))
rows = []
for label, content, a_them, mh_them, invR in T1:
    v = [0] * 8
    for rep, e, ep, mult in content:
        v[REPS.index((rep, e, ep))] += mult
    rows.append((label, v))
W104 = [0] * 8
W104[0], W104[3], W104[4], W104[5], W104[6] = 16, 1, 1, 4, 1
rows.append(("witness (104,1)", W104))
for label, v in rows:
    c = cont(v)
    a = numeric_min(c)
    if a is None:
        P("  %-40s %13s %13s %10s" % (label, "--", "--", "no global min"))
        continue
    x = math.pi * a
    mo = moments(c)
    h = 2e-4
    f = F(c, np.array([a - h, a, a + h]))
    fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
    mh = KK * math.sqrt(fpp) / a if fpp > 0 else None
    mu = MU(mh) if mh else float("nan")
    ex = float(F(c, np.array([a]))[0] - F(c, np.array([1e-9]))[0])
    fo = depth_formula(mo["D"], mu, x)
    P("  %-40s %13.6f %13.6f %10.2f" % (label, fo, ex, 100 * (fo - ex) / abs(ex)))

P("")
P("=" * 100)
P("2 -- THE TWO SCALES, SIDE BY SIDE")
P("=" * 100)
X0 = math.pi * 2 * MW / 9220.0
P("  a counterexample must beat 9.22 TeV, so alpha < 2 m_W/9220 = %.6f, i.e. x < %.6f"
  % (2 * MW / 9220.0, X0))
P("  the quantised gap at alpha = 1 : |F(1) - F(0)| >= (31/32) zeta(5) = %.6f" % float(QUANTUM))
P("")
P("  %-8s %16s %16s %s" % ("8D", "well depth", "quantum gap", "well shallower?"))
mu_hi = MU(MH_HI)
for k in (1, 3, 9, 51, 501, 5001, 17817, 20000):
    d = abs(depth_formula(k / 8.0, mu_hi, X0))
    P("  %-8d %16.6f %16.6f %s" % (k, d, float(QUANTUM), d < float(QUANTUM)))
P("")
kmax = 8 * (float(QUANTUM) - mu_hi * X0 ** 4 / 16) * 8 / (Z3 * X0 ** 2)
P("  the argument reaches every rung with 8D < %.0f." % kmax)
P("  And it does not need to reach further: the certificate bounds rung 3 at 6.27 TeV and every")
P("  rung above it lower still, so no rung other than 8D = 1 can beat 9.22 TeV at all.")

P("")
P("=" * 100)
P("3 -- THE ARGUMENT, PUT TOGETHER")
P("=" * 100)
d1 = abs(depth_formula(1 / 8.0, mu_hi, X0))
P("  Let C be a content with 1/R5 > 9.22 TeV whose electroweak point is the true vacuum.")
P("  (a) Only the rung 8D = 1 can exceed 9.22 TeV, so D = 1/8.")
P("  (b) Its well is at most %.6f deep, against a quantum gap of %.6f -- a factor %.0f."
  % (d1, float(QUANTUM), float(QUANTUM) / d1))
P("  (c) If W < 0 then F(1) <= F(0) - %.6f < F(0) - %.6f <= F(x*), so the electroweak point is"
  % (float(QUANTUM), d1))
P("      NOT the global minimum, contradicting the assumption.")
P("  (d) W = 0 is arithmetically impossible (2W is odd).")
P("  (e) Hence W > 0.  And the integer program over W > 0 finds no content above 9.22 TeV.")
P("")
P("  ==> W > 0 IS NECESSARY IN THE ONLY REGION WHERE IT MATTERS, so 1/R5 <= 9.22 TeV holds for")
P("      ANY bulk content whose electroweak point is its true vacuum -- not merely for the")
P("      screened family.  The screen was never a restriction where the bound is decided.")
P("")
P("  CONTROLS, and both can fail:")
P("     the depth formula must reproduce the exact well : see section 1")
P("     the factor of safety must be large, not marginal : %.0f" % (float(QUANTUM) / d1))
P("")
P("  WHAT IS STILL ASSUMED, said plainly: step (a) uses the per-rung certificate, which is")
P("  itself global only under the measured monotonicity of t*(k)/k.  This argument removes the")
P("  W > 0 caveat from 9.22 TeV; it does not remove that one.")

import json
(HERE / "outputs" / "necessity_of_W.json").write_text(json.dumps(dict(
    x0=X0, alpha0=2 * MW / 9220.0, quantum=float(QUANTUM),
    depth_at_rung1=d1, safety=float(QUANTUM) / d1, kmax=kmax), indent=2), encoding="utf-8")
P("")
P("  data written: outputs/necessity_of_W.json")
