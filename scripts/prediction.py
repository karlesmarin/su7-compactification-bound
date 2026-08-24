#!/usr/bin/env python3
"""What the closed form PREDICTS, and against what it can be falsified.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  Collects the falsifiable content of Part VII: the compactification-scale ceiling with its
  anchor band, the m_h - 1/R_5 correlation the closed form forces, and the arithmetic tests any
  future table of this kind must pass.

Three kinds of prediction, in decreasing order of how hard they are to wriggle out of:

  A  ARITHMETIC, and anchor-free.  8D is an odd integer over 8, and 8D + A_4 = 0 (mod 3).  Both are
     properties of the CONTENT alone -- no normalisation, no loop order, no scheme.  Any table
     published for this model class must satisfy them, so they are a check anyone can run.

  B  CORRELATION.  alpha_min and m_h are both functions of the same two moments, so they are not
     independent data: pinning m_h to its measured value FIXES 1/R_5 from the content, with no
     minimisation.  That is a prediction per content, and it is what makes C possible.

  C  A CEILING, with experimental teeth.  1/R_5 <= 10.03 TeV for ANY bulk content.  The model class
     therefore cannot hide above that scale.  Reported here with the Part VI anchor band attached,
     because every absolute number in this series inherits it.
"""
import json
import math
import pathlib

import numpy as np

exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW = 80.4
CEIL_OURS = 10034.0                      # ceiling_ilp.py, m_h <= 127 GeV, g4 = 0.63

P("")
P("=" * 100)
P("A -- THE ARITHMETIC TESTS, which no normalisation can touch")
P("=" * 100)
P("%-5s %8s %8s %14s %14s %10s" % ("row", "8D", "A4", "8D odd?", "8D+A4 = 0 (3)?", "verdict"))
for label, content, a_them, mh_them, invR in T1:
    mo = moments(content)
    k, t = round(8 * mo["D"]), round(mo["A4"])
    ok = (k % 2 == 1) and ((k + t) % 3 == 0)
    P("%-5s %8d %8d %14s %14s %10s" % (label, k, t, k % 2 == 1, (k + t) % 3 == 0,
                                       "passes" if ok else "*** FAILS ***"))
P("")
P("Their own five rows pass both.  A future table that failed either would be arithmetically")
P("impossible in this model class, whatever its authors' conventions -- which is the point.")

P("")
P("=" * 100)
P("B -- THE CORRELATION: m_h and 1/R_5 are not independent")
P("=" * 100)
P("With m_h pinned and mu = m_h^2/(K pi^2)^2,   x^2 = 12 z3 D/(6 mu + A_4),  so")
P("")
P("     1/R_5 = 2 pi m_W sqrt( (6 mu + A_4) / (12 z3 D) )     -- content in, scale out")
P("")
P("%-5s %8s %8s %13s %13s %13s" % ("row", "D", "A4", "1/R5 theirs", "1/R5 ours", "ratio"))
rows = []
for label, content, a_them, mh_them, invR in T1:
    a = numeric_min(content)
    mo = moments(content)
    rows.append((label, a_them, a, 2 * MW / a_them, 2 * MW / a))
    P("%-5s %8.4f %8.0f %13.0f %13.0f %13.3f" %
      (label, mo["D"], mo["A4"], 2 * MW / a_them, 2 * MW / a, a / a_them))
r_lo = min(r[2] / r[1] for r in rows)
r_hi = max(r[2] / r[1] for r in rows)
P("")
P("The anchor band, stated exactly: the ROW-WISE range of alpha_ours/alpha_theirs is %.2f-%.2f." %
  (r_lo, r_hi))
_w48 = [r for r, (lab, ct, *_ ) in zip(rows, T1) if any(rep == "48" for rep, _, _, _ in ct)]
_wo = [r for r, (lab, ct, *_ ) in zip(rows, T1) if not any(rep == "48" for rep, _, _, _ in ct)]
P("Part VI section 7 quotes the GROUP MEANS, which are a different statistic: %.2f with a 48 (n=%d)"
  % (sum(r[2] / r[1] for r in _w48) / len(_w48), len(_w48)))
P("and %.2f without (n=%d).  Every ABSOLUTE scale below carries this factor; it is the one number"
  % (sum(r[2] / r[1] for r in _wo) / len(_wo), len(_wo)))
P("in this paper that is not settled, and the bands quoted use the row-wise range, not the means.")

P("")
P("=" * 100)
P("C -- THE CEILING, and what it costs the model to survive")
P("=" * 100)
P("  ceiling in OUR normalisation, any content, m_h <= 127 GeV : %8.0f GeV" % CEIL_OURS)
P("  the same, rescaled by the anchor band                     : %8.0f - %.0f GeV" %
  (CEIL_OURS * r_lo, CEIL_OURS * r_hi))
P("  ceiling at their own content size (at most 6 multiplets)  :     2182 GeV  (ceiling_ilp.py)")
P("  their own published 1/R_5 column                          : %s GeV" %
  ", ".join("%.0f" % (2 * MW / r[1]) for r in rows))
P("")
P("Read against the collider bounds this scale controls -- the KK excitations of the gauge bosons")
P("sit at 1/R_5 -- the statement is sharp in both directions:")
P("")
P("  * ATLAS/CMS exclude a KK gluon below about 4 TeV and a sequential Z' below about 5 TeV.")
P("    A content of their own size gives at most %s GeV in our normalisation and %s GeV even at" %
  ("2182", "%.0f" % (2182 * r_hi)))
P("    the top of the anchor band -- BELOW the published exclusions.  The small-content corner of")
P("    the model class is therefore already disfavoured, and that is a consequence of the closed")
P("    form, not of any new experiment.")
P("  * the ceiling caps the whole class at %.1f TeV (%.1f TeV at the top of the band).  There is no" %
  (CEIL_OURS / 1000, CEIL_OURS * r_hi / 1000))
P("    corner of the multiplet lattice that puts the compactification scale beyond that, so a")
P("    machine that reaches it settles the class.  FCC-hh does; HL-LHC does not.")
P("")
P("  >> The prediction is bounded ABOVE, which is the useful direction: the model cannot retreat")
P("     to higher scales by adding matter.  Adding matter raises A_4, and A_4 raises 1/R_5 only")
P("     until the Higgs window closes it off -- which is exactly what the integer program measures.")

P("")
P("=" * 100)
P("D -- AND A PREDICTION IN SOMEONE ELSE'S FRAMEWORK")
P("=" * 100)
P("Cacciapaglia et al. (arXiv:2409.16137) decide orbifold stability by comparing the potential at")
P("a = 0 and a = 1/2, with weight eta(5)/zeta(5) = 15/16.  Our rung k=0 is that criterion, and in")
P("this model it reads  32 F(0)/zeta(5) = 9 + even, always ODD.  Hence")
P("")
P("     |F(0)| >= zeta(5)/32   for every bulk content, and the bound is attained.")
P("")
P("Their criterion, applied to a model they did not consider, is NEVER marginal.  That is a")
P("statement in their language, about their question, that they can check against their own code.")

OUT.mkdir(exist_ok=True)
(OUT / "prediction.json").write_text(json.dumps(dict(
    ceiling_ours=CEIL_OURS, anchor_band=[r_lo, r_hi],
    ceiling_band=[CEIL_OURS * r_lo, CEIL_OURS * r_hi],
    ceiling_their_size=2182.0,
    rows=[dict(row=r[0], alpha_theirs=r[1], alpha_ours=r[2],
               invR_theirs=r[3], invR_ours=r[4]) for r in rows]), indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "prediction.json"))
