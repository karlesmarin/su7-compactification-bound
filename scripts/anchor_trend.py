#!/usr/bin/env python3
"""Where does the Komori-Maru residual live?  It is not the 48, and it is not a constant.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part VI section 7 reports the anchor residual as two GROUP MEANS -- 1.94 for the rows carrying a
48 and 1.20 for those without -- and Part VII inherits it as a band.  Read that way it invites
the diagnosis "the 48 is mistranscribed".  This script asks the five rows directly, and that
diagnosis does not survive them:

  H1  the residual is a global normalisation      -> all five ratios equal.
  H2  the residual is the 48's charge/parity      -> the three rows WITHOUT a 48 agree.
  H3  the residual is a function of the content   -> ratios move with a moment.

H1 and H2 are falsifiable here and both die: the no-48 rows are 1.286, 1.025, 1.286, which are
not equal, and no 48 appears in any of them.  H3 stands -- the observed pairing is one of only 2
monotone pairings out of 120, p = 0.017.

But which moment is NOT settled, and the control says so: A_4 orders the five rows exactly as well
as D does.  Over these rows the two moments are nearly collinear, and five points cannot separate
them.  Reporting "monotone in D" would be reporting the half of the result that happens to be
interesting.

What does not depend on which moment it is, and is the reason this was worth doing: the certified
ceiling of CEILING.md sits at D = 1/8, NINE TIMES below the smallest anchored D, at the end of the
trend where the disagreement is largest and still growing.  The band quoted with it (1.20-1.94)
was measured at D in [1.125, 3.625].  So the ceiling is quoted outside the range that anchors it.
That is scope, not error: nothing here says the ceiling is wrong.

It does, however, kill a claim of ours: a residual that were a global constant would cancel in the
RATIO 10.03/3.97 = 2.53, so that ratio would be anchor-free.  H1 is dead, so it does not cancel,
and 2.53 inherits whatever the residual does between D = 1/8 and D = 11/8.  Both endpoints are
below the anchored range.

Controls, both of which must fail:
  K1  a shuffled pairing of ratios to rows must destroy the monotonicity
  K2  A_4 as the explanatory moment instead of D

K1 caught a broken criterion before it caught anything else: see the tie tolerance in monotone().
"""
import itertools
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

P("=" * 100)
P("THE FIVE ROWS, and what the ratio is a function of")
P("=" * 100)

rows = []
for label, content, a_them, mh_them, invR in T1:
    a = numeric_min(content)
    mo = moments(content)
    n48 = sum(m for rep, m, _, _ in
              [(r, 1, 0, 0) for r, _, _, _ in content]) if False else \
        sum(1 for rep, _, _, _ in content if rep == "48")
    rows.append(dict(label=label, D=mo["D"], A4=mo["A4"], ratio=a / a_them, n48=n48))

P("%-6s %8s %8s %6s %10s" % ("row", "D", "A4", "has 48", "ratio"))
for r in rows:
    P("%-6s %8.4f %8.0f %6s %10.3f" % (r["label"], r["D"], r["A4"], "yes" if r["n48"] else "no",
                                       r["ratio"]))

P("")
P("=" * 100)
P("H1 -- a global normalisation: all five ratios equal")
P("=" * 100)
lo, hi = min(r["ratio"] for r in rows), max(r["ratio"] for r in rows)
P("    spread %.3f - %.3f, i.e. a factor %.2f across the five rows" % (lo, hi, hi / lo))
P("    VERDICT: %s" % ("holds" if hi / lo < 1.01 else "*** DEAD *** -- the residual is not a constant"))

P("")
P("=" * 100)
P("H2 -- the 48: then the rows WITHOUT one must agree")
P("=" * 100)
wo = [r for r in rows if not r["n48"]]
w = [r for r in rows if r["n48"]]
lo2, hi2 = min(r["ratio"] for r in wo), max(r["ratio"] for r in wo)
P("    with a 48    (n=%d): %s   mean %.2f" % (len(w), ", ".join("%.3f" % r["ratio"] for r in w),
                                               sum(r["ratio"] for r in w) / len(w)))
P("    without      (n=%d): %s   mean %.2f" % (len(wo), ", ".join("%.3f" % r["ratio"] for r in wo),
                                               sum(r["ratio"] for r in wo) / len(wo)))
P("    the no-48 rows spread %.3f - %.3f, a factor %.2f, with no 48 anywhere in them" %
  (lo2, hi2, hi2 / lo2))
P("    VERDICT: %s" % ("holds" if hi2 / lo2 < 1.01 else
                       "*** DEAD *** -- the residual moves where the 48 cannot reach"))

P("")
P("=" * 100)
P("H3 -- a function of the content: sort by D and by A_4 and see which one orders the ratios")
P("=" * 100)


def monotone(pairs, tol=0.01):
    """pairs sorted by the candidate moment; ties in the moment must not break it.

    The tie tolerance is not cosmetic and it was wrong once: rows (1) and (5) carry the SAME D and
    ratios that differ in the fourth decimal, so an exact tie test called the true ordering a
    violation -- and the control then reported 0 of 120 pairings monotone, which is the signature
    of a criterion no hypothesis can pass rather than of a hypothesis that failed.  A tie in x is
    honoured if the two y agree to `tol` relative.
    """
    s = sorted(pairs)
    bad = 0
    for (xa, ya), (xb, yb) in zip(s, s[1:]):
        if abs(xa - xb) < 1e-12:
            bad += abs(ya - yb) / max(abs(ya), abs(yb)) > tol
        else:
            bad += not (yb <= ya * (1 + tol))    # decreasing
    return bad


for name, key in (("D", "D"), ("A_4", "A4")):
    pairs = [(r[key], r["ratio"]) for r in rows]
    bad = monotone(pairs)
    P("    by %-4s: %s" % (name, "  ".join("%.4g->%.3f" % p for p in sorted(pairs))))
    P("    %-9s %d violation(s) of 'decreasing, ties included'   %s"
      % ("", bad, "MONOTONE" if bad == 0 else "not monotone"))

P("")
P("=" * 100)
P("CONTROLS -- each must fail, or the test above measures nothing")
P("=" * 100)
d_pairs = [(r["D"], r["ratio"]) for r in rows]
vals = [r["ratio"] for r in rows]
shuffles = [p for p in itertools.permutations(vals)]
bad_counts = [monotone(list(zip([r["D"] for r in rows], p))) for p in shuffles]
clean = sum(1 for b in bad_counts if b == 0)
P("    K1  all %d re-pairings of the five ratios to the five D: %d are monotone."
  % (len(shuffles), clean))
P("        the observed pairing is one of them, so the odds of this by chance are %d/%d = %.3f"
  % (clean, len(shuffles), clean / len(shuffles)))
P("    K2  A_4 as the explanatory moment: %d violation(s) -- %s"
  % (monotone([(r["A4"], r["ratio"]) for r in rows]),
     "it does NOT order them" if monotone([(r["A4"], r["ratio"]) for r in rows]) else
     "*** it orders them too, so D is not singled out ***"))

P("")
P("=" * 100)
P("WHERE THE CEILING SITS RELATIVE TO THE ANCHORED RANGE")
P("=" * 100)
Ds = sorted(r["D"] for r in rows)
P("    the five rows anchor D in [%.3f, %.3f]" % (Ds[0], Ds[-1]))
P("    the certified ceiling of CEILING.md sits at D = 1/8 = 0.125, which is %.1fx BELOW the"
  % (Ds[0] / 0.125))
P("    smallest anchored D, and the trend is still rising there.")
P("")
P("    ratio at the smallest anchored D (%.3f) : %.3f" % (Ds[0], max(vals)))
P("    ratio at the largest             (%.3f) : %.3f" % (Ds[-1], min(vals)))
P("")
P("    So the ceiling is quoted at the one corner of the content space where the anchor is worst")
P("    and unmeasured.  That is a statement about scope, not a correction: nothing here says the")
P("    ceiling is wrong.  It says the band attached to it was read off the wrong statistic.")
