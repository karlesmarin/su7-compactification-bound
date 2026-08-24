#!/usr/bin/env python3
"""anchor_delta.py -- the residual should live where the working anchor cannot see, so look there.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE PRIORITY QUESTION.  Every absolute number in this paper carries the same caveat: our F does not
reproduce \\cite{KM25}'s published alpha_min column, by a factor of 1.94 with a 48 and 1.20 without.
Everything else measured here -- a 0.39 % coupling, a 2 % matching scale, a 1.4 % Higgs window --
lives inside that.  So the question that decides what the next session is for is not "what is g_4"
but "what would have to be true for the absolute numbers to mean anything".

AND THERE IS A DEDUCTION AVAILABLE THAT NOBODY HAS DRAWN.  Two anchors have been tried:

  vGIQ, hep-th/0204223 -- five-dimensional, SINGLE PARITY.  It WORKS.  Their criterion (5.3) is
  our D > 0; the three critical N_f come out 3, 9/2, 3/4 exactly by two routes; the SU(2) adjoint
  minimum is 1/4 exactly; the whole group-theory bracket of their m_h^2 including the 3:4.  So F,
  the expansion, D and the minimiser are validated against an independent published computation.

  \\cite{KM25} -- two parities.  It MISSES, by 1.94 and 1.20.

Single parity means Delta_k = sum m s c^{2k} does not mix signs; the working anchor sits where the
parity-weighted moment cannot do anything.  If the machinery is validated there and fails here,
THE RESIDUAL MUST LIVE IN SOMETHING THAT VANISHES WITH Delta.  That is a localisation, and it is
testable on the five published rows.

WHAT anchor_trend.py DID AND DID NOT DO.  It ordered the five ratios by D and by A_4, found both
monotone, and said honestly that five points cannot separate two nearly collinear moments.  It never
tried Delta.  Delta is a different combination of the same two traces -- D = A_2 - (3/4) B_2 against
Delta_1 = A_2 - B_2 -- so it is not obviously collinear with either, and it is the one the vGIQ
result points at.

Sections:
  1  the five rows with Delta at every rung, beside D and A_4
  2  does Delta order the ratios, and is it collinear with D over these rows?
  3  the extrapolation the localisation would license, and why it is not licensed yet
  4  controls

Run:  python anchor_delta.py
"""
import itertools
import json
import math
import pathlib
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]


def SD(content, k):
    """S_k, the ladder's Delta_k, and B_k -- and the three are not two.

    THE PAPER USES Delta FOR TWO DIFFERENT OBJECTS, which is what sent the first version of this
    file after the wrong variable.  Eq. (rung) defines Delta_k = sum m s c^{2k}, the
    parity-weighted moment.  The anchor sections use Delta for the ANTIPERIODIC CONTENT: section
    parity says "in the single-parity corner every tower is integer-moded, Delta = 0, and
    8D = 8A_2", and section anchor says "D = A_2 - (3/4) B_2 collapses to D = A_2".  Both of
    those are B = 0, not sum m s c^{2k} = 0 -- with a single parity the ladder's Delta_k equals
    S_k and is as far from zero as it can be.

    The variable the vGIQ localisation is about is therefore B_k, the antiperiodic trace.
    """
    tab = table(content)
    S = sum(m * float(c) ** (2 * k) for m, s, c in tab)
    Dl = sum(m * s * float(c) ** (2 * k) for m, s, c in tab)
    B = sum(m * float(c) ** (2 * k) for m, s, c in tab if s < 0)
    return S, Dl, B


# ================================================================= 1
P("=" * 100)
P("1 -- THE FIVE PUBLISHED ROWS, WITH Delta BESIDE THE MOMENTS ALREADY TRIED")
P("=" * 100)
rows = []
for label, content, a_them, mh_them, invR in T1:
    a_ours, mo = closed_form(content)
    n48 = sum(mult for rep, e, ep, mult in content if rep == "48")
    r = dict(label=label, ratio=a_ours / a_them, D=mo["D"], A4=mo["A4"], n48=n48,
             a_ours=a_ours, a_them=a_them)
    for k in (0, 1, 2):
        s, d, b = SD(content, k)
        r["S%d" % k], r["Dl%d" % k], r["B%d" % k] = s, d, b
    rows.append(r)

P("  %-5s %9s %9s %8s %7s %10s %10s %10s" %
  ("row", "a ours", "a theirs", "ratio", "has 48", "D", "A_4", "Delta_1"))
for r in rows:
    P("  %-5s %9.5f %9.5f %8.3f %7s %10.4f %10.0f %10.1f"
      % (r["label"], r["a_ours"], r["a_them"], r["ratio"], "yes" if r["n48"] else "no",
         r["D"], r["A4"], r["Dl1"]))
P("")
P("  %-5s %11s %11s %11s   %11s %11s %11s"
  % ("row", "ladder D_0", "ladder D_1", "ladder D_2", "B_0", "B_1", "B_2"))
for r in rows:
    P("  %-5s %11.1f %11.1f %11.1f   %11.1f %11.1f %11.1f"
      % (r["label"], r["Dl0"], r["Dl1"], r["Dl2"], r["B0"], r["B1"], r["B2"]))
P("")
P("  The left block is the ladder's Delta_k = sum m s c^{2k}; the right is B_k, the antiperiodic")
P("  trace.  B is the one that vanishes in the corner where the anchor works -- see SD() -- and")
P("  it is the one the localisation is about.  Both are shown because the paper writes Delta for")
P("  both, which is a collision worth fixing there.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- DOES Delta ORDER THE RATIOS?  and is it collinear with D over these rows?")
P("=" * 100)


def monotone(pairs, tol=0.01):
    """anchor_trend.py's test, in behaviour and for its reason.

    The first version of this file rewrote it with an exact tie test, and rows (1) and (5) --
    same D, ratios differing in the fourth decimal -- then made EVERY variable come out
    non-monotone, D included.  That is the signature of a broken criterion, not of a failed
    hypothesis; anchor_trend.py had already been bitten by it and left the fix in a docstring.
    A tie in x is honoured if the two y agree to `tol` relative.
    """
    s = sorted(pairs)
    bad = 0
    for (xa, ya), (xb, yb) in zip(s, s[1:]):
        if abs(xa - xb) < 1e-12:
            bad += abs(ya - yb) / max(abs(ya), abs(yb)) > tol
        else:
            bad += not (yb <= ya * (1 + tol))
    return bad == 0


def perm_p(key):
    """the exact permutation test anchor_trend.py uses: how many of the 120 re-pairings are
    monotone against this ordering variable."""
    xs = [r[key] for r in rows]
    ys = [r["ratio"] for r in rows]
    good = sum(1 for p in itertools.permutations(ys) if monotone(list(zip(xs, p))))
    return good, len(list(itertools.permutations(ys)))


P("  %-10s %10s %12s %14s %10s" % ("variable", "monotone", "of 120", "p", "corr with D"))
res = {}
for name, key in (("D", "D"), ("A_4", "A4"), ("B_0", "B0"), ("B_1", "B1"), ("B_2", "B2"),
                  ("ladder D_1", "Dl1"), ("S_1", "S1")):
    g, n = perm_p(key)
    c = np.corrcoef([r[key] for r in rows], [r["D"] for r in rows])[0, 1]
    res[name] = dict(monotone=monotone([(r[key], r["ratio"]) for r in rows]), good=g, p=g / n,
                     corr_D=float(c))
    P("  %-10s %10s %12d %14.4f %10.3f"
      % (name, res[name]["monotone"], g, g / n, c))
P("")
P("  The p column is the fraction of the 120 re-pairings that would also come out monotone, so")
P("  it is the honest significance of 'ordered by' on five points.  A variable that orders the")
P("  ratios AND is uncorrelated with D would be information; one that orders them because it is")
P("  D in disguise would not.")
P("")
_dl = [n for n in ("B_0", "B_1", "B_2") if res[n]["monotone"]]
if _dl:
    best = min(_dl, key=lambda n: (res[n]["p"], abs(res[n]["corr_D"])))
    P("  Delta orders the ratios at rungs: %s" % ", ".join(_dl))
    P("  the least collinear of them is %s, with correlation %.3f against D."
      % (best, res[best]["corr_D"]))
    if abs(res[best]["corr_D"]) > 0.95:
        P("")
        P("  AND IT ORDERS THEM FOR THE SAME REASON D DOES.  Over these five rows B is collinear")
        P("  with D to better than 0.95, so it separates nothing.  The vGIQ deduction stands as a")
        P("  deduction -- the residual must vanish with the antiperiodic content, because the")
        P("  anchor that works has none -- but the five published rows cannot confirm it, for")
        P("  exactly the reason anchor_trend.py already gave about D and A_4: five points on a")
        P("  nearly one-dimensional family.")
    else:
        P("")
        P("  AND B IS NOT COLLINEAR WITH D over these rows, correlation %.3f.  So it orders the"
          % res[best]["corr_D"])
        P("  ratios for a reason of its own, and that is the first piece of evidence FOR the")
        P("  localisation that the published table can supply.  It is five points and it is not a")
        P("  proof; what it does is tell the next session where to spend its anchor.")
else:
    P("  B does NOT order the ratios at any rung.  That is evidence AGAINST the localisation and")
    P("  it is worth taking at face value: if the residual were a function of the antiperiodic")
    P("  content alone, the antiperiodic content should order it, and it does not.")

# ================================================================= 2b
P("")
P("=" * 100)
P("2b -- AND THE ORDERING WAS THE WRONG INSTRUMENT.  ONE PAIR SAYS MORE THAN FIVE POINTS.")
P("=" * 100)
P("  Rows (1) and (2) have the SAME parity-blind moments to the last digit:")
P("")
P("  %-5s %10s %10s %10s %10s %10s %10s" % ("row", "S_0", "S_1", "S_2", "D", "A_4", "ratio"))
for r in rows:
    if r["label"] in ("(1)", "(2)"):
        P("  %-5s %10.1f %10.1f %10.1f %10.4f %10.0f %10.3f"
          % (r["label"], r["S0"], r["S1"], r["S2"], r["D"], r["A4"], r["ratio"]))
r1 = [r for r in rows if r["label"] == "(1)"][0]
r2 = [r for r in rows if r["label"] == "(2)"][0]
same = all(abs(r1["S%d" % k] - r2["S%d" % k]) < 1e-9 for k in (0, 1, 2))
P("")
P("  identical at every rung : %s        their ratios differ by %.0f %%"
  % (same, 100 * abs(r1["ratio"] - r2["ratio"]) / min(r1["ratio"], r2["ratio"])))
assert same, "rows (1) and (2) do not share their moments -- this section's premise is wrong"
P("")
P("  So the residual is NOT a function of the parity-blind moments.  Whatever it is, it reads")
P("  the parity SPLIT -- and that is the localisation, established by one comparison instead of")
P("  by an ordering over a collinear family.  The ordering test above was the wrong instrument,")
P("  not a negative result.")
P("")
P("  AND THE SPLIT IS A NEAR-CANCELLATION, which is why a small error there is a large error")
P("  here.  D = A_2 - (3/4) B_2 on these two rows:")
P("")
P("  %-5s %12s %12s %14s %12s" % ("row", "A_2", "(3/4) B_2", "D = difference", "D / A_2"))
for r, lab in ((r1, "(1)"), (r2, "(2)")):
    b2 = r["B1"]
    a2 = r["S1"] - b2
    P("  %-5s %12.3f %12.3f %14.4f %11.2f %%" % (lab, a2, 0.75 * b2, a2 - 0.75 * b2,
                                                 100 * (a2 - 0.75 * b2) / a2))
P("")
P("  D is a two-percent residue of A_2, so the one-percent difference in the split between these")
P("  two rows moves D by a factor of %.2f.  And x ~ sqrt(D), so it moves alpha_min by %.0f %% --"
  % (r2["D"] / r1["D"], 100 * (math.sqrt(r2["D"] / r1["D"]) - 1)))
P("  which is the size of the residual itself.  A small error in the antiperiodic trace would")
P("  therefore reproduce this, and that is a testable repair rather than a hope:")
P("")


def repair_lambda(content, a_them):
    """the single factor on every ANTIPERIODIC multiplicity that would put our alpha_min on
    theirs.  One parameter, applied to the whole antiperiodic sector, gauge included."""
    def amin(lam):
        """(status, alpha).  The two failure modes are OPPOSITE ends and must not share a code:
        'dneg' means lambda is too LARGE (D has gone negative), 'nofp' means it is too SMALL
        (D so large that the fixed point runs off).  The first version returned None for both
        and the bisection collapsed onto whichever bracket it happened to touch first, printing
        alpha = nan beside a confident verdict -- once as a false positive and once as a false
        negative.  [[a-parse-failure-is-not-a-verdict]]"""
        tab = [(m * (lam if s < 0 else 1.0), s, c) for m, s, c in table(content)]
        A2 = sum(m * float(c) ** 2 for m, s, c in tab if s > 0)
        B2 = sum(m * float(c) ** 2 for m, s, c in tab if s < 0)
        A4 = sum(m * float(c) ** 4 for m, s, c in tab if s > 0)
        B4 = sum(m * float(c) ** 4 for m, s, c in tab if s < 0)
        A4L = sum(m * float(c) ** 4 * math.log(float(c)) for m, s, c in tab if s > 0)
        D_, G_ = A2 - 0.75 * B2, A4 * H4 - A4L - LN2 * B4
        if D_ <= 0:
            return "dneg", None
        x = 0.1
        for _ in range(500):
            den = 4 * G_ - A4 * (4 * math.log(x) + 1)
            if den <= 0:
                return "nofp", None
            x = math.sqrt(24 * Z3 * D_ / den)
        return "ok", x / math.pi
    # Raising lambda raises B_2, which LOWERS D = A_2 - (3/4) B_2 and so lowers alpha.  It also
    # has a hard wall: D goes negative at lambda = A_2 / ((3/4) B_2), a little above 1, and there
    # amin returns None.  None therefore means OVERSHOT and must move the upper bracket, not the
    # lower -- the first version pushed it the wrong way and every row came back pinned at the
    # bracket with alpha = nan, which is the signature of a search that hit a wall, not of a
    # repair that worked.  [[a-control-that-cannot-fail]]
    lo, hi = 0.5, 3.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        st, v = amin(mid)
        if st == "dneg":            # lambda too large
            hi = mid
        elif st == "nofp":          # lambda too small
            lo = mid
        elif v > a_them:            # alpha still above theirs: raise lambda
            lo = mid
        else:
            hi = mid
    lam = 0.5 * (lo + hi)
    st, v = amin(lam)
    return lam, (v if st == "ok" else None)


P("  %-5s %12s %14s %14s" % ("row", "lambda needed", "alpha at it", "theirs"))
lams = []
for label, content, a_them, mh_them, invR in T1:
    lam, av = repair_lambda(content, a_them)
    lams.append(lam)
    P("  %-5s %12.4f %14.5f %14.5f" % (label, lam, av if av else float("nan"), a_them))
P("")
P("  and the quantity to judge is lambda MINUS ONE, not lambda.  With lambda near unity a")
P("  spread of 0.012 looks tight and is not: the perturbation itself is what varies.")
P("")
eps = [l - 1.0 for l in lams]
P("  %-5s %14s" % ("row", "lambda - 1"))
for (label, content, a_them, mh_them, invR), e in zip(T1, eps):
    P("  %-5s %14.4f" % (label, e))
P("")
P("  lambda      : %.4f to %.4f, a spread of %.4f" % (min(lams), max(lams), max(lams) - min(lams)))
P("  lambda - 1  : %.4f to %.4f, a ratio of %.1f between the largest and the smallest"
  % (min(eps), max(eps), max(eps) / min(eps)))
P("")
if max(eps) / min(eps) < 1.5:
    P("  ONE NUMBER FIXES ALL FIVE -- a genuine one-parameter repair.")
else:
    P("  SO IT IS NOT A COMMON RESCALING.  The required perturbation varies by a factor of %.0f"
      % (max(eps) / min(eps)))
    P("  across the five rows, so a single factor on the antiperiodic sector does not fix them,")
    P("  and this joins the one-parameter repairs anchor_trend.py already excluded.")
P("")
P("  BUT WHAT IT DOES ESTABLISH IS THE SIZE, AND THAT IS THE USEFUL PART.  Every row needs a")
P("  perturbation of the SAME SIGN and of order ONE PERCENT -- %.1f %% to %.1f %% -- to land on"
  % (100 * min(eps), 100 * max(eps)))
P("  their published alpha_min.  A residual that presents as a factor of two in alpha is")
P("  therefore consistent with a sub-percent discrepancy in one trace, amplified by the")
P("  near-cancellation.  That changes what the caveat means: not 'our potential is wrong by")
P("  a factor of two' but 'one trace may be wrong in its second digit, where D has no room'.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- WHAT WOULD SETTLE IT, STATED AS A REQUIREMENT AND NOT AS A HOPE")
P("=" * 100)
P("  The deduction is: the machinery is validated at the vGIQ point and fails at \\cite{KM25}'s,")
P("  so the residual lives in what differs between them.  To turn that into a measurement one")
P("  needs rows that move Delta WITHOUT moving D, and the five published ones do not:")
P("")
P("  %-28s %12s %12s" % ("quantity", "min", "max"))
for name, key in (("D", "D"), ("A_4", "A4"), ("Delta_1", "Dl1"), ("ratio", "ratio")):
    v = [r[key] for r in rows]
    P("  %-28s %12.4f %12.4f" % (name, min(v), max(v)))
P("")
P("  The ratio spans %.2f to %.2f -- a factor of %.2f -- across a family whose three candidate"
  % (min(r["ratio"] for r in rows), max(r["ratio"] for r in rows),
     max(r["ratio"] for r in rows) / min(r["ratio"] for r in rows)))
P("  explanatory variables are mutually collinear.  So the requirement is not another statistic")
P("  on these five points.  It is ONE published potential, by anyone, at two parities, against")
P("  which our F can be evaluated -- the anchor for Delta != 0 that section open already names as")
P("  the sharpest thing missing.  What this file adds is WHY it is the sharpest: not because the")
P("  gap is unexplained, but because the working anchor is blind exactly where the gap lives, so")
P("  no amount of re-reading the vGIQ agreement can ever bound it.")

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- CONTROLS")
P("=" * 100)
P("  (i) THE TEST MUST BE ABLE TO SAY NO.  A variable that cannot order the ratios has to come")
P("      back non-monotone; shuffle the ratios against D and count how many of the 120")
P("      re-pairings survive, which is anchor_trend.py's own control:")
g, n = perm_p("D")
P("      re-pairings of the five ratios to the five D that are monotone : %d of %d" % (g, n))
assert 0 < g < n, "the monotonicity test is degenerate on five points"
P("")
P("  (ii) AND THE ROWS MUST BE THE PUBLISHED ONES.  ratio = our alpha_min over theirs.  The")
P("       paper's 1.94 and 1.20 are the MEANS OF THE TWO GROUPS, not the extremes -- I compared")
P("       them against the extremes first and read a discrepancy that was not there.  It has to")
P("       be the means that come back, or this file is not reading the paper's five rows:")
_w48 = [r["ratio"] for r in rows if r["n48"]]
_wo48 = [r["ratio"] for r in rows if not r["n48"]]
_m48, _mo = sum(_w48) / len(_w48), sum(_wo48) / len(_wo48)
P("       mean with a 48    : %.3f over %d rows   (paper: 1.94)" % (_m48, len(_w48)))
P("       mean without a 48 : %.3f over %d rows   (paper: 1.20)" % (_mo, len(_wo48)))
P("       the extremes, for contrast : %.2f and %.2f -- neither is what the paper quotes"
  % (max(r["ratio"] for r in rows), min(r["ratio"] for r in rows)))
assert abs(_m48 - 1.94) < 0.005 and abs(_mo - 1.20) < 0.005, \
    "the group means are not the paper's 1.94 / 1.20 -- these are not the published rows"
P("")
P("  (iii) AND THE LADDER'S Delta MUST VANISH ON A PARITY PAIR while B does NOT -- which is the")
P("        whole reason the two symbols must not share a letter:")
P("")
P("        %-24s %16s %14s" % ("rep at both parities", "ladder Delta_1", "B_1"))
_g0 = SD([], 1)
okd = okb = True
for rep in ("7", "28", "48", "84"):
    pair = [(rep, 1, 1, 1), (rep, 1, -1, 1)]
    d1 = SD(pair, 1)[1] - _g0[1]
    b1 = SD(pair, 1)[2] - _g0[2]
    okd &= abs(d1) < 1e-9
    okb &= abs(b1) > 1e-9
    P("        %-24s %16.1f %14.1f" % (rep, d1, b1))
P("")
P("        CONTROL -- the ladder's Delta cancels on every pair : %s" % okd)
P("        CONTROL -- B does NOT cancel on any pair            : %s" % okb)
assert okd and okb, "the two Deltas are not distinguishable -- the whole section is confused"
P("        So they are different objects, and this file tests the second one.")

P("")
P("  (iv) AND THE PARITY LABEL p IS NOT THE PERIODICITY s.  This control exists because I tried to")
P("       sharpen the \\cite{PGY24} paragraph with a number and the number was wrong.  I argued")
P("       that a parity pair adds x to A_2 and x to B_2 so it moves D = A_2 - (3/4) B_2 by x/4,")
P("       while an unpaired multiplet adds x to A_2 alone and moves D by x -- a clean factor of")
P("       four.  That reasoning silently assumes p = +1 means every tower of the multiplet is")
P("       periodic.  It does not: (rep, e, p) fixes the orbifold assignment, and the resulting")
P("       charge table carries BOTH signs of s.  The measured ratios are not 4 and are not even")
P("       of one sign, which is the assumption failing rather than the arithmetic:")
P("")
P("       %-8s %14s %15s %10s %12s" % ("rep", "pair moves D", "the p=+1 moves D", "ratio", "B_2 at p=+1"))
_D0 = moments([])["D"]
_B0 = SD([], 2)[2]
anyneg = False
for rep in ("7", "28", "48", "84"):
    dp = moments([(rep, 1, 1, 1), (rep, 1, -1, 1)])["D"] - _D0
    ds = moments([(rep, 1, 1, 1)])["D"] - _D0
    b = SD([(rep, 1, 1, 1)], 2)[2] - _B0
    anyneg |= ds < 0
    P("       %-8s %14.4f %15.4f %10.4f %12.1f" % (rep, dp, ds, ds / dp, b))
P("")
P("       CONTROL -- a p = +1 multiplet has B_2 != 0, so p is not s : %s"
  % all(abs(SD([(r, 1, 1, 1)], 2)[2] - _B0) > 1e-9 for r in ("7", "28", "48", "84")))
P("       CONTROL -- and one of them moves D the WRONG WAY          : %s" % anyneg)
assert all(abs(SD([(r, 1, 1, 1)], 2)[2] - _B0) > 1e-9 for r in ("7", "28", "48", "84")), \
    "a p = +1 multiplet is purely periodic after all -- then the factor-4 claim was right"
P("       So the paper's answer to \\cite{PGY24} has to stay the vanishing of the ladder's Delta,")
P("       which control (iii) verifies, and must not be dressed up as a ratio of curvatures.")

(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "anchor_delta.json").write_text(json.dumps(dict(
    rows=[{k: (float(v) if isinstance(v, (int, float)) else v) for k, v in r.items()}
          for r in rows], tests=res), indent=1), encoding="utf-8")
P("")
P("archived: outputs/anchor_delta.json")
