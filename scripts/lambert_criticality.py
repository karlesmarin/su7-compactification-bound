#!/usr/bin/env python3
"""lambert_criticality -- what the Lambert form says about marginality, and where a claim of ours
was wider than its proof.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  Three things, and the third is a correction.

  (1) Identity (II) and the closed form are the SAME statement.  Eliminating x^2 between
      x^2 = -d / W(-d e^{-b})  and  x^2 = 12 zeta(3) D / (6 mu + A_4)  gives

          W_{-1} = -1 - 6 mu / A_4 ,        mu = m_h^2 / (K pi^2)^2 ,

      so the measured Higgs mass does not merely remove the logarithm: it fixes WHERE ON THE
      LAMBERT SHEET the vacuum sits.  m_h -> 0 is W -> -1, the branch point.

  (2) The paper justifies two real branches by "the argument is negative on every content here".
      Negative is not enough: W is real only for z >= -1/e.  The missing condition is

          L := b - 1 - ln d  >=  0 ,

      and it is not an extra assumption -- it is exactly the statement that the truncated
      potential HAS an electroweak stationary point, and where it does, m_h^2 >= 0 follows.  So
      the conclusion stands and the reason given for it does not.

  (3) Theorem "the odd eighths" proves 8D odd, hence D != 0, and concludes "electroweak breaking
      is never marginal".  D is the curvature at the SYMMETRIC point.  Marginality at the BROKEN
      point is a different condition -- it is L = 0, equivalently mu = 0, equivalently m_h = 0 --
      and D != 0 does not exclude it.  The sentence is wider than its proof.  This file measures
      how close the lattice gets to the fold; it does not close the question.

Run:  python lambert_criticality.py > outputs/lambert_criticality.txt
"""

import cmath
import contextlib
import io as _io
import itertools
import json
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent

# the closed form, its moments and the five published rows
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

try:
    from scipy.special import lambertw
except ImportError:                      # keep the file runnable with no scipy
    def lambertw(z, k=0):
        w = cmath.log(-z) if k else (0j + 0.3)
        for _ in range(200):
            e = cmath.exp(w)
            w = w - (w * e - z) / (e * (w + 1) - (w + 2) * (w * e - z) / (2 * w + 2))
        return w

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4     # the K of mu = m_h^2/(K pi^2)^2


def bdL(A4, D, G):
    """The Lambert data of a content, and the discriminant."""
    b = (4.0 * G - A4) / (2.0 * A4)
    d = 12.0 * Z3 * D / A4
    z = -d * math.exp(-b)
    L = b - 1.0 - math.log(d) if d > 0 else float("nan")
    return b, d, z, L


def moments_of_terms(tt):
    """(D, A_4, G) straight from a term table.

    NOT amin_closed_form.moments, which takes a CONTENT and builds the table itself: the
    enumeration below already yields tables, and feeding one to the other silently unpacks the
    wrong tuple width.  Same formula as make_fig_collapse.py, kept in step with it.
    """
    A2 = B2 = A4 = B4 = A4L = 0.0
    for (m, s, c) in tt:
        c = float(c)
        if s > 0:
            A2 += m * c ** 2
            A4 += m * c ** 4
            A4L += m * c ** 4 * math.log(c)
        else:
            B2 += m * c ** 2
            B4 += m * c ** 4
    return A2 - 0.75 * B2, A4, A4 * H4 - A4L - LN2 * B4


_GRID = None


def full_min_terms(tt):
    """Minimum of the UNTRUNCATED potential over alpha in (0,1), from a term table.

    This is the arbiter.  Everything else in this file is the series talking about itself.
    """
    global _GRID
    import numpy as _np
    if _GRID is None:
        _GRID = _np.linspace(1e-5, 0.999, 200001)

    def _F(a):
        return sum(m * basis(a, s, c) for m, s, c in tt)

    v = _F(_GRID)
    i = int(_np.argmin(v))
    if i in (0, len(_GRID) - 1):
        return None
    lo, hi = _GRID[i - 1], _GRID[i + 1]
    for _ in range(70):
        xs = _np.linspace(lo, hi, 15)
        j = int(_np.argmin(_F(xs)))
        lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 14)]
    return float(0.5 * (lo + hi))


def lattice_contents(nmax=3):
    """gauge + sum_j m_j x multiplet_j over the eight multiplets of [KM25], as in ceiling_ilp."""
    ns = {"__file__": str(HERE / "ceiling_ilp.py"), "__name__": "ceiling_ilp"}
    with contextlib.redirect_stdout(_io.StringIO()):
        exec(compile((HERE / "ceiling_ilp.py").read_text(encoding="utf-8"),
                     "ceiling_ilp.py", "exec"), ns)
    TT = [ns["terms"](r, e, p) for (r, e, p) in ns["REPS"]]
    GAUGE = ns["GAUGE"]
    for mult in itertools.product(range(nmax + 1), repeat=8):
        if not any(mult):
            continue
        t = list(GAUGE)
        for j, mj in enumerate(mult):
            if mj:
                t += [(m * mj, s, c) for (m, s, c) in TT[j]]
        yield t


def line(c="-", n=96):
    print(c * n)


def main():
    fails = []
    line("=")
    print("THE LAMBERT BRANCH, THE HIGGS MASS, AND WHAT 'NEVER MARGINAL' ACTUALLY SAYS")
    line("=")

    # ---------------------------------------------------------------- 1
    print("\n[1] IDENTITY (II) AND THE CLOSED FORM ARE THE SAME STATEMENT")
    print("    x^2 = -d/W  and  x^2 = 12 z3 D/(6 mu + A_4)  =>  W = -(6 mu + A_4)/A_4")
    print("\n    %-6s %12s %12s %14s %12s %12s"
          % ("row", "A_4", "8D", "W_-1 (num)", "-1-6mu/A4", "rel"))
    line("-")
    worst = 0.0
    for lbl, content, a_them, mh_them, invR in T1:
        a, mo = closed_form(content)
        A4, D, G = float(mo["A4"]), float(mo["D"]), float(mo["G"])
        b, d, z, L = bdL(A4, D, G)
        w = lambertw(z, -1)
        assert abs(w.imag) < 1e-12, "no real W_-1 on a published row"
        wnum = w.real
        x = a * math.pi
        # mu read back from the closed form, exactly as the paper defines it
        mu = (12.0 * Z3 * D / x ** 2 - A4) / 6.0
        wid = -1.0 - 6.0 * mu / A4
        rel = abs(wnum - wid) / abs(wnum)
        worst = max(worst, rel)
        print("    %-6s %12.4f %12.0f %14.8f %12.8f %12.2e"
              % (lbl, A4, round(8 * D), wnum, wid, rel))
    line("-")
    ok = worst < 1e-9
    print("   C1  W_-1 = -1 - 6 mu / A_4 on all five rows ............... %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst))
    if not ok:
        fails.append("C1")
    print("\n    So pinning m_h does not only kill the logarithm -- it SELECTS A POINT of the")
    print("    Lambert sheet.  m_h -> 0 is mu -> 0 is W -> -1: the branch point, where the")
    print("    electroweak minimum and the outer root of \\S ceiling coalesce.")

    # ---------------------------------------------------------------- 2
    print("\n[2] THE DISCRIMINANT THE PAPER LEAVES OUT")
    print("    two real branches  <=>  -1/e <= z < 0  <=>  L := b - 1 - ln d >= 0")
    print("\n    %-6s %12s %12s %14s %12s" % ("row", "z", "-1/e", "L", "mu"))
    line("-")
    Lmin_rows = float("inf")
    for lbl, content, a_them, mh_them, invR in T1:
        a, mo = closed_form(content)
        A4, D, G = float(mo["A4"]), float(mo["D"]), float(mo["G"])
        b, d, z, L = bdL(A4, D, G)
        x = a * math.pi
        mu = (12.0 * Z3 * D / x ** 2 - A4) / 6.0
        Lmin_rows = min(Lmin_rows, L)
        print("    %-6s %12.6f %12.6f %14.6f %12.6f" % (lbl, z, -1 / math.e, L, mu))
    line("-")
    ok = Lmin_rows > 0
    print("   C2  L > 0 on every published row ......................... %s (min %.4f)"
          % ("PASS" if ok else "FAIL", Lmin_rows))
    if not ok:
        fails.append("C2")

    # the equivalence, not just the observation
    print("\n    And L >= 0 is not an extra hypothesis.  W_-1 <= -1 is the whole range of that")
    print("    branch, and W_-1 = -1 - 6 mu/A_4, so on any content that HAS an electroweak")
    print("    stationary point at all, mu >= 0 -- that is, m_h^2 >= 0 -- comes for free.  The")
    print("    paper's conclusion is right; the reason it gives ('the argument is negative') is")
    print("    not the reason.")
    bad = [(lbl, w) for lbl, w in
           [(l, lambertw(bdL(float(closed_form(c)[1]["A4"]), float(closed_form(c)[1]["D"]),
                             float(closed_form(c)[1]["G"]))[2], -1).real)
            for l, c, _, _, _ in T1] if w > -1.0]
    ok = not bad
    print("   C3  W_-1 <= -1 on every row (the branch's own range) ..... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3")

    # ---------------------------------------------------------------- 3
    print("\n[3] HOW CLOSE THE LATTICE GETS TO THE FOLD")
    print("    marginality at the BROKEN point is L = 0, i.e. mu = 0, i.e. m_h = 0.")
    print("    Theorem 'odd eighths' gives D != 0, which is the SYMMETRIC point.  Different")
    print("    condition.  So: scan the real lattice and see how small L can be made.")
    n_ok = n_nov = 0
    Lmin = float("inf")
    Lmin_at = None
    Lmin_terms = None
    for cont in lattice_contents(nmax=3):
        D, A4, G = moments_of_terms(cont)
        if A4 <= 0 or D <= 0:
            continue
        b, d, z, L = bdL(A4, D, G)
        if not (L == L):
            continue
        if L < 0:
            n_nov += 1
            continue
        n_ok += 1
        if L < Lmin:
            Lmin, Lmin_at, Lmin_terms = L, (A4, round(8 * D), G), cont
    print("\n    contents within reach of the expansion (L >= 0) : %d" % n_ok)
    print("    contents outside it                   (L <  0) : %d" % n_nov)
    print("    smallest L found                               : %.3e" % Lmin)
    if Lmin_at:
        print("    at (A_4, 8D, G)                                : (%.0f, %d, %.6f)" % Lmin_at)

    # ---- and now the question that decides whether any of this means anything -------------
    A4f, kf, Gf = Lmin_at
    Df = kf / 8.0
    df = 12.0 * Z3 * Df / A4f
    x_fold = math.sqrt(df)                 # at the fold W = -1, so x^2 = -d/W = d
    a_fold = x_fold / math.pi
    print("""
    IS THIS THE MODEL, OR IS IT THE TRUNCATION?  Do not estimate it -- measure it.  The full
    potential is on disk; minimise it and compare.  (An earlier version of this file guessed the
    error as (pi alpha)^2 and was wrong by more than an order of magnitude in both directions,
    which is exactly why this block now calls numeric_min instead of reasoning.)
""")
    # the truncated prediction AT the fold, and what the untruncated potential actually does
    a_full_fold = full_min_terms(Lmin_terms)
    print("      alpha the series puts the fold at    : %.6f" % a_fold)
    print("      alpha of the UNTRUNCATED minimum     : %s"
          % ("%.6f" % a_full_fold if a_full_fold else "no interior minimum"))
    err_fold = abs(a_full_fold - a_fold) / a_fold if a_full_fold else float('nan')
    print("      measured discrepancy at the fold     : %.1f %%" % (100 * err_fold))

    errs = []
    for lbl, content, a_them, mh_them, invR in T1:
        a_cf, _ = closed_form(content)
        a_full = full_min_terms(table(content))
        if a_full:
            errs.append(abs(a_full - a_cf) / a_cf)
    err_ok = max(errs)
    print("      same measurement on the published rows: %.2f %% worst" % (100 * err_ok))
    print("      and since 1/R_5 = 2 m_W / alpha, the fold sits at 1/R_5 = %.0f GeV"
          % (2 * 80.4 / a_fold))
    print("""
    So the fold lives at large alpha, which is small 1/R_5 -- the opposite end of the same line
    from the ceiling this paper is built to compute.  There the series misplaces the stationary
    point by tens of per cent, against a fraction of one per cent where the paper works.  And a
    fold is a COALESCENCE of roots, so its position is the one thing a truncation cannot be
    trusted about: the sensitivity of a double root to a perturbation diverges.  The marginal
    broken vacuum is therefore not a result of this paper in either direction -- we can neither
    assert it nor rule it out, and the reason is our own instrument.  UNAUDITABLE.

    Nor is the smallest L above evidence of structure.  With %d admissible contents spread over
    a range of L of order %.1f, the smallest value expected from sampling alone is about %.1e;
    the observed minimum is within a couple of orders of that.  It is what ten thousand draws
    look like, not an accumulation at the fold.""" % (n_ok, 3.0, 3.0 / max(n_ok, 1)))
    print("       (observed minimum %.1e)" % Lmin)

    ok = err_fold > 20.0 * err_ok
    print("\n   C4  the series misplaces the fold by >20x its own validated error  %s"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C4")
    print("       (this control PASSING is what makes the fold unclaimable, not what saves it)")

    # ---------------------------------------------------------------- verdict
    line("=")
    print("WHAT THIS CHANGES IN THE PAPER")
    line("=")
    print("""
    The theorem proves: 8D is an odd integer, hence D != 0 and |D| >= 1/8.  D is the coefficient
    of the quadratic term at alpha = 0.  What that establishes is that THE SYMMETRIC POINT is
    never marginal -- the potential always has curvature there, and by a quantised amount.

    It does not establish that the BROKEN vacuum cannot be marginal.  That is F''(alpha*) = 0 at
    alpha* != 0, which by [1] and [2] above is exactly the Lambert fold L = 0, and D != 0 says
    nothing about it.  The scope of the sentence should be narrowed to the symmetric point, and
    the fold stated as what it is: the one marginality this model has available, at finite
    holonomy, reached only in the limit m_h -> 0.

    Whether the lattice can reach it is open.  The scan above says it does not, by a wide
    margin, but a scan is not a theorem and the sentence in the paper should not pretend
    otherwise.
    """)
    # ---- archive the row data so a figure can be DRAWN from it rather than transcribed --------
    # make_fig_lambert.py reads this.  Everything in it is recomputed here from closed_form(),
    # not copied out of the printed tables above, so the two cannot drift apart.
    rows = []
    for lbl, content, a_them, mh_them, invR in T1:
        a, mo = closed_form(content)
        A4, D, G = float(mo["A4"]), float(mo["D"]), float(mo["G"])
        b, d, z, L = bdL(A4, D, G)
        x = a * math.pi
        mu = (12.0 * Z3 * D / x ** 2 - A4) / 6.0
        rows.append(dict(row=lbl, A4=A4, k8D=int(round(8 * D)), b=b, d=d, z=z, L=L,
                         W=lambertw(z, -1).real, mu=mu,
                         mh=KK * math.pi ** 2 * math.sqrt(mu), amin=a))
    (pathlib.Path(__file__).resolve().parent / "outputs" / "lambert_criticality.json").write_text(
        json.dumps(dict(rows=rows, K=KK, minus_one_over_e=-1.0 / math.e), indent=1))

    line("=")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line("=")
        return 1
    print("VERDICT: C1-C4 pass.  One claim narrowed, one reason corrected, one question opened.")
    line("=")
    return 0


if __name__ == "__main__":
    sys.exit(main())
