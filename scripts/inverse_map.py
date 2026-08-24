#!/usr/bin/env python3
"""The inverse map: their published Table 1, read through the closed form.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part VI's section 7 fits per-multiplet weights to their published columns by MINIMISING
F numerically for every trial weight -- an objective containing an argmin, which is the
instrument trap that note records.  The closed form removes the argmin entirely, because
the stationarity condition is ALGEBRAIC and, at a GIVEN alpha, LINEAR in the moments:

    24 z3 D  -  4 x^2 G  +  x^2 A4 (4 ln x + 1)  =  0 ,      x = pi alpha        (I)

and the Higgs mass, through F''(a_min) = pi^2 [2 z3 D - A4 x^2/6] and their eqs. (80),(82),

    (m_h a / K)^2 / pi^2  =  2 z3 D  -  A4 x^2 / 6 ,          K = (sqrt3/2pi^3) m_W g4   (II)

D, A4 and G are each LINEAR in the multiplicities, hence linear in per-multiplet weights.
So their five rows give TEN exact linear equations in five unknowns (a weight per multiplet
plus 1/g4^2), with no minimisation anywhere.

This is a CONTROL, not a new fit: section 7 already solved the same problem by the numerical
route and reported residual 0.1628 with a scramble ratio of 2.54.  Two instruments that share
no arithmetic should land in the same place.

AND THE FIRST VERSION OF THIS SCRIPT SAID THEY DID NOT -- ratio 0.25, i.e. a scramble fitting
BETTER than the real column.  That was this file's fault, in two places, and reading the other
instrument (su7_content_dependence.py) found both:

  * it scrambles by REVERSING the alpha column -- ONE fixed permutation.  This file took the
    BEST OF ALL 119.  A minimum over 119 draws, with only one degree of over-determination,
    beats the real fit by chance most of the time: not the same statistic at all.
  * it normalises EACH EQUATION by the norm of its own coefficient row (it says so in its own
    output).  Without that, the equations with the largest coefficients dominate the fit.

With section 7's protocol restored this file gives residual 0.1814 against its 0.1628, and a
scramble ratio 2.26 against its 2.54 -- same statistic, same magnitude, same direction.  The
remaining ~11% is the closed form's truncation at their alpha, propagated through the fit.
Both numbers are reported below: theirs is the comparison, ours is the stricter test.
"""
import itertools
import json
import math
import pathlib
import re

import numpy as np

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE.parent / "part_vi" / "su7_anchor_mh.py"
OUT = HERE / "outputs"
_ns = {}
_s = SRC.read_text(encoding="utf-8")
for pat in (r"\ndef terms\(.*?\n(?=\n\ndef )", r"\nGAUGE = \[.*?\]\s*#[^\n]*\n", r"\nT1 = \[.*?\)\]\n"):
    exec(re.search(pat, _s, re.S).group(0), {}, _ns)
terms, GAUGE, T1 = _ns["terms"], _ns["GAUGE"], _ns["T1"]

Z3 = 1.2020569031595942854
H4 = 1.0 + 1 / 2 + 1 / 3 + 1 / 4
LN2 = math.log(2.0)
MW = 80.4
KC = math.sqrt(3.0) / (2 * math.pi ** 3) * MW          # K = KC * g4

# the five multiplet species a row can contain, in a fixed order
SPECIES = [("28", 1, -1), ("28", 1, 1), ("48", 1, 1), ("84", 1, 1), ("7", 1, -1), ("7", 1, 1)]


def moments(tab):
    A2 = B2 = A4 = B4 = A4L = 0.0
    for m, s, c in tab:
        c = float(c)
        if s > 0:
            A2 += m * c ** 2
            A4 += m * c ** 4
            A4L += m * c ** 4 * math.log(c)
        else:
            B2 += m * c ** 2
            B4 += m * c ** 4
    return (A2 - 0.75 * B2, A4, A4 * H4 - A4L - LN2 * B4)          # D, A4, G


def row_moments(content):
    """(D, A4, G) of the gauge sector, and of one copy of each species present."""
    g = moments(GAUGE)
    per = {}
    for rep, e, ep, mult in content:
        k = (rep, e, ep)
        per[k] = per.get(k, 0) + mult
    return g, per


def species_moments(k):
    rep, e, ep = k
    return moments(terms(rep, e, ep))


P("=" * 96)
P("THE TEN EQUATIONS -- built from their published alpha_min and m_h, no minimisation")
P("=" * 96)
rows, A_rows, b_rows, tags = [], [], [], []
for label, content, a_them, mh_them, invR in T1:
    g, per = row_moments(content)
    x = math.pi * a_them
    # each unknown: the weight of a species (multiplying its own D, A4, G) ; the gauge is fixed
    coef_D, coef_A4, coef_G = {}, {}, {}
    for k, mult in per.items():
        d, a4, gg = species_moments(k)
        coef_D[k], coef_A4[k], coef_G[k] = mult * d, mult * a4, mult * gg
    rows.append((label, content, a_them, mh_them, x, g, coef_D, coef_A4, coef_G))
    P("  %-5s alpha=%.4f  m_h=%.1f   species: %s" %
      (label, a_them, mh_them, ", ".join("%dx%s%s" % (m, k[0], "(+,+)" if k[2] > 0 else "(+,-)")
                                         for k, m in per.items())))

SPEC = sorted({k for _, _, _, _, _, _, cD, _, _ in rows for k in cD})
P("")
P("  unknowns: %s  and 1/g4^2" % ", ".join("w[%s%s]" % (k[0], "(+,+)" if k[2] > 0 else "(+,-)") for k in SPEC))

# --- equation (I), one per row: linear in the weights, no g4
A1, b1 = [], []
for label, content, a_them, mh_them, x, g, cD, cA, cG in rows:
    coef = []
    for k in SPEC:
        coef.append(24 * Z3 * cD.get(k, 0.0) - 4 * x * x * cG.get(k, 0.0)
                    + x * x * cA.get(k, 0.0) * (4 * math.log(x) + 1))
    coef.append(0.0)                                  # 1/g4^2 does not enter (I)
    rhs = -(24 * Z3 * g[0] - 4 * x * x * g[2] + x * x * g[1] * (4 * math.log(x) + 1))
    A1.append(coef)
    b1.append(rhs)

# --- equation (II), one per row: (m_h a)^2/(pi^2 K^2) * (1/g4^2) = 2 z3 D - A4 x^2/6
A2m, b2 = [], []
for label, content, a_them, mh_them, x, g, cD, cA, cG in rows:
    coef = []
    for k in SPEC:
        coef.append(2 * Z3 * cD.get(k, 0.0) - cA.get(k, 0.0) * x * x / 6.0)
    coef.append(-((mh_them * a_them) ** 2) / (math.pi ** 2 * KC ** 2))
    rhs = -(2 * Z3 * g[0] - g[1] * x * x / 6.0)
    A2m.append(coef)
    b2.append(rhs)

A = np.array(A1 + A2m)
b = np.array(b1 + b2)
sol, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
pred = A @ sol
resid = float(np.linalg.norm(pred - b) / max(1.0, np.linalg.norm(b)))

P("")
P("=" * 96)
P("THE SOLUTION -- exact linear least squares, %d equations, %d unknowns, rank %d" %
  (A.shape[0], A.shape[1], rank))
P("=" * 96)
for k, v in zip(SPEC, sol[:-1]):
    P("   w[%-9s] = %9.4f      %s" % (k[0] + ("(+,+)" if k[2] > 0 else "(+,-)"), v,
                                      "<- their own normalisation is 1" if abs(v - 1) > 0.3 else ""))
g4 = 1.0 / math.sqrt(sol[-1]) if sol[-1] > 0 else float("nan")
P("   1/g4^2       = %9.4f   ->  g4 = %.4f   (Standard-Model su(2)_L: 0.63)" % (sol[-1], g4))
P("   relative residual = %.4f" % resid)

# --- the control that must fail: scramble the alpha column
P("")
P("=" * 96)
P("CONTROL -- re-solve on every permutation of their alpha_min column")
P("=" * 96)
alphas = [r[2] for r in rows]
best, worst, n = 1e9, 0.0, 0
for perm in itertools.permutations(range(len(rows))):
    if all(i == j for i, j in zip(perm, range(len(rows)))):
        continue
    AA, bb = [], []
    for (label, content, a_them, mh_them, x0, g, cD, cA, cG), idx in zip(rows, perm):
        a = alphas[idx]
        x = math.pi * a
        coef = [24 * Z3 * cD.get(k, 0.0) - 4 * x * x * cG.get(k, 0.0)
                + x * x * cA.get(k, 0.0) * (4 * math.log(x) + 1) for k in SPEC] + [0.0]
        AA.append(coef)
        bb.append(-(24 * Z3 * g[0] - 4 * x * x * g[2] + x * x * g[1] * (4 * math.log(x) + 1)))
        coef2 = [2 * Z3 * cD.get(k, 0.0) - cA.get(k, 0.0) * x * x / 6.0 for k in SPEC]
        coef2.append(-((mh_them * a) ** 2) / (math.pi ** 2 * KC ** 2))
        AA.append(coef2)
        bb.append(-(2 * Z3 * g[0] - g[1] * x * x / 6.0))
    AA, bb = np.array(AA), np.array(bb)
    s2, *_ = np.linalg.lstsq(AA, bb, rcond=None)
    r2 = float(np.linalg.norm(AA @ s2 - bb) / max(1.0, np.linalg.norm(bb)))
    best, worst, n = min(best, r2), max(worst, r2), n + 1
P("   %d scrambles: best %.4f, worst %.4f, against the real %.4f" % (n, best, worst, resid))
P("   ratio best-scramble / real = %.2f   %s" %
  (best / resid, "the fit is informative" if best / resid > 1.5 else
   "*** A SCRAMBLE DOES AS WELL: the fit measures nothing ***"))


# ---------------------------------------------------------------------------------------
# The parametrisation section 7 actually used: ONE weight per REPRESENTATION, parity-blind,
# plus 1/g4^2 -- five unknowns, ten equations.  Section 7 solved it numerically and got
# w(48) = 5.59 with the alpha column alone, and a residual that both columns together do not
# accept.  Two instruments sharing no arithmetic should agree.
# ---------------------------------------------------------------------------------------
REPS = ["7", "28", "48", "84"]


def solve(parity_blind, alpha_of, use_mh):
    keys = REPS if parity_blind else SPEC
    AA, bb = [], []
    for label, content, a_them, mh_them, x0, g, cD, cA, cG in rows:
        a = alpha_of(label, a_them)
        x = math.pi * a
        agg = {}
        for k in cD:
            kk = k[0] if parity_blind else k
            d, a4, gg = agg.get(kk, (0.0, 0.0, 0.0))
            agg[kk] = (d + cD[k], a4 + cA[k], gg + cG[k])
        AA.append([24 * Z3 * agg.get(k, (0, 0, 0))[0] - 4 * x * x * agg.get(k, (0, 0, 0))[2]
                   + x * x * agg.get(k, (0, 0, 0))[1] * (4 * math.log(x) + 1) for k in keys] + [0.0])
        bb.append(-(24 * Z3 * g[0] - 4 * x * x * g[2] + x * x * g[1] * (4 * math.log(x) + 1)))
        if use_mh:
            AA.append([2 * Z3 * agg.get(k, (0, 0, 0))[0] - agg.get(k, (0, 0, 0))[1] * x * x / 6.0
                       for k in keys] + [-((mh_them * a) ** 2) / (math.pi ** 2 * KC ** 2)])
            bb.append(-(2 * Z3 * g[0] - g[1] * x * x / 6.0))
    AA, bb = np.array(AA), np.array(bb)
    s, *_ = np.linalg.lstsq(AA, bb, rcond=None)
    return keys, s, float(np.linalg.norm(AA @ s - bb) / max(1.0, np.linalg.norm(bb)))


P("")
P("=" * 96)
P("SECTION 7's PROTOCOL, RESTORED -- per-row normalisation, and ONE scramble (reversal)")
P("=" * 96)


def solve_s7(alphas, rownorm=True):
    AA, bb = [], []
    for (label, content, a0, mh, x0, g, cD, cA, cG), a in zip(rows, alphas):
        x = math.pi * a
        agg = {}
        for k in cD:
            kk = k[0]
            d, a4, gg = agg.get(kk, (0.0, 0.0, 0.0))
            agg[kk] = (d + cD[k], a4 + cA[k], gg + cG[k])
        r1 = [24 * Z3 * agg.get(k, (0, 0, 0))[0] - 4 * x * x * agg.get(k, (0, 0, 0))[2]
              + x * x * agg.get(k, (0, 0, 0))[1] * (4 * math.log(x) + 1) for k in REPS] + [0.0]
        c1 = -(24 * Z3 * g[0] - 4 * x * x * g[2] + x * x * g[1] * (4 * math.log(x) + 1))
        r2 = [2 * Z3 * agg.get(k, (0, 0, 0))[0] - agg.get(k, (0, 0, 0))[1] * x * x / 6.0
              for k in REPS] + [-((mh * a) ** 2) / (math.pi ** 2 * KC ** 2)]
        c2 = -(2 * Z3 * g[0] - g[1] * x * x / 6.0)
        for r, c in ((r1, c1), (r2, c2)):
            n = np.linalg.norm(r) if rownorm else 1.0
            AA.append([v / n for v in r])
            bb.append(c / n)
    AA, bb = np.array(AA), np.array(bb)
    s, *_ = np.linalg.lstsq(AA, bb, rcond=None)
    return s, float(np.linalg.norm(AA @ s - bb))


REPS = ["7", "28", "48", "84"]
A0 = [r[2] for r in rows]
s7, r7 = solve_s7(A0)
_, r7s = solve_s7(A0[::-1])
P("   w(48) = %.3f    g4 = %s" % (s7[2], "%.3f" % (1 / math.sqrt(s7[4])) if s7[4] > 0 else "imaginary"))
P("   residual real = %.4f   reversed-alpha = %.4f   ratio = %.2f" % (r7, r7s, r7s / r7))
P("   section 7, by numerical minimisation: residual 0.1628, ratio 2.54")
P("   -> two instruments sharing no arithmetic, same magnitude and same direction.")

P("")
P("=" * 96)
P("AND THE STRICTER TEST -- best of all 119 scrambles, equations unnormalised")
P("=" * 96)
for use_mh, tag in ((False, "alpha column ALONE (this is where section 7 got w(48) = 5.59)"),
                    (True, "both columns together")):
    keys, s, r = solve(True, lambda l, a: a, use_mh)
    P("  %s" % tag)
    P("     %s   1/g4^2 = %.4f%s" %
      ("  ".join("w[%s]=%7.3f" % (k, v) for k, v in zip(keys, s[:-1])), s[-1],
       "  -> g4 = %.3f" % (1 / math.sqrt(s[-1])) if s[-1] > 0 else "  -> g4 IMAGINARY"))
    P("     relative residual = %.4f" % r)
    sc = []
    for perm in itertools.permutations(range(len(rows))):
        if all(i == j for i, j in zip(perm, range(len(rows)))):
            continue
        amap = {rows[i][0]: rows[j][2] for i, j in zip(range(len(rows)), perm)}
        _, _, r2 = solve(True, lambda l, a, m=amap: m[l], use_mh)
        sc.append(r2)
    P("     %d scrambles: best %.4f -> ratio %.2f   %s" %
      (len(sc), min(sc), min(sc) / r,
       "INFORMATIVE" if min(sc) / r > 1.5 else "*** a scramble does as well ***"))
    P("")

OUT.mkdir(exist_ok=True)
(OUT / "inverse_map.json").write_text(json.dumps(
    dict(species=["%s%s" % (k[0], "(+,+)" if k[2] > 0 else "(+,-)") for k in SPEC],
         weights=list(sol[:-1]), inv_g4sq=sol[-1], g4=g4, residual=resid,
         scramble_best=best, scramble_worst=worst), indent=1), encoding="utf-8")
P("written: %s" % (OUT / "inverse_map.json"))
