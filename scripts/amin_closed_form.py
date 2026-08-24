#!/usr/bin/env python3
"""alpha_min in closed form, run on the five published rows against OUR F.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

This tests the FORMULA, not the anchor.  The anchor question is "our F against their
published column"; that is untouched here and stays open.  What is measured here is
"our closed form against our own numerical minimisation of the same F" -- so a
disagreement would be a defect of the expansion and nothing else.

The expansion (Wood-Robinson, see ../part_vi/BRANCH_POINT.md).  With x = pi*alpha and
F = sum_terms m * Re Li_5(s e^{i c pi alpha}) -- which is exactly what
../part_vi/su7_anchor_mh.py's basis() computes:

    Re Li_5( e^{i c pi a}) = z(5) - z(3)(c pi a)^2/2 + (c pi a)^4 [H4 - ln(c pi a)]/24 + ...
    Re Li_5(-e^{i c pi a}) = -eta(5) + eta(3)(c pi a)^2/2 - ln2 (c pi a)^4/24 + ...

so, writing  A2 = sum_{s=+1} m c^2,  B2 = sum_{s=-1} m c^2,  A4 = sum_{s=+1} m c^4,
B4 = sum_{s=-1} m c^4,  A4L = sum_{s=+1} m c^4 ln c,

    D = A2 - (3/4) B2                 (3/4 = eta(3)/z(3), their eq. (67)'s 1/n^5)
    G = A4 H4 - A4L - ln2 B4
    F(a) = const - (z(3) D / 2) x^2 + (x^4/24) [ G - A4 ln x ] + O(x^6)

and dF/dx = 0 gives the fixed point

    x^2 = 24 z(3) D / [ 4 G - A4 (4 ln x + 1) ] .

terms(), GAUGE and T1 are EXTRACTED FROM ../part_vi/su7_anchor_mh.py rather than
retyped, so a change there cannot silently desynchronise this file.
"""
import json
import math
import os
import pathlib
import re

import numpy as np

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE.parent / "part_vi" / "su7_anchor_mh.py"
OUT = HERE / "outputs"

# ---------------------------------------------------------------- the source of truth
_src = SRC.read_text(encoding="utf-8")
_ns = {}
for pat, what in ((r"\ndef terms\(.*?\n(?=\n\ndef )", "terms"),
                  (r"\nGAUGE = \[.*?\]\s*#[^\n]*\n", "GAUGE"),
                  (r"\nT1 = \[.*?\)\]\n", "T1")):
    m = re.search(pat, _src, re.S)
    if not m:
        raise SystemExit("FATAL: could not extract %s from %s" % (what, SRC))
    exec(m.group(0), {}, _ns)
terms, GAUGE, T1 = _ns["terms"], _ns["GAUGE"], _ns["T1"]
P("extracted from %s: terms(), GAUGE (%d entries), T1 (%d rows)" % (SRC.name, len(GAUGE), len(T1)))

# ---------------------------------------------------------------- our F, verbatim
# THE NAMES ARE LONG ON PURPOSE.  This file is exec'd INTO the caller's namespace, so every
# global it defines is a global of the caller.  They used to be _n, _nf, _w, _sgn, and on
# 2026-08-22 vacuum_constraint.py wrote  _w = Wof(...)  in a control loop: from that line on,
# basis() multiplied by a Fraction instead of by n^-5 and every F it returned was garbage --
# F(0) = 900900 where it should be -18.76 -- with no exception anywhere.  The archived output
# predated the collision, so the gates that read the archive saw nothing.  The guard below is
# the second half of the repair: a clobber now dies on the next call instead of printing
# plausible numbers.  Ver [[my-fix-introduces-the-next-bug]], [[edits-must-fail-loudly]].
NMAX = 600
_AMIN_N = np.arange(1, NMAX + 1)
_AMIN_NF = _AMIN_N.astype(float)
_AMIN_W = _AMIN_NF ** -5.0
_AMIN_SGN = {1: np.ones(NMAX), -1: (-1.0) ** _AMIN_N}
# The guard's OWN copy of the length.  The first version compared against NMAX, and NMAX is a
# public name in this namespace that a caller may reasonably want for its own purpose -- and two
# of them do: ceiling_ilp.py uses NMAX for the largest multiplet count it enumerates (14) and
# rigorous_ceiling2.py for its winding cutoff (4000).  Both then died on a guard that was
# reporting a clobber of _AMIN_W that had not happened.  A guard against caller collisions must
# not itself stand on a name the caller can collide with.  [[the-guard-itself-can-be-the-liar]]
_AMIN_NMAX = NMAX


def basis(alpha, s, c):
    if not (isinstance(_AMIN_W, np.ndarray) and _AMIN_W.shape == (_AMIN_NMAX,)
            and isinstance(_AMIN_NF, np.ndarray) and _AMIN_NF.shape == (_AMIN_NMAX,)):
        raise RuntimeError("amin_closed_form: _AMIN_W/_AMIN_NF were overwritten by the caller "
                           "-- F would silently return garbage; rename the caller's variable")
    a = np.atleast_1d(np.asarray(alpha, dtype=float))
    return (np.cos(np.outer(a, c * math.pi * _AMIN_NF)) * (_AMIN_W * _AMIN_SGN[s])).sum(axis=1)


def table(content):
    """the full (m, s, c) term table of a row: gauge sector plus its content."""
    out = list(GAUGE)
    for rep, eta, etap, mult in content:
        for m, s, c in terms(rep, eta, etap):
            out.append((m * mult, s, c))
    return out


def F(content, alpha):
    return sum(m * basis(alpha, s, c) for m, s, c in table(content))


GRID = np.linspace(1e-6, 1.0, 200001)


def numeric_min(content):
    v = F(content, GRID)
    i = int(np.argmin(v))
    if i in (0, len(GRID) - 1):
        return None
    lo, hi = GRID[i - 1], GRID[i + 1]
    for _ in range(60):
        xs = np.linspace(lo, hi, 15)
        j = int(np.argmin(F(content, xs)))
        lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, 14)]
    return float(0.5 * (lo + hi))


# ---------------------------------------------------------------- the closed form
Z3 = 1.2020569031595942854
H4 = 1.0 + 1 / 2 + 1 / 3 + 1 / 4
LN2 = math.log(2.0)


def moments(content):
    A2 = B2 = A4 = B4 = A4L = 0.0
    for m, s, c in table(content):
        c = float(c)
        if s > 0:
            A2 += m * c ** 2
            A4 += m * c ** 4
            A4L += m * c ** 4 * math.log(c)
        else:
            B2 += m * c ** 2
            B4 += m * c ** 4
    return dict(A2=A2, B2=B2, A4=A4, B4=B4, A4L=A4L,
                D=A2 - 0.75 * B2, G=A4 * H4 - A4L - LN2 * B4)


def closed_form(content):
    mo = moments(content)
    D, A4, G = mo["D"], mo["A4"], mo["G"]
    if D <= 0:
        return None, mo
    x = 0.1
    for _ in range(500):
        den = 4 * G - A4 * (4 * math.log(x) + 1)
        if den <= 0:
            return None, mo
        xn = math.sqrt(24 * Z3 * D / den)
        if abs(xn - x) < 1e-14:
            x = xn
            break
        x = xn
    return x / math.pi, mo


# ---------------------------------------------------------------- run
P("")
P("=" * 92)
P("THE FIVE PUBLISHED ROWS -- closed form against OUR OWN numerical minimum of the same F")
P("=" * 92)
P("%-5s %8s %8s %11s %11s %8s %11s %7s" %
  ("row", "D", "A4", "ours(num)", "closed form", "err %", "theirs", "ours/th"))
rec = []
for label, content, a_them, mh_them, invR in T1:
    a_num = numeric_min(content)
    a_cf, mo = closed_form(content)
    has48 = any(rep == "48" for rep, _, _, _ in content)
    err = 100 * (a_cf - a_num) / a_num if (a_cf and a_num) else float("nan")
    ratio = a_num / a_them if a_num else float("nan")
    P("%-5s %8.4f %8.2f %11.6f %11.6f %8.3f %11.4f %7.3f%s" %
      (label, mo["D"], mo["A4"], a_num, a_cf, err, a_them, ratio, "   <- has a 48" if has48 else ""))
    rec.append(dict(row=label, D=mo["D"], A4=mo["A4"], G=mo["G"], alpha_ours=a_num,
                    alpha_closed=a_cf, err_pct=err, alpha_theirs=a_them,
                    ratio_ours_theirs=ratio, has_48=has48))

P("")
errs = [abs(r["err_pct"]) for r in rec if r["err_pct"] == r["err_pct"]]
P("closed form vs our own minimisation:  max |err| = %.3f %%   median = %.3f %%" %
  (max(errs), sorted(errs)[len(errs) // 2]))

P("")
P("CONTROL -- the published ratio pattern must come out (paper section 7: 1.94 with a 48, 1.20 without):")
for tag in (True, False):
    rs = [r["ratio_ours_theirs"] for r in rec if r["has_48"] is tag]
    P("   %-12s n=%d   mean ratio ours/theirs = %.3f   [%s]" %
      ("with a 48" if tag else "without", len(rs), sum(rs) / len(rs),
       ", ".join("%.2f" % v for v in rs)))

P("")
P("CONTROL THAT MUST FAIL -- drop the log term (use G alone, x^2 = 24 z3 D / 4G):")
for label, content, a_them, _, _ in T1:
    mo = moments(content)
    if mo["D"] <= 0 or mo["G"] <= 0:
        P("   %-5s n/a" % label)
        continue
    x = math.sqrt(24 * Z3 * mo["D"] / (4 * mo["G"]))
    a_num = next(r["alpha_ours"] for r in rec if r["row"] == label)
    P("   %-5s no-log = %.6f   vs ours %.6f   err %+8.2f %%" %
      (label, x / math.pi, a_num, 100 * (x / math.pi - a_num) / a_num))

OUT.mkdir(exist_ok=True)
(OUT / "amin_closed_form.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "amin_closed_form.json"))
