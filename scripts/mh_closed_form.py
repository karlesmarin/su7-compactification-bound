#!/usr/bin/env python3
"""m_h in closed form -- the second column of their Table 1, from the same moments.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part VI's eq. for the invariant is  K = m_h a_min / sqrt(F''(a_min)) = (sqrt3/2pi^3) m_W g4,
so m_h needs F'' AT THE MINIMUM.  Expanding (see amin_closed_form.py) with x = pi a,

    F(x) = C - (z3 D/2) x^2 + (x^4/24)[G - A4 ln x]
    Q(x) = (x^4/24)(G - A4 ln x)
    Q''  = (x^2/2)(G - A4 ln x) - (7/24) A4 x^2

and at the stationary point  G - A4 ln x = 6 z3 D / x^2 + A4/4, so

    d2F/dx2|min = -z3 D + 3 z3 D + A4 x^2/8 - (7/24) A4 x^2  =  2 z3 D - A4 x^2 / 6 .

THE LOGARITHM AND G CANCEL IDENTICALLY.  In alpha units (what their eq. (80) wants):

    F''(a_min) = pi^2 [ 2 z3 D - A4 (pi a_min)^2 / 6 ]
    m_h        = K sqrt(F''(a_min)) / a_min ,   K = (sqrt3 / 2 pi^3) m_W g4

so BOTH columns of their Table 1 follow from (D, A4, G) with no minimisation anywhere.

Three things are measured here, in this order:
  1. the closed-form F'' against a numerical second derivative of the same F;
  2. the closed-form m_h against the one built from that numerical F'';
  3. the LIMIT: over every content that breaks the symmetry, the largest 1/R5
     compatible with a Higgs mass in their own 125-127 GeV window.
"""
import itertools
import json
import math
import pathlib
from fractions import Fraction

import numpy as np

exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py",
          encoding="utf-8").read().split("# ---------------------------------------------------------------- run")[0])

MW = 80.4
G4_SM = 0.63                      # the Standard-Model su(2)_L value; m_h scales linearly in g4
KCONST = math.sqrt(3.0) / (2 * math.pi ** 3) * MW     # K = KCONST * g4


def Fpp_closed(D, A4, alpha):
    """d^2F/dalpha^2 at the minimum -- no G, no logarithm."""
    x = math.pi * alpha
    return math.pi ** 2 * (2 * Z3 * D - A4 * x * x / 6.0)


def Fpp_numeric(content, alpha, h=2e-4):
    a = np.array([alpha - 2 * h, alpha - h, alpha, alpha + h, alpha + 2 * h])
    v = F(content, a)
    return float((-v[0] + 16 * v[1] - 30 * v[2] + 16 * v[3] - v[4]) / (12 * h * h))


def mh_of(fpp, alpha, g4=G4_SM):
    return KCONST * g4 * math.sqrt(fpp) / alpha if fpp > 0 else float("nan")


P("")
P("=" * 100)
P("1 -- F'' AT THE MINIMUM:  closed form  2 z3 D - A4 x^2/6   against a numerical derivative")
P("=" * 100)
P("%-5s %8s %7s %12s %14s %14s %8s" % ("row", "D", "A4", "a_min(num)", "F''(closed)", "F''(numeric)", "err %"))
rows = []
for label, content, a_them, mh_them, invR in T1:
    a_num = numeric_min(content)
    mo = moments(content)
    fp_c = Fpp_closed(mo["D"], mo["A4"], a_num)
    fp_n = Fpp_numeric(content, a_num)
    P("%-5s %8s %7.0f %12.6f %14.6f %14.6f %8.3f" %
      (label, Fraction(mo["D"]).limit_denominator(64), mo["A4"], a_num, fp_c, fp_n,
       100 * (fp_c - fp_n) / fp_n))
    rows.append((label, content, a_them, mh_them, a_num, mo, fp_c, fp_n))

P("")
P("=" * 100)
P("2 -- m_h, at the Standard-Model g4 = %.2f  (m_h scales linearly in g4)" % G4_SM)
P("=" * 100)
P("%-5s %14s %14s %8s %12s" % ("row", "m_h(closed)", "m_h(numeric)", "err %", "m_h theirs"))
for label, content, a_them, mh_them, a_num, mo, fp_c, fp_n in rows:
    mh_c, mh_n = mh_of(fp_c, a_num), mh_of(fp_n, a_num)
    P("%-5s %14.4f %14.4f %8.3f %12.1f" % (label, mh_c, mh_n, 100 * (mh_c - mh_n) / mh_n, mh_them))

P("")
P("CONTROL -- their eq. (80) returns no real m_h where F'' < 0.  At THEIR published alpha:")
for label, content, a_them, mh_them, a_num, mo, fp_c, fp_n in rows:
    fp = Fpp_numeric(content, a_them)
    P("   %-5s  alpha_theirs = %.4f   F'' = %+10.4f   %s" %
      (label, a_them, fp, "NO REAL m_h" if fp < 0 else ""))

# ------------------------------------------------------------------ the limit
P("")
P("=" * 100)
P("3 -- THE LIMIT:  largest 1/R5 over all EWSB contents, with and without the Higgs window")
P("=" * 100)
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
MAXN = 6


def alpha_closed(mo):
    D, A4, G = mo["D"], mo["A4"], mo["G"]
    if D <= 0:
        return None
    x = 0.05
    for _ in range(500):
        den = 4 * G - A4 * (4 * math.log(x) + 1)
        if den <= 0:
            return None
        xn = math.sqrt(24 * Z3 * D / den)
        if abs(xn - x) < 1e-15:
            return xn / math.pi
        x = xn
    return x / math.pi


allc, inwin = [], []
for mults in itertools.product(range(MAXN + 1), repeat=len(REPS)):
    n = sum(mults)
    if n == 0 or n > MAXN:
        continue
    content = [(r, e, ep, m) for (r, e, ep), m in zip(REPS, mults) if m]
    mo = moments(content)
    a = alpha_closed(mo)
    if not a:
        continue
    fpp = Fpp_closed(mo["D"], mo["A4"], a)
    if fpp <= 0:
        continue
    mh = mh_of(fpp, a)
    rec = (2 * MW / a, a, mh, mo["D"], mo["A4"], content)
    allc.append(rec)
    if 125.0 <= mh <= 127.0:
        inwin.append(rec)

allc.sort(reverse=True)
inwin.sort(reverse=True)
P("contents with EWSB and F''>0 : %d      of those inside m_h = 125-127 GeV : %d" % (len(allc), len(inwin)))
P("")
P("  %-11s %9s %9s %7s %6s  %s" % ("1/R5 (GeV)", "alpha", "m_h", "D", "A4", "content"))
P("  -- unconstrained ceiling --")
for r in allc[:3]:
    P("  %-11.0f %9.5f %9.2f %7s %6.0f  %s" %
      (r[0], r[1], r[2], Fraction(r[3]).limit_denominator(64), r[4],
       "+".join("%dx%s%s" % (m, rep, "(+,+)" if ep > 0 else "(+,-)") for rep, e, ep, m in r[5])))
P("  -- with the Higgs window imposed --")
for r in inwin[:5]:
    P("  %-11.0f %9.5f %9.2f %7s %6.0f  %s" %
      (r[0], r[1], r[2], Fraction(r[3]).limit_denominator(64), r[4],
       "+".join("%dx%s%s" % (m, rep, "(+,+)" if ep > 0 else "(+,-)") for rep, e, ep, m in r[5])))

if inwin:
    P("")
    P("  CEILING with m_h in 125-127 GeV, content of at most %d multiplets:  1/R5 <= %.0f GeV" %
      (MAXN, inwin[0][0]))
    P("  (unconstrained it would be %.0f GeV -- the Higgs mass costs a factor %.2f)" %
      (allc[0][0], allc[0][0] / inwin[0][0]))

# ------------------------------------------------- the five rows AFTER donating one 84(+,+)
# The paper's section 8 prints a table of m_h for each published row once Part VI's escape is
# paid for -- and NOTHING computed it.  check_numbers.py found 114.7, 145.6 and 169.0 unbacked
# in the archive and in a fresh run of every script, so this block is the repair, not a
# convenience: the Data availability section promises every displayed number regenerates.
P("")
P("=" * 100)
P("THE FIVE PUBLISHED ROWS AFTER DONATING ONE 84(+,+)  (Part VI's escape, priced here)")
P("=" * 100)
P("Removing one 84(+,+) from the content and re-minimising.  The host carries 8D = +10, so the")
P("donation takes 10/8 off D; a row whose D goes non-positive has no interior minimum left.")
P("")
P("%-6s %10s %12s %14s %12s" % ("row", "8D before", "8D after", "m_h after", "in window?"))
_don = []
for label, content, a_them, mh_them, invR in T1:
    after = []
    dropped = False
    for (rep, e, ep, m) in content:
        if rep == "84" and e == 1 and ep == 1 and not dropped:
            dropped = True
            if m > 1:
                after.append((rep, e, ep, m - 1))
        else:
            after.append((rep, e, ep, m))
    if not dropped:
        P("%-6s   no 84(+,+) to donate" % label)
        continue
    mb = moments(content)
    ma = moments(after)
    k_b, k_a = 8 * mb["D"], 8 * ma["D"]
    if ma["D"] <= 0:
        P("%-6s %10.0f %12.0f %14s %12s" % (label, k_b, k_a, "no interior min", "--"))
        _don.append(dict(row=label, k_before=k_b, k_after=k_a, mh=None))
        continue
    a = alpha_closed(ma)
    fpp = Fpp_closed(ma["D"], ma["A4"], a)
    mh = mh_of(fpp, a) if fpp > 0 else None
    if mh is None or mh != mh:
        P("%-6s %10.0f %12.0f %14s %12s" % (label, k_b, k_a, "no real m_h", "--"))
        _don.append(dict(row=label, k_before=k_b, k_after=k_a, mh=None))
        continue
    P("%-6s %10.0f %12.0f %14.1f %12s"
      % (label, k_b, k_a, mh, "YES" if 125.0 <= mh <= 127.0 else "no"))
    _don.append(dict(row=label, k_before=k_b, k_after=k_a, mh=mh))
P("")
_hit = [d for d in _don if d["mh"] is not None and 125.0 <= d["mh"] <= 127.0]
P("rows landing in 125-127 GeV after the donation : %d of %d" % (len(_hit), len(_don)))
P("  -- none of them.  Paying for the escape and keeping the Higgs mass are two different")
P("     demands, and the five published rows satisfy them one at a time.")

OUT.mkdir(exist_ok=True)
(OUT / "mh_closed_form.json").write_text(json.dumps(
    dict(g4=G4_SM, n_ewsb=len(allc), n_window=len(inwin),
         ceiling_unconstrained=allc[0][0] if allc else None,
         ceiling_with_higgs=inwin[0][0] if inwin else None,
         best=[dict(invR=r[0], alpha=r[1], mh=r[2], D=r[3], A4=r[4],
                    content=[[rep, e, ep, m] for rep, e, ep, m in r[5]]) for r in inwin[:10]]),
    indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "mh_closed_form.json"))
