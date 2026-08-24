#!/usr/bin/env python3
"""higgs_window.py -- the ceiling is quoted at the top of a window, and the top is not measured.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE QUESTION.  Every ceiling in this paper is evaluated at m_h = 127 GeV, the top of \\cite{KM25}'s
window [125, 127].  That is the right choice for an UPPER BOUND, because the ceiling increases with
m_h -- section 1 measures that rather than assuming it.  But 127 is not a value the Higgs mass can
have: the measurement is 125.20 +- 0.11 GeV, so the top of the window sits sixteen standard
deviations away.  The bound is therefore conservative, which is safe, and the PHYSICAL number --
the ceiling at the mass the Higgs actually has -- is a different one and is not quoted anywhere.

WHAT COMES OUT.  At the measured mass the vertex moves down a whole lattice step, from A_4 = 104 to
101, and the ceiling from 9218 to 9087 GeV.  The one-sigma band does NOT move it -- 125.09 and
125.31 both give 101 -- so the sharper statement is stable against the measurement error and
unstable only against the width of a window nobody has to use.

Both numbers are worth having and they say different things:
    1/R5 <= 9.22 TeV   for any content whose Higgs mass lies anywhere in [125, 127]
    1/R5 <= 9.09 TeV   for any content whose Higgs mass is the measured one
and the second is the one a reader wanting a number about the world should take.

Run:  python higgs_window.py
"""
import json
import math
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_PDG, MH_ERR = 125.20, 0.11            # PDG 2024
MH_LO, MH_HI = 125.0, 127.0              # the window of \cite{KM25}

REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])
LN2, LN3 = math.log(2), math.log(3)
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
WV = [1, -1, 3, -3, 6, -6, 9, -9]
U_OFF, W_OFF = -19.5, -1.5


def best_at(Atgt, Ktgt, need_W):
    """the largest logarithmic budget at (A_4, 8D), optionally with W > 0.  Exhaustive."""
    best = [None]
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, u, v, w):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            if need_W and w + m * WV[0] + W_OFF <= 0:
                return
            c = (u + m * GQ_SYM[0][2] + U_OFF) * LN2 + v * LN3
            if best[0] is None or c > best[0]:
                best[0] = c
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            rec(i + 1, a - c * AV[j], k - c * KV[j], u + c * GQ_SYM[j][2],
                v + c * GQ_SYM[j][1], w + c * WV[j])

    rec(0, Atgt, Ktgt, 0, 0, 0)
    return best[0]


def x_of(t, mh):
    return math.sqrt(12 * Z3 * (1 / 8.0) / (6 * MU(mh) + t))


def vertex(mh, need_W=True, top=280):
    """largest A_4 on 8D = 1 whose budget covers what (I) demands."""
    mu = MU(mh)
    for t in range(top, 59, -1):
        if (t + 1) % 3:
            continue
        x = x_of(t, mh)
        need = float(25 * t) / 12 - (t * (math.log(x) + 0.75) + 3 * mu)
        have = best_at(t - A0, 1 - K0, need_W)
        if have is not None and have >= need:
            return t, 2 * math.pi * MW / x
    return None, None


# ================================================================= 1
P("=" * 100)
P("1 -- THE CEILING RISES WITH m_h, SO THE TOP OF THE WINDOW IS THE CONSERVATIVE CHOICE")
P("=" * 100)
P("  mu ~ m_h^2, and (II) is x^2 = 12 zeta(3) D / (A_4 + 6 mu), so a larger m_h means a larger")
P("  denominator, a smaller x, and 1/R5 = 2 pi m_W / x LARGER.  Measured rather than asserted,")
P("  at the fixed vertex A_4 = 104 so that only m_h moves:")
P("")
P("  %10s %12s %16s" % ("m_h", "mu", "1/R5 at 104"))
prev, mono = None, True
for mh in (124.0, 125.0, 125.20, 126.0, 127.0, 128.0):
    r = 2 * math.pi * MW / x_of(104, mh)
    if prev is not None and r <= prev:
        mono = False
    prev = r
    P("  %10.2f %12.4f %16.1f" % (mh, MU(mh), r))
P("")
P("  CONTROL -- strictly increasing in m_h : %s" % mono)
assert mono, "the ceiling is not monotone in m_h -- the top of the window is not the safe end"
P("  So quoting the ceiling at 127 is right for an upper bound over the window.  What follows is")
P("  not a correction to that; it is the other number, which the paper does not give.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- BUT 127 IS NOT A VALUE THE HIGGS MASS CAN HAVE")
P("=" * 100)
P("  measurement (PDG): m_h = %.2f +- %.2f GeV" % (MH_PDG, MH_ERR))
P("  the window used  : [%.0f, %.0f], which is \\cite{KM25}'s" % (MH_LO, MH_HI))
P("")
for lab, v in (("the window's top, 127", MH_HI), ("the window's floor, 125", MH_LO)):
    P("     %-26s is %5.1f sigma from the measurement" % (lab, abs(v - MH_PDG) / MH_ERR))
P("")
P("  %-30s %10s %12s %14s" % ("m_h", "value", "A_4 vertex", "1/R5 (GeV)"))
rows = []
for lab, mh in (("PDG, -1 sigma", MH_PDG - MH_ERR), ("PDG central", MH_PDG),
                ("PDG, +1 sigma", MH_PDG + MH_ERR), ("PDG, +5 sigma", MH_PDG + 5 * MH_ERR),
                ("the window's floor, 125", MH_LO), ("the window's top, 127", MH_HI)):
    t, r = vertex(mh)
    rows.append((lab, mh, t, r))
    P("  %-30s %10.2f %12s %14.1f" % (lab, mh, t, r))
P("")
band = {t for lab, mh, t, r in rows[:3]}
P("  CONTROL -- the one-sigma band moves the vertex : %s" % (len(band) > 1))
P("             (it must not, or the sharper number is not stable against the measurement)")
assert len(band) == 1, "the measurement error alone moves the vertex -- the sharper bound is unstable"
P("")
P("  the vertex at the measured mass : A_4 = %d" % rows[1][2])
P("  the vertex at the window's top  : A_4 = %d" % rows[5][2])
P("  they differ by %d, which is a whole step of the mod-6 lattice." % (rows[5][2] - rows[1][2]))

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- THE TWO NUMBERS, AND WHAT EACH ONE MEANS")
P("=" * 100)
P("  Both are true and they quantify different statements:")
P("")
P("     1/R5 <= %.2f TeV   for any bulk content whose Higgs mass lies ANYWHERE in [125, 127]"
  % (rows[5][3] / 1000))
P("     1/R5 <= %.2f TeV   for any bulk content whose Higgs mass is the MEASURED one"
  % (max(r for _, _, _, r in rows[:3]) / 1000))
P("")
P("  a %.1f %% difference, and the second is the one a reader wanting a number about the world"
  % (100 * (rows[5][3] - rows[1][3]) / rows[5][3]))
P("  should take.")
P("")
P("  And the unconstrained ceiling of the same section, for completeness -- no W > 0 imposed:")
P("")
P("  %-30s %10s %12s %14s" % ("m_h", "value", "A_4 vertex", "1/R5 (GeV)"))
free = []
for lab, mh in (("PDG central", MH_PDG), ("the window's top, 127", MH_HI)):
    t, r = vertex(mh, need_W=False)
    free.append((lab, mh, t, r))
    P("  %-30s %10.2f %12s %14.1f" % (lab, mh, t, r))

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- CONTROLS")
P("=" * 100)
P("  (i) THE SEARCH MUST BE ABLE TO RETURN SOMETHING ELSE.  Push m_h well outside the window in")
P("      both directions and the vertex has to move; if it does not, the scan is not measuring")
P("      m_h at all:")
P("")
P("      %10s %12s %14s" % ("m_h", "A_4 vertex", "1/R5 (GeV)"))
seen = set()
for mh in (110.0, 120.0, 125.20, 130.0, 140.0):
    t, r = vertex(mh)
    seen.add(t)
    P("      %10.2f %12s %14.1f" % (mh, t, r))
P("      distinct vertices over that range : %d   <-- must be more than one" % len(seen))
assert len(seen) > 1, "the vertex does not respond to m_h -- the scan is inert"

P("")
P("  (ii) AND THE PAPER'S OWN NUMBER MUST COME BACK OUT.  At the top of the window the vertex")
P("       and the ceiling must be exactly what section ceiling prints:")
P("       A_4 = %s (paper: 104)   1/R5 = %.0f GeV (paper: 9218)" % (rows[5][2], rows[5][3]))
assert rows[5][2] == 104 and abs(rows[5][3] - 9218) < 1.0, \
    "the reconstruction does not reproduce the paper's vertex -- something else has drifted"

P("")
P("  (iii) THE NUMBERS AS THE PAPER WOULD PRINT THEM:")
P("       ceiling at the measured Higgs mass : %.0f GeV = %.2f TeV"
  % (rows[1][3], rows[1][3] / 1000))
P("       its vertex                         : A_4 = %d" % rows[1][2])
P("       at +1 sigma                        : %.0f GeV" % rows[2][3])
P("       at the window's top                : %.0f GeV = %.2f TeV" % (rows[5][3], rows[5][3] / 1000))
P("       the window's top, in sigma         : %.1f" % ((MH_HI - MH_PDG) / MH_ERR))
P("       unconstrained, measured mass       : %.0f GeV = %.2f TeV" % (free[0][3], free[0][3] / 1000))
P("       unconstrained, window's top        : %.0f GeV = %.2f TeV" % (free[1][3], free[1][3] / 1000))

(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "higgs_window.json").write_text(json.dumps(dict(
    mh_pdg=MH_PDG, mh_err=MH_ERR, window=[MH_LO, MH_HI],
    rows=[dict(label=l, mh=m, A4=t, invR=r) for l, m, t, r in rows],
    unconstrained=[dict(label=l, mh=m, A4=t, invR=r) for l, m, t, r in free]),
    indent=1), encoding="utf-8")
P("")
P("archived: outputs/higgs_window.json")
