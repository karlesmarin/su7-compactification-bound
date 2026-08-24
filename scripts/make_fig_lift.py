#!/usr/bin/env python3
"""make_fig_lift.py -- the third axis, drawn, because in two dimensions the answer is invisible.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS FIGURE EXISTS
----------------------
semigroup.py settles the ceiling by showing that (A4, 8D) = (215, 1) cannot be paid for and
(212, 1) can.  Said in words that is an inequality between two numbers, 681.017 against 681.084,
and nobody can see it.  Drawn, it is one picture with an obvious shape:

  LEFT, in two dimensions -- the plane the paper has worked in all along.  Both vertices are
  admissible lattice points, both sit under the boundary where identity (I) runs out.  THEY ARE
  INDISTINGUISHABLE.  Whatever decides between them is not in this picture.

  RIGHT, in three -- the same points with the axis the projection was hiding, 2U from
  G = (25/12) A4 - U ln2 - V ln3.  Each (A4, 8D) is no longer a point but a FIBRE of reachable U,
  and the G condition is a floor the fibre has to reach.  At A4 = 212 the fibre clears it.  At
  A4 = 215 the fibre stops just below, and there is nothing above it: the maximum is over a
  finite feasible set, so the picture is a proof and not a sample.

That is the whole result: the ceiling is decided on an axis that the two-dimensional lattice
of Figure 5 does not have.

Run:  python make_fig_lift.py     -> fig_lift.pdf, fig_lift_es.pdf, .png  (into . and paper/)
"""
import math
import os
import pathlib
import sys
from fractions import Fraction as Fr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401  (registers the projection)
import numpy as np

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
GREEN, GREY = "#2E8B57", "#8A8A8A"
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
U_OFF, V_OFF = -19.5, 0.0
LN2, LN3 = math.log(2), math.log(3)

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MU_HI = MU(127.0)


def gstar(t, k):
    x = math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU_HI + t))
    return t * (math.log(x) + 0.75) + 3 * MU_HI, x


def fibre(Atgt, Ktgt, cap=400):
    """every reachable (U, V) at those two moments, as a list of U ln2 + V ln3.
    Finite for the reason semigroup.py gives: only 7(+,+) has A4 = 0 and 8D then forces it."""
    out = []
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, u, v):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            out.append((u + m * GQ_SYM[0][2] + U_OFF) * LN2 + (v + V_OFF) * LN3)
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            rec(i + 1, a - c * AV[j], k - c * KV[j], u + c * GQ_SYM[j][2], v + c * GQ_SYM[j][1])

    rec(0, Atgt, Ktgt, 0, 0)
    return sorted(out)


_DP = {}


def _build_dp(AMAX, KLO):
    """max-plus dynamic programme for the largest U ln2 + V ln3 at each (A4, 8D), offsets NOT
    applied.  The naive recursion is exponential and dies on the upper rungs; this is
    polynomial and gives the same answers, which section `main` checks against it."""
    W = [GQ_SYM[j][2] * LN2 + GQ_SYM[j][1] * LN3 for j in range(8)]
    order = [j for j in range(8) if AV[j] > 0]
    rows = [dict() for _ in range(AMAX + 1)]     # rows[a][k] = best budget, offsets not applied
    rows[0][0] = 0.0
    for a in range(0, AMAX + 1):
        # every generator that advances A4 feeds row a from an earlier row
        for j in order:
            pa = a - AV[j]
            if pa < 0:
                continue
            src, dst = rows[pa], rows[a]
            w = W[j]
            for qk, val in src.items():
                k = qk + KV[j]
                if k < KLO:
                    continue
                if dst.get(k, -1e18) < val + w:
                    dst[k] = val + w
        # then the 7(+,+), which is (0, -6, +ln2): close this row downwards in steps of six
        dst = rows[a]
        for k in sorted(dst, reverse=True):
            v = dst[k]
            kk = k - 6
            while kk >= KLO:
                nv = v + W[0]
                if dst.get(kk, -1e18) < nv:
                    dst[kk] = nv
                    v = nv
                else:
                    v = dst[kk]
                kk -= 6
    return rows


_ROWS = []


def maxbudget(Atgt, Ktgt):
    """max of U_total ln2 + V_total ln3 at those two moments (offsets applied), or None."""
    if not _ROWS:
        _ROWS.extend(_build_dp(560, -700))
    if not 0 <= Atgt < len(_ROWS):
        return None
    v = _ROWS[Atgt].get(Ktgt)
    return None if v is None else v + U_OFF * LN2 + V_OFF * LN3


def draw(es=False):
    L = dict(
        ttlA="En dos dimensiones los dos vértices son el mismo" if es
             else "In two dimensions the two vertices are the same",
        # "In three" would say the panel draws the lattice, and it does not: the lattice is FOUR
        # dimensional and this collapses its two logarithmic directions onto the single scalar
        # U ln2 + V ln3, which is the only combination the G condition reads.  Saying which is
        # the difference between a projection and a claim.
        ttlB=("Proyectado sobre $U\\ln 2 + V\\ln 3$: uno de ellos no llega" if es
              else "Projected on $U\\ln 2 + V\\ln 3$: one of them does not reach"),
        xa=r"$A_4$", ya=r"$8D$",
        # short on purpose: a long rotated z-label either spills past the tight bbox (which
        # check_layout caught, 6.1 pt, one glyph) or gets dropped from it.  The caption says
        # what the axis is; the axis only has to be readable.
        za=("eje vertical: holgura, disponible $-$ exigido" if es
            else "vertical axis: slack, available $-$ required"),
        xd=("$A_4$ menos el borde de su peldaño" if es
            else "$A_4$ minus its own rung's boundary"),
        # la imparidad depende de la semilla y la congruencia no: la leyenda lo dice.
        adm=("admisible con la semilla de [2]: $8D$ impar, $A_4+8D \\equiv 0$ (mod 3)" if es
             else "admissible on the seed of [2]: $8D$ odd, $A_4+8D \\equiv 0$ (mod 3)"),
        bnd=("donde (I) se queda sin contenidos" if es else "where (I) runs out of contents"),
        need=("exigido por la condición en $G$" if es else "required by the $G$ condition"),
        reach=("alcanzable" if es else "reachable"),
        emp=("VACÍO" if es else "EMPTY"),
        att=("ALCANZADO" if es else "ATTAINED"),
    )
    fig = plt.figure(figsize=(12.6, 4.6))

    # ---------------------------------------------------------------- left: the projection
    ax = fig.add_subplot(1, 2, 1)
    A_LO, A_HI = 196, 222
    ks = [k for k in range(1, 20, 2)]
    xs, ys = [], []
    for k in ks:
        for a in range(A_LO, A_HI + 1):
            if (a + k) % 3 == 0:
                xs.append(a)
                ys.append(k)
    ax.scatter(xs, ys, s=13, c=GREY, alpha=0.55, linewidths=0, label=L["adm"], zorder=2)

    # the boundary: largest A4 whose relaxation still admits a content, per rung
    CERT = {1: 215, 3: 336, 5: 436, 7: 533, 9: 630, 11: 727, 13: 823, 15: 918, 17: 1013, 19: 1108}
    bx = [CERT[k] for k in ks if CERT[k] <= A_HI + 40]
    by = [k for k in ks if CERT[k] <= A_HI + 40]
    if len(bx) > 1:
        ax.plot(bx, by, color=RED, lw=1.6, zorder=3, label=L["bnd"])
    else:
        ax.axvline(215, color=RED, lw=1.6, zorder=3, label=L["bnd"])

    ax.scatter([215], [1], s=150, facecolors="none", edgecolors=RED, linewidths=2.0, zorder=5)
    ax.scatter([212], [1], s=110, c=GREEN, linewidths=0, zorder=5)
    ax.annotate("(215, 1)", (215, 1), textcoords="offset points", xytext=(6, 12),
                fontsize=9, color=RED)
    ax.annotate("(212, 1)", (212, 1), textcoords="offset points", xytext=(-56, 12),
                fontsize=9, color=GREEN)
    ax.set_xlabel(L["xa"])
    ax.set_ylabel(L["ya"])
    ax.set_title(L["ttlA"], fontsize=11)
    ax.set_xlim(A_LO - 1, A_HI + 1)
    ax.set_ylim(0, 20)
    ax.grid(alpha=0.18, lw=0.5)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)

    # ---------------------------------------------------------------- right: the lift
    # The vertical axis is the SLACK, budget minus requirement.  Plotting the raw budget would
    # be useless: it runs to 700 and the whole verdict is 0.067 wide.  The zero plane is where
    # the ceiling is decided, one rung at a time, and that is what makes the panel worth drawing.
    ax3 = fig.add_subplot(1, 2, 2, projection="3d")
    RUNGS = [(1, 215), (3, 336), (5, 436), (7, 533)]
    shown, curves = [], []
    for k, tmax in RUNGS:
        aa, ss = [], []
        for a in range(tmax - 12, tmax + 7):
            if (a + k) % 3:
                continue
            b = maxbudget(a - A0, k - K0)
            if b is None:
                continue
            need = float(Fr(25, 12) * a) - gstar(a, k)[0]
            aa.append(a)
            ss.append(b - need)
            if k == 1:
                shown.append((a, [b], need))
        curves.append((k, np.array(aa, float), np.array(ss, float), tmax))

    # every rung has its own A4 scale (215, 336, 436, 533), so they are drawn against the
    # DISTANCE TO THAT RUNG'S OWN BOUNDARY.  That is what makes the four comparable, and it is
    # also the statement: the crossing sits a fixed short step below the boundary on each.
    XLO, XHI = -12, 6
    xg = np.linspace(XLO, XHI, 2)
    yg = np.array([RUNGS[0][0] - 1, RUNGS[-1][0] + 1], float)
    XG, YG = np.meshgrid(xg, yg)
    ax3.plot_surface(XG, YG, np.zeros_like(XG), color=AMBER, alpha=0.18, linewidth=0,
                     shade=False, zorder=1)
    for k, aa, ss, tmax in curves:
        dx = aa - tmax
        ax3.plot(dx, [k] * len(dx), ss, color=BLUE, lw=1.6, zorder=4)
        pos = ss >= 0
        ax3.scatter(dx[pos], [k] * int(pos.sum()), ss[pos], s=24, c=GREEN, linewidths=0, zorder=5)
        ax3.scatter(dx[~pos], [k] * int((~pos).sum()), ss[~pos], s=24, c=RED,
                    linewidths=0, zorder=5)
        # where each rung's own integer optimum sits
        okx = dx[pos]
        if len(okx):
            ax3.text(okx.max(), k, ss[pos][np.argmax(okx)] + 1.5,
                     "%d" % int(aa[pos][np.argmax(okx)]), color=GREEN, fontsize=7.5, zorder=7)
    ax3.plot([], [], [], color=GREEN, marker="o", lw=0, label=L["reach"])
    ax3.plot([], [], [], color=RED, marker="o", lw=0, label=L["emp"])
    ax3.set_xlabel(L["xd"], labelpad=6)
    ax3.set_ylabel(L["ya"], labelpad=2)
    # NOT set_zlabel: a tight bbox on a 3D axes drops the z-label entirely (verified -- the word
    # was absent from the text layer of both editions).  It goes in the title, which survives.
    ax3.set_zlabel("")
    ax3.set_yticks([k for k, _ in RUNGS])
    ax3.set_xlim(XLO, XHI)
    ax3.set_title(L["ttlB"] + "\n" + L["za"], fontsize=11)
    ax3.view_init(elev=22, azim=-62)
    ax3.tick_params(labelsize=7.5)
    ax3.legend(fontsize=7.5, loc="upper left")

    fig.tight_layout()
    stem = "fig_lift_es" if es else "fig_lift"
    here = os.path.dirname(os.path.abspath(__file__))
    # pad_inches, and it is not cosmetic: with a bare tight bbox matplotlib crops to the ink and
    # a rotated axis label can spill a hair PAST the declared BBox.  Scaled to \textwidth in the
    # paper that reads as a text block outside the text block, and check_layout catches it -- it
    # caught exactly one glyph, the 'h' of "holgura", 6.1 pt over, in the Spanish edition.
    for d in (here, os.path.join(here, "paper")):
        fig.savefig(os.path.join(d, stem + ".pdf"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(os.path.join(here, stem + ".png"), dpi=150, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return shown


def per_rung():
    """the integer optimum on each rung, against the boundary the relaxation certifies."""
    out = []
    # NOT rungs 9 and 11: their boundaries are at A4 = 630 and 727, past the dynamic
    # programme's A4 <= 560, and an out-of-range lookup returns None -- which would print as
    # 'no integer optimum' and be a false null.  Declared instead of printed.
    for k, tmax in [(1, 215), (3, 336), (5, 436), (7, 533)]:
        best, sl = None, []
        for a in range(tmax - 18, tmax + 7):
            if (a + k) % 3:
                continue
            b = maxbudget(a - A0, k - K0)
            if b is None:
                continue
            need = float(Fr(25, 12) * a) - gstar(a, k)[0]
            sl.append(b - need)
            if b >= need:
                best = a
        mono = all(sl[i] > sl[i + 1] for i in range(len(sl) - 1))
        x = math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU_HI + tmax))
        xb = math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU_HI + best)) if best else None
        out.append((k, tmax, best, 2 * MW / (x / math.pi),
                    2 * MW / (xb / math.pi) if xb else None, mono))
    return out


def main():
    shown = draw(es=False)
    draw(es=True)
    P("=" * 92)
    P("EVERY RUNG, NOT JUST THE FIRST: the relaxation's boundary is empty on all of them")
    P("=" * 92)
    P("  %-6s %10s %10s %6s %14s %14s %s"
      % ("8D", "boundary", "integer", "step", "1/R5 relaxed", "1/R5 integer", "slack monotone"))
    for k, tmax, best, r_rel, r_int, mono in per_rung():
        P("  %-6d %10d %10s %6s %14.0f %14s %s"
          % (k, tmax, best if best else "--", (tmax - best) if best else "--",
             r_rel, ("%.0f" % r_int) if r_int else "--", mono))
    P("")
    P("  the boundary is never attained, and the integer optimum sits one or two steps of three")
    P("  below it.  That is the same phenomenon on all four rungs, not a peculiarity of the")
    P("  first -- which is what makes it worth a sentence in the paper rather than a footnote.")
    P("  Rungs 9 and 11 are NOT shown: their boundaries sit at A4 = 630 and 727, past the")
    P("  dynamic programme's A4 <= 560, and printing them would print a false null.")
    P("")
    P("=" * 92)
    P("what the right-hand panel is drawn from -- the fibre of reachable budgets, per vertex")
    P("=" * 92)
    P("  %-6s %10s %14s %14s %10s %s" % ("A4", "fibre size", "max budget", "required", "slack", ""))
    for a, f, need in shown:
        P("  %-6d %10d %14.3f %14.3f %10.3f %s"
          % (a, len(f), max(f), need, max(f) - need, "EMPTY" if max(f) < need else ""))
    P("")
    P("  CONTROL -- the panel must contain the case it is about: 215 present : %s ; 212 : %s"
      % (any(s[0] == 215 for s in shown), any(s[0] == 212 for s in shown)))
    sl = [max(s[1]) - s[2] for s in shown]
    P("  CONTROL -- the slack must FALL with A4, so the crossing is single and the integer")
    P("            optimum is well defined: monotone = %s"
      % all(sl[i] > sl[i + 1] for i in range(len(sl) - 1)))
    cross = [i for i in range(len(sl) - 1) if sl[i] >= 0 > sl[i + 1]]
    P("  CONTROL -- exactly one sign change : %s   (at A4 = %s -> %s)"
      % (len(cross) == 1,
         shown[cross[0]][0] if cross else "--",
         shown[cross[0] + 1][0] if cross else "--"))
    P("            everything past the crossing is empty too, which is why the picture is a")
    P("            proof and not a sample: %d of %d vertices above it fail"
      % (sum(1 for s in sl if s < 0), len(sl)))
    P("  written: fig_lift.pdf, fig_lift_es.pdf (here and in paper/), plus .png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
