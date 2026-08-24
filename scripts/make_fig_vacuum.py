#!/usr/bin/env python3
"""make_fig_vacuum.py -- the false vacuum, drawn, because one curve says it and no sentence does.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS FIGURE EXISTS
----------------------
vacuum_constraint.py established that the content attaining the 10.01 TeV ceiling has its
electroweak minimum at alpha = 0.0161 -- correct, and predicted by the closed form to 0.02 % --
while the potential at the OTHER symmetric point, alpha = 1, is far below it.  Said in numbers
that is -204.99 against +111.44 and it means nothing to the eye.  Plotted over the whole
fundamental domain it is one glance:

  LEFT.  The exact polylogarithmic F(alpha) on 0 < alpha <= 1, for a published row and for the
  ceiling's witness, each shifted so F(0) = 0 so the two are comparable.  The published row rises
  to the right: its electroweak minimum is the deepest point of the domain.  The witness falls off
  a cliff: alpha = 1 is deeper, and the minimum the ceiling was built on is a false one.  The
  inset is the same witness near the origin, where the minimum is perfectly real.

  RIGHT.  What that costs.  The largest logarithmic budget available at each vertex of the rung
  8D = 1, with and without W > 0, against the budget the Higgs window demands.  The unconstrained
  curve crosses at A4 = 212; the constrained one does not catch up until A4 = 104.  The distance
  between the two crossings is 0.8 TeV of ceiling.

Run:  python make_fig_vacuum.py   -> fig_vacuum.pdf, fig_vacuum_es.pdf, .png  (into . and paper/)
"""
import math
import os
import pathlib
import sys
from fractions import Fraction as Fr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
GREEN, GREY = "#2E8B57", "#8A8A8A"
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

Z5 = 1.0369277551433699
MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_HI = 127.0
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]

_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])
U_OFF, LN2, LN3 = -19.5, math.log(2), math.log(3)
WV = [1, -1, 3, -3, 6, -6, 9, -9]
W0 = -1

WIT212 = [0] * 8
WIT212[0], WIT212[1], WIT212[3] = 17, 2, 57
WIT104 = [0] * 8
WIT104[0], WIT104[3], WIT104[4], WIT104[5], WIT104[6] = 16, 1, 1, 4, 1
ROW1 = None


def gstar(t, k):
    x = math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU(MH_HI) + t))
    return t * (math.log(x) + 0.75) + 3 * MU(MH_HI), x


def best_at(Atgt, Ktgt, need_W):
    best = [None]
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, u, v, w):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            ww = w + m * WV[0] + W0
            if need_W and ww <= 0:
                return
            cand = (u + m * GQ_SYM[0][2] + U_OFF) * LN2 + v * LN3
            if best[0] is None or cand > best[0]:
                best[0] = cand
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            rec(i + 1, a - c * AV[j], k - c * KV[j], u + c * GQ_SYM[j][2],
                v + c * GQ_SYM[j][1], w + c * WV[j])

    rec(0, Atgt, Ktgt, 0, 0, 0)
    return best[0]


def draw(es, data):
    L = dict(
        tA=("El vacío que el techo daba por bueno" if es
            else "The vacuum the ceiling took for granted"),
        tB=("Lo que cuesta exigir que sea el verdadero" if es
            else "What requiring the true one costs"),
        xa=r"$\alpha$", ya=(r"$F(\alpha)-F(0)$"),
        row=("fila publicada (2)" if es else "published row (2)"),
        w212=("testigo del techo, $(212,1)$" if es else "ceiling's witness, $(212,1)$"),
        w104=("testigo con vacío verdadero, $(104,1)$" if es
              else "true-vacuum witness, $(104,1)$"),
        ins=("el testigo del techo, cerca del origen" if es else "the ceiling's witness, near the origin"),
        xb=r"$A_4$ en el peldaño $8D=1$" if es else r"$A_4$ on the rung $8D=1$",
        yb=("presupuesto menos exigido" if es else "budget minus required"),
        free=("sin la condición" if es else "no condition"),
        wpos=("con $W>0$" if es else "with $W>0$"),
    )
    fig, ax = plt.subplots(1, 2, figsize=(12.8, 4.5))

    # ---- left: the potential over the whole domain
    xs = np.linspace(1e-6, 1.0, 4000)
    for v, col, lab in ((data["row"], BLUE, L["row"]),
                        (WIT104, GREEN, L["w104"]),
                        (WIT212, RED, L["w212"])):
        c = cont(v)
        y = F(c, xs)
        y0 = float(F(c, np.array([1e-9]))[0])
        ax[0].plot(xs, y - y0, color=col, lw=1.8, label=lab)
    ax[0].axhline(0, color=GREY, lw=0.7, ls=":")
    ax[0].axvline(1.0, color=GREY, lw=0.7, ls=":")
    ax[0].set_xlabel(L["xa"])
    ax[0].set_ylabel(L["ya"])
    ax[0].set_title(L["tA"], fontsize=11)
    ax[0].legend(fontsize=7.5, loc="lower left", framealpha=0.92)
    ax[0].grid(alpha=0.16, lw=0.5)

    # the inset is the ceiling's witness ALONE.  Drawing all three here was useless: the
    # published row's curvature is two orders larger and flattens the other two into the axis,
    # and the whole point of the inset is that the red curve's minimum is perfectly real.
    axi = ax[0].inset_axes([0.50, 0.58, 0.46, 0.38])
    xi = np.linspace(1e-6, 0.032, 2000)
    c = cont(WIT212)
    y0 = float(F(c, np.array([1e-9]))[0])
    yi = F(c, xi) - y0
    axi.plot(xi, yi, color=RED, lw=1.4)
    axi.axhline(0, color=GREY, lw=0.6, ls=":")
    j = int(np.argmin(yi))
    axi.scatter([xi[j]], [yi[j]], s=22, c=RED, zorder=5)
    # %.4f y no %.5f: xi es una rejilla de 2000 puntos sobre [0, 0.032], asi que su paso es
    # 1.6e-5 y el QUINTO decimal no esta resuelto.  Con %.5f el rotulo decia 0.01614 mientras
    # el cuerpo del articulo dice 0.016129 ---el minimizador de verdad, no el de la rejilla---
    # y el lector veia discrepar la figura con el texto en un digito que la figura no puede
    # saber.  Cuatro decimales es lo que esta rejilla sostiene, y coinciden.
    axi.annotate(r"$\alpha=%.4f$" % xi[j], (xi[j], yi[j]), textcoords="offset points",
                 xytext=(6, -12), fontsize=6.5, color=RED)
    axi.set_title(L["ins"], fontsize=7.5)
    axi.tick_params(labelsize=6.5)

    # ---- right: the two budget curves
    aa, s_free, s_w = [], [], []
    for A4t in range(215, 89, -3):
        Gs, _ = gstar(A4t, 1)
        need = float(Fr(25, 12) * A4t) - Gs
        a = best_at(A4t - A0, 1 - K0, False)
        b = best_at(A4t - A0, 1 - K0, True)
        if a is None or b is None:
            continue
        aa.append(A4t)
        s_free.append(a - need)
        s_w.append(b - need)
    aa = np.array(aa, float)
    ax[1].axhline(0, color=AMBER, lw=1.8, zorder=2)
    ax[1].plot(aa, s_free, color=GREY, lw=1.8, label=L["free"], zorder=3)
    ax[1].plot(aa, s_w, color=GREEN, lw=1.8, label=L["wpos"], zorder=3)
    for x0, col in ((212, GREY), (104, GREEN)):
        ax[1].axvline(x0, color=col, lw=1.0, ls="--", alpha=0.8, zorder=1)
        ax[1].annotate("%d" % x0, (x0, ax[1].get_ylim()[0]), textcoords="offset points",
                       xytext=(3, 6), fontsize=8, color=col)
    ax[1].set_xlabel(L["xb"])
    ax[1].set_ylabel(L["yb"])
    ax[1].set_title(L["tB"], fontsize=11)
    ax[1].legend(fontsize=8, loc="upper left", framealpha=0.92)
    ax[1].grid(alpha=0.16, lw=0.5)
    ax[1].invert_xaxis()

    fig.tight_layout()
    stem = "fig_vacuum_es" if es else "fig_vacuum"
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (here, os.path.join(here, "paper")):
        fig.savefig(os.path.join(d, stem + ".pdf"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(os.path.join(here, stem + ".png"), dpi=150, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return aa, s_free, s_w


def main():
    v = [0] * 8
    for rep, e, ep, mult in T1[1][1]:
        v[REPS.index((rep, e, ep))] += mult
    data = {"row": v}
    aa, sf, sw = draw(False, data)
    draw(True, data)
    P("=" * 92)
    P("what the panels are drawn from")
    P("=" * 92)
    P("  %-30s %14s %14s %10s" % ("content", "F(0)", "F(1)", "verdict"))
    for lab, vv in (("published row (2)", v), ("witness (212,1)", WIT212),
                    ("witness (104,1)", WIT104)):
        c = cont(vv)
        f = F(c, np.array([1e-9, 1.0]))
        P("  %-30s %14.3f %14.3f %10s"
          % (lab, f[0], f[1], "true" if f[1] > f[0] else "FALSE"))
    P("")
    P("  CONTROL -- the panel must contain both crossings, and only the constrained curve may")
    P("  cross late:")
    cf = [aa[i] for i in range(len(aa)) if sf[i] >= 0]
    cw = [aa[i] for i in range(len(aa)) if sw[i] >= 0]
    P("     unconstrained crosses at A4 = %s" % (max(cf) if cf else "--"))
    P("     with W > 0        crosses at A4 = %s" % (max(cw) if cw else "--"))
    P("     the two crossings are %s, which is the whole point of the panel"
      % ("different" if (cf and cw and max(cf) != max(cw)) else "THE SAME -- redraw"))
    P("  written: fig_vacuum.pdf, fig_vacuum_es.pdf (here and in paper/), plus .png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
