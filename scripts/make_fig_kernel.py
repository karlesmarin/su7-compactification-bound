#!/usr/bin/env python3
"""make_fig_kernel.py -- one multiplet and forty-two of them, drawn on top of each other.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS FIGURE EXISTS
----------------------
Theorem (complete invariants) says two bulk contents have the same one-loop potential exactly when
they agree on five integers.  Said in words that is a rank argument, and a reader is entitled to
disbelieve it.  Drawn, it is one picture that settles it and one that explains it.

  LEFT -- the degeneracy, on the exact polylogarithmic potential.  A single 48(+,+) against
  24 x 7(+,+) + 18 x 7(+,-): one multiplet against forty-two, different representations,
  different dimensions, and ONE CURVE.  The inset is the residual, which sits at machine zero
  across the whole fundamental domain, so this is an identity and not a near-miss.

  RIGHT -- why it happens.  The six functions Re Li_5(s e^{i c pi a}), c in {1,2,3}, s = +-1, look
  like six directions and are five: their singular values are five of order thirty and one at
  1e-12.  The missing one is the duplication formula at s = 5, which is CCD24's eq. (2.11).  A
  content is invisible to the potential exactly when it is invisible to those five.

The two panels are the two halves of the theorem: the left is what it asserts, the right is why.

Run:  python make_fig_kernel.py    -> fig_kernel.pdf, fig_kernel_es.pdf, .png (here and paper/)
"""
import math
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
GREEN, GREY, INK = "#2E8B57", "#8A8A8A", "#222222"
MUTED = "#6E6E6E"
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]
KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]

# the pair: 48(+,+) = 24 x 7(+,+) + 18 x 7(+,-), the second of the three relations of eq. fivegen
LEFT_V = [0, 0, 0, 0, 1, 0, 0, 0]
RIGHT_V = [24, 18, 0, 0, 0, 0, 0, 0]


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#BBBBBB")
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)


def draw(es=False):
    T = lambda en, sp: sp if es else en
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.25),
                                   gridspec_kw=dict(width_ratios=[1.32, 1.0]))

    # ---------------------------------------------------------------- left: one curve, two contents
    a = np.linspace(0.0, 1.0, 1600)
    fA = F(cont(LEFT_V), a)
    fB = F(cont(RIGHT_V), a)
    axL.plot(a, fA, color=BLUE, lw=3.4, alpha=0.9,
             label=T(r"$\mathbf{48}^{(+,+)}$  — 1 multiplet",
                     r"$\mathbf{48}^{(+,+)}$  — 1 multiplete"))
    axL.plot(a, fB, color=RED, lw=1.5, ls=(0, (5, 3)),
             label=T(r"$24\times\mathbf{7}^{(+,+)}+18\times\mathbf{7}^{(+,-)}$  — 42 multiplets",
                     r"$24\times\mathbf{7}^{(+,+)}+18\times\mathbf{7}^{(+,-)}$  — 42 multipletes"))
    axL.set_xlabel(T(r"Wilson-line phase  $\alpha$", r"fase de Wilson  $\alpha$"), fontsize=9.5)
    axL.set_ylabel(T(r"$F(\alpha)$", r"$F(\alpha)$"), fontsize=9.5)
    axL.set_title(T("Different content. One potential.",
                    "Contenido distinto. Un solo potencial."),
                  loc="left", fontsize=10.5, color=INK, fontweight="bold", pad=9)
    axL.legend(fontsize=8.4, loc="upper left", framealpha=0.92)
    axL.set_xlim(0, 1)
    tidy(axL)

    # The residual sits at 1e-11, NOT at double precision, and saying "machine zero" here would
    # suggest the identity is only nearly true.  It is exactly true: dimension_five.sage proves
    # it over QQ.  What the inset shows is the floor of THIS evaluation, the 600-winding
    # truncation of the polylogarithm, and the label says so.
    # np.clip, not np.maximum-then-ylim: a point BELOW the axis floor still enters the PDF as a
    # path far outside the axes, and check_layout measures ink, not axes.  The first version of
    # this figure put 7345 pt of ink on an 842 pt page for exactly that reason.
    ins = axL.inset_axes([0.50, 0.13, 0.46, 0.30])
    d = np.abs(fA - fB)
    ins.semilogy(a, np.clip(d, 1.2e-15, 8e-8), color=GREY, lw=1.1)
    ins.set_ylim(1e-15, 1e-7)
    ins.set_xlim(0, 1)
    ins.set_title(T(r"residual $|\Delta F|$  —  the truncation floor",
                    r"residuo $|\Delta F|$  —  el suelo de truncación"), fontsize=7.2,
                  color=MUTED, pad=3)
    ins.tick_params(labelsize=6.4, colors=MUTED)
    for s in ("top", "right"):
        ins.spines[s].set_visible(False)
    ins.grid(alpha=0.15, linewidth=0.5)

    # ---------------------------------------------------------------- right: six functions, rank five
    ys = np.linspace(0.0, 1.0, 2001)
    G = np.array([basis(ys, s, c) for s, c in KEYS])
    sv = np.linalg.svd(G, compute_uv=False)
    idx = np.arange(1, 7)
    # bottom=FLOOR, explicitly.  A log-scale bar drawn from the default bottom=0 has its base at
    # log(0) = -inf, and the rectangle that reaches the PDF is unbounded -- the other half of the
    # 7345 pt of ink.  Give it a floor and the bar is a bar.
    FLOOR = 3e-16
    axR.bar(idx, np.clip(sv, FLOOR, None) - FLOOR, bottom=FLOOR,
            color=[BLUE] * 5 + [RED], width=0.62, alpha=0.9)
    axR.set_yscale("log")
    axR.set_ylim(FLOOR, 3e2)
    axR.set_xticks(idx)
    axR.set_xticklabels([str(i) for i in idx])
    axR.set_xlabel(T("singular value, ordered", "valor singular, ordenado"), fontsize=9.5)
    axR.set_ylabel(T("magnitude", "magnitud"), fontsize=9.5)
    axR.set_title(T("Six functions, five directions.", "Seis funciones, cinco direcciones."),
                  loc="left", fontsize=10.5, color=INK, fontweight="bold", pad=9)
    axR.annotate(T("the sixth is not there", "la sexta no está"),
                 xy=(6, max(sv[-1], 1e-15) * 4), xytext=(4.5, 3e-8),
                 fontsize=8.6, color=RED, ha="center",
                 arrowprops=dict(arrowstyle="->", color=RED, lw=1.1))
    # ONE line, and in the empty wedge above bars 4-6.  The two-line version sat under bars 1
    # and 2 and read as a smudge at page width -- which only showed up on the rendered page, not
    # in the PNG.  The full duplication formula belongs in the caption, and is there.
    axR.text(3.32, 5.0, T(r"$g^{(+,1)}+g^{(-,1)}=g^{(+,2)}/16$",
                          r"$g^{(+,1)}+g^{(-,1)}=g^{(+,2)}/16$"),
             fontsize=9.0, color=MUTED, va="bottom", ha="left")
    tidy(axR)

    fig.tight_layout()
    stem = "fig_kernel_es" if es else "fig_kernel"
    here = os.path.dirname(os.path.abspath(__file__))
    # pad_inches is not cosmetic here either: a bare tight bbox crops to the ink and an inset
    # tick label can spill past the declared BBox, which check_layout then reports as text
    # outside the text block.
    for d_ in (here, os.path.join(here, "paper")):
        fig.savefig(os.path.join(d_, stem + ".pdf"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(os.path.join(here, stem + ".png"), dpi=150, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return float(np.max(d)), sv


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    worst, sv = draw(es=False)
    draw(es=True)
    P("")
    P("=" * 92)
    P("WHAT THE FIGURE ASSERTS, as numbers")
    P("=" * 92)
    P("  left panel  : 48(+,+)  against  24 x 7(+,+) + 18 x 7(+,-)")
    P("     multiplets            : 1 against 42")
    P("     max |F difference| over 1600 points of [0,1] : %.2e" % worst)
    P("     double precision epsilon                     : %.2e" % np.finfo(float).eps)
    assert worst < 1e-9, "the pair is not degenerate -- the figure would be a lie"
    P("")
    P("  right panel : singular values of the six basis functions")
    P("     %s" % "  ".join("%.3e" % s for s in sv))
    P("     rank at 1e-10 relative : %d" % int(np.sum(sv > sv[0] * 1e-10)))
    assert int(np.sum(sv > sv[0] * 1e-10)) == 5, "the rank is not five -- the figure would be a lie"
    P("")
    P("  written: fig_kernel.pdf, fig_kernel_es.pdf (here and in paper/), plus .png")
