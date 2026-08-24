#!/usr/bin/env python3
"""make_figures_vii.py -- the four figures of Part VII.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  Draws the four figures of the Part VII paper from the archived outputs; nothing is recomputed
  from scratch and nothing is transcribed by hand.

Numbers are READ from outputs/*.json, each produced by the script of the same name.  The term
tables come from amin_closed_form.py, which is itself extracted from ../part_vi/su7_anchor_mh.py,
so a change upstream cannot silently desynchronise a figure.  [[save-the-outputs-not-just-the-scripts]]

  fig_ladder     the expansion as ONE geometric ladder.  Every term of the potential walks it as
                 sigma^k with sigma = c^2/4, so charge two is the fixed point -- the flat spine of
                 the figure -- and everything else fans away from it.  The right panel is the
                 gauge sector's own three rungs and the exact relation between them.
  fig_ceiling    the headline.  The best 1/R_5 against the size of the content, the certified
                 ceiling as a band, the collider exclusions underneath, and the per-D ceiling
                 showing that the maximum sits on the quantum D = 1/8.
  fig_parity     why there is no polynomial: the same Fourier sum at an even and an odd index,
                 and the arithmetic tell -- zeta(p)/pi^p rational exactly when p is even.
  fig_surface    the ceiling's geometry in 3D: 1/R_5 over the (A_4, 8D) plane once m_h is pinned,
                 the feasibility boundary drawn on the surface, and the certified maximum.

Palette inherited from Part VI so the series reads as one, and re-validated rather than trusted:
  node scripts/validate_palette.js "#3A86C8,#E0A030,#C0392B" --mode light
  -> lightness band PASS, chroma floor PASS, CVD separation dE 20.7 (deutan) / 22.7 (tritan) PASS,
     normal-vision floor 24.3 PASS, contrast WARN on the amber -- discharged, as in Part VI, by a
     direct label on every mark it fills.
Sequential use (height, magnitude) is ONE hue light-to-dark.  Diverging use (the sign of D) is the
same blue/red pair with a neutral midpoint, consistently with Parts V and VI.
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ES = "--es" in sys.argv
T = lambda en, es: es if ES else en
SUF = "_es" if ES else ""
OUT = os.path.join(HERE, "paper")
os.makedirs(OUT, exist_ok=True)
J = lambda n: json.load(open(os.path.join(HERE, "outputs", n)))

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
STROKE, INK, MUTED = "#1F4E79", "#1F2933", "#6B7280"
GRID, SURF = "#E5E7EB", "#FCFCFB"
# Every italic caption below sits inside the axes, where the data is.  Four of them were printing
# straight through curves, dot rows and step lines -- legible on the screen where they were placed
# and not legible in the paper.  A caption is a claim; if it cannot be read it is not made.
# The box alone does NOT fix it: matplotlib draws text at zorder 3 by default, under every series
# here, so the box goes behind the data and the glyphs get overprinted anyway.  Both are needed.
CAPBOX = dict(fc=SURF, ec="none", alpha=0.88, pad=1.6)
CAPZ = 12
SEQ = LinearSegmentedColormap.from_list("seq", ["#EAF2FA", "#9CC4E4", BLUE, STROKE, "#12314C"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.facecolor": SURF, "figure.facecolor": SURF,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
})


def tidy(ax):
    ax.tick_params(length=2.5, width=0.6)
    for s in ax.spines.values():
        s.set_linewidth(0.7)


# =========================================================== 1. the ladder
def fig_ladder():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.5),
                                   gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.28))

    # -- left: the geometric modes.  sigma = c^2/4 for every term, and additionally c^2 if s = -1.
    modes = [(0.25, "$c=1$", MUTED, "-"), (1.0, "$c=2$", STROKE, "-"),
             (2.25, "$c=3$", AMBER, "-"), (4.0, "$c=2$, antiper.", BLUE, "--"),
             (9.0, "$c=3$, antiper.", RED, "--")]
    ks = np.array([0, 1, 2])
    for sig, lab, col, ls in modes:
        y = sig ** ks
        lw = 2.6 if sig == 1.0 else 1.7
        axL.plot(ks, y, ls, color=col, lw=lw, marker="o", ms=6 if sig == 1.0 else 4.5,
                 mfc=col, mec=SURF, mew=1.1, zorder=4 if sig == 1.0 else 3,
                 solid_capstyle="round")
        axL.annotate(r"$\sigma=%s$  %s" % (("1/4" if sig == .25 else
                                            ("9/4" if sig == 2.25 else "%g" % sig)), lab),
                     xy=(2, y[-1]), xytext=(6, 0), textcoords="offset points",
                     va="center", ha="left", fontsize=8,
                     color=col, fontweight="bold" if sig == 1.0 else "normal")
    axL.axhline(1.0, color=STROKE, lw=0.7, ls=":", zorder=1)
    axL.set_yscale("log")
    axL.set_xticks(ks)
    axL.set_xticklabels([T("$k=0$\nthe value\n$p=5$, $\\zeta(5)$", "$k=0$\nel valor\n$p=5$, $\\zeta(5)$"),
                         T("$k=1$\nthe curvature\n$p=3$, $\\zeta(3)$", "$k=1$\nla curvatura\n$p=3$, $\\zeta(3)$"),
                         T("$k=2$\nthe fourth moment\n$p=1$, the pole", "$k=2$\ncuarto momento\n$p=1$, el polo")],
                        fontsize=7.6)
    axL.set_xlim(-0.28, 2.95)
    axL.set_ylabel(T("weight along the ladder,  $\\sigma^{\\,k}$",
                     "peso a lo largo de la escalera,  $\\sigma^{\\,k}$"))
    axL.set_title(T("Every term walks one geometric ladder", "Cada término recorre una sola escalera geométrica"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    axL.text(-0.22, 1.55, T("charge two is the FIXED POINT:\nsame weight on all three rungs",
                            "la carga dos es el PUNTO FIJO:\nel mismo peso en los tres peldaños"),
             fontsize=7.6, color=STROKE, va="bottom", ha="left", style="italic", bbox=CAPBOX, zorder=CAPZ)
    tidy(axL)

    # -- right: the gauge sector's own three rungs, from ladder_closed_form.json
    LC = J("ladder_closed_form.json")["reference"]["rungs"]
    vals = [LC["5"], LC["3"], LC["1"]]
    v = [float(x) for x in vals]
    cols = [BLUE if z > 0 else RED for z in v]
    axR.bar(range(3), v, color=cols, width=0.55, edgecolor=SURF, linewidth=1.6, zorder=3)
    for i, z in enumerate(v):
        axR.annotate("%+d" % z, xy=(i, z), xytext=(0, 7 if z > 0 else -13),
                     textcoords="offset points", ha="center", fontsize=10,
                     fontweight="bold", color=cols[i])
    axR.axhline(0, color=MUTED, lw=0.8, zorder=2)
    axR.set_xticks(range(3))
    axR.set_xticklabels([T("value\n$p=5$", "valor\n$p=5$"), T("curvature\n$p=3$", "curvatura\n$p=3$"),
                         T("4th moment\n$p=1$", "4º momento\n$p=1$")], fontsize=8)
    axR.set_ylabel(T("gauge contribution  $O(p)$", "aportación gauge  $O(p)$"))
    axR.set_ylim(-52, 30)
    d1, d2 = v[0] - v[1], v[1] - v[2]
    axR.annotate("", xy=(0.5, 22), xytext=(0.5, 22), )
    axR.plot([0, 1], [20, 20], color=STROKE, lw=1.0)
    # the bracket used to run the full [1, 2] and its right end came out underneath the bold
    # "-36" the bar carries; shortened, it now stops clear of it.
    axR.plot([1.12, 1.88], [-44, -44], color=STROKE, lw=1.0)
    axR.text(0.5, 21.5, "$O(5)-O(3)=%+d$" % d1, ha="center", fontsize=8, color=STROKE)
    axR.text(1.5, -42.5, "$O(3)-O(1)=%+d$" % d2, ha="center", fontsize=8, color=STROKE,
             bbox=CAPBOX, zorder=CAPZ)
    axR.text(1.0, -50, T("and $%d = 4\\times %d$ identically, on all 4096 assignments" % (d1, d2),
                         "y $%d = 4\\times %d$ idénticamente, en las 4096 asignaciones" % (d1, d2)),
             ha="center", fontsize=8.2, color=INK, style="italic")
    axR.set_title(T("Odd, odd, even --- and two of the three are free",
                    "Impar, impar, par --- y sólo dos de los tres son libres"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    tidy(axR)

    fig.savefig(os.path.join(OUT, "fig_ladder%s.pdf" % SUF), bbox_inches="tight")
    plt.close(fig)


# =========================================================== 2. the ceiling
def fig_ceiling():
    C = J("ceiling_ilp.json")
    PR = J("prediction.json")
    curve = C["size_curve"]
    ceil = C["ceiling_GeV"] / 1000.0
    band_hi = PR["ceiling_band"][1] / 1000.0

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.4),
                                   gridspec_kw=dict(width_ratios=[1.45, 1], wspace=0.26))

    N = [r["N"] for r in curve]
    y = [r["invR"] / 1000.0 for r in curve]
    axL.axhspan(ceil, band_hi, color=BLUE, alpha=0.10, zorder=1, lw=0)
    axL.axhline(ceil, color=STROKE, lw=1.8, zorder=4)
    axL.axhspan(0, 5.0, color=RED, alpha=0.085, zorder=1, lw=0)
    axL.axhline(5.0, color=RED, lw=1.0, ls=":", zorder=3)
    axL.step(N, y, where="post", color=BLUE, lw=2.0, zorder=5, solid_joinstyle="round")
    axL.plot(N, y, "o", ms=5.5, mfc=BLUE, mec=SURF, mew=1.2, zorder=6)
    for r in curve:
        if r["N"] in (6, 8):
            col = RED if r["N"] == 6 else STROKE
            axL.plot(r["N"], r["invR"] / 1000.0, "o", ms=9, mfc=col, mec=SURF, mew=1.6, zorder=7)
            axL.annotate(T("$N=%d$:  %.2f TeV\n$D=%s$" % (r["N"], r["invR"] / 1000.0,
                                                          "%d/8" % round(8 * r["D"])),
                           "$N=%d$:  %.2f TeV\n$D=%s$" % (r["N"], r["invR"] / 1000.0,
                                                          "%d/8" % round(8 * r["D"]))),
                         xy=(r["N"], r["invR"] / 1000.0),
                         # the N=6 label used to be offset to the RIGHT, where it ran straight
                         # through the riser at N=7 and through the ATLAS/CMS caption; both are
                         # to its right, so it goes left instead.  The boxes are what make the
                         # numbers readable where they must sit over the step line at all.
                         xytext=(-11, 8) if r["N"] == 6 else (11, -34),
                         ha="right" if r["N"] == 6 else "left",
                         textcoords="offset points",
                         fontsize=8.4, color=col, fontweight="bold",
                         bbox=dict(fc=SURF, ec="none", alpha=0.86, pad=1.4))
    axL.text(13.8, ceil + 0.35, T("certified ceiling  %.2f TeV" % ceil,
                                  "techo certificado  %.2f TeV" % ceil),
             ha="right", fontsize=8.6, color=STROKE, fontweight="bold")
    # NO "anchor band".  La S2 demuestra que el residuo del anclaje NO es una banda: es una
    # TENDENCIA monotona en D, con p = 0.017, y llamarla banda invita justo a lo que esa seccion
    # prohibe -- propagarla como un multiplicador.  El texto se corrigio y la figura se quedo con
    # la palabra vieja, que es la que el lector recuerda.  [[a-figure-is-an-unaudited-document]]
    axL.text(13.8, (ceil + band_hi) / 2, T("anchor extrapolation", "extrapolación del anclaje"),
             ha="right", fontsize=7.8, color=STROKE, style="italic", alpha=0.9)
    # this caption sat at y = 2.2, which is exactly where the step line and the N=6 marker are.
    # It printed through both.  The band below y = 1.5 is empty across the whole panel.
    # NOT "excluded by ATLAS/CMS resonance searches", which is what this said and which the paper
    # itself refutes two sections later: those limits are set on a NARROW benchmark resonance and
    # this object has Gamma/M = 0.16 and an unfixed BR(G_1 -> jj), so they do not transport to it.
    # A figure that says "excluded" where the text says "cannot be carried across" is the text
    # losing an argument to its own picture.  [[a-figure-is-an-unaudited-document]]
    axL.text(4.0, 0.78, T("generic benchmark resonance limits\n(not an $SU(7)$-specific exclusion)",
                          "límites genéricos de resonancia de referencia\n(no una exclusión propia de $SU(7)$)"),
             fontsize=7.8, color=RED, va="center", ha="left")
    axL.set_xlabel(T("$N$, multiplets in the bulk content", "$N$, multipletes en el contenido"))
    axL.set_ylabel(T("best $1/R_5$ inside the Higgs window  (TeV)",
                     "mejor $1/R_5$ dentro de la ventana de Higgs  (TeV)"))
    axL.set_xlim(3.6, 14.4)
    axL.set_ylim(0, band_hi * 1.04)
    axL.set_title(T("The cap on the content was doing the work",
                    "La tapa sobre el contenido era la que trabajaba"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    tidy(axL)

    # -- right: the per-D ceiling.  the maximum sits on the quantum.
    pk = C["per_k"]
    kk = [r["k8D"] for r in pk]
    vv = [r["invR"] / 1000.0 for r in pk]
    axR.plot(kk, vv, "-", color=BLUE, lw=1.8, zorder=3)
    axR.plot(kk, vv, "o", ms=5, mfc=BLUE, mec=SURF, mew=1.1, zorder=4)
    axR.plot(kk[0], vv[0], "o", ms=10, mfc=AMBER, mec=STROKE, mew=1.5, zorder=5)
    axR.annotate(T("$D=1/8$, the quantum\n%.2f TeV" % vv[0], "$D=1/8$, el cuanto\n%.2f TeV" % vv[0]),
                 xy=(kk[0], vv[0]), xytext=(20, -20), textcoords="offset points",
                 fontsize=8.6, color=STROKE, fontweight="bold")
    # The dotted line used to be the LAST SWEPT rung, i.e. "where the scan happened to stop".
    # It is now the asymptotic ray: r_infinity solves r(A - (1/2)ln r) = nu in closed form and
    # M_infinity = sqrt(prefactor * r_inf).  Drawing the limit instead of the last sample is the
    # difference between a flattening and a value.  Read from asymptotic_ray.json so the picture
    # cannot drift from the run.
    _ray = json.load(open(os.path.join(HERE, "outputs", "asymptotic_ray.json"), encoding="utf-8"))
    minf = _ray["M_infinity"] / 1000.0
    axR.axhline(minf, color=RED, lw=1.1, ls="--", zorder=2)
    # left of the curve's descent, above the line: the only wedge on this panel that is empty in
    # both editions.  Placed by eye against the rendered PNG, not guessed.
    axR.text(1.35, minf + 0.62, T(r"$M_\infty = %.3f$ TeV, a Lambert $W$" % minf,
                                  r"$M_\infty = %.3f$ TeV, una Lambert $W$" % minf),
             ha="left", fontsize=8.2, color=RED, style="italic", bbox=CAPBOX, zorder=CAPZ)
    axR.text(1.35, minf + 0.24, T("never reached", "nunca alcanzada"),
             ha="left", fontsize=7.4, color=MUTED, style="italic", bbox=CAPBOX, zorder=CAPZ)
    axR.set_xscale("log")
    # el Teorema 1 es CONDICIONAL, asi que "by the theorem" a secas lee como
    # incondicional: es el Corolario, con la semilla publicada, el que lo cierra.
    axR.set_xlabel(T("$8D$  (odd on the seed of [2])", "$8D$  (impar con la semilla de [2])"))
    axR.set_ylabel(T("ceiling at that $D$  (TeV)", "techo a ese $D$  (TeV)"))
    axR.set_title(T("Monotone in $D$: the maximum is on the quantum",
                    "Monótono en $D$: el máximo está en el cuanto"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    tidy(axR)

    fig.savefig(os.path.join(OUT, "fig_ceiling%s.pdf" % SUF), bbox_inches="tight")
    plt.close(fig)


# =========================================================== 3. the parity dichotomy
def fig_parity():
    import mpmath as mp
    mp.mp.dps = 25
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.1),
                                   gridspec_kw=dict(width_ratios=[1.4, 1], wspace=0.28))

    xs = np.linspace(0.001, 0.999, 320)
    s4 = np.array([float(mp.nsum(lambda n: mp.cos(2 * mp.pi * n * x) / n ** 4, [1, mp.inf]))
                   for x in xs])
    s5 = np.array([float(mp.nsum(lambda n: mp.cos(2 * mp.pi * n * x) / n ** 5, [1, mp.inf]))
                   for x in xs])
    b4 = np.array([float(-(2 * mp.pi) ** 4 * mp.bernpoly(4, x) / (2 * mp.factorial(4))) for x in xs])

    axL.plot(xs, s4, color=BLUE, lw=2.4, zorder=3,
             label=T("$p=4$ (even) --- the thermal case", "$p=4$ (par) --- el caso térmico"))
    axL.plot(xs, b4, color=STROKE, lw=1.1, ls=(0, (5, 3)), zorder=4,
             label=T("$-(2\\pi)^4 B_4(x)/48$ --- exact", "$-(2\\pi)^4 B_4(x)/48$ --- exacto"))
    axL.plot(xs, s5, color=RED, lw=2.4, zorder=3,
             label=T("$p=5$ (odd) --- the compactified case", "$p=5$ (impar) --- el caso compactificado"))
    axL.axhline(0, color=MUTED, lw=0.7)
    axL.set_xlabel(T("holonomy phase $x$", "fase de holonomía $x$"))
    axL.set_ylabel(T("$\\sum_n \\cos(2\\pi n x)/n^{\\,p}$", "$\\sum_n \\cos(2\\pi n x)/n^{\\,p}$"))
    leg = axL.legend(frameon=False, fontsize=8, loc="upper center", ncol=1,
                     handlelength=2.2, borderaxespad=0.2)
    for t in leg.get_texts():
        t.set_color(INK)
    axL.set_title(T("One index apart: a polynomial, or not",
                    "A un índice de distancia: polinomio, o no"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    # room made BELOW the trough for the caption: boxing it at s4.min()*0.92 put it on top of the
    # two minima, which is the one place on this panel a reader compares the curves.
    axL.set_ylim(s4.min() * 1.42, None)
    # NOT "no closed form": a Clausen function IS one, and this paper goes on to find a
    # Lambert-W closed form for the stationarity condition three sections later.  What the
    # x^4 log x branch point rules out is a POLYNOMIAL -- specifically a Bernoulli-polynomial
    # representation like the even case has -- and that is what the panel should say.
    axL.text(0.5, s4.min() * 1.20, T("the even curve IS a quartic --- Gross--Pisarski--Yaffe;\n"
                                     "the odd one is a Clausen function: no Bernoulli polynomial",
                                     "la curva par ES una cuártica --- Gross--Pisarski--Yaffe;\n"
                                     "la impar es de Clausen: ningún polinomio de Bernoulli"),
             ha="center", va="top", fontsize=7.8, color=INK, style="italic", bbox=CAPBOX, zorder=CAPZ)
    tidy(axL)

    # -- right: the arithmetic tell.  zeta(p)/pi^p rational exactly when p is even.
    ps = [2, 3, 4, 5, 6]
    rat = {2: "1/6", 4: "1/90", 6: "1/945"}
    for i, p in enumerate(ps):
        r = float(mp.zeta(p) / mp.pi ** p)
        even = p % 2 == 0
        col = BLUE if even else RED
        axR.plot([i], [r], "o", ms=13, mfc=col, mec=SURF, mew=1.6, zorder=4)
        axR.annotate(rat[p] if even else T("no relation", "sin relación"),
                     xy=(i, r), xytext=(0, 15), textcoords="offset points",
                     ha="center", fontsize=8.6, color=col,
                     fontweight="bold" if even else "normal")
        if not even:
            axR.annotate(T("$\\zeta(%d)$ --- rung of\nOUR ladder" % p,
                           "$\\zeta(%d)$ --- peldaño\nde NUESTRA escalera" % p),
                         xy=(i, r), xytext=(0, -34), textcoords="offset points",
                         ha="center", fontsize=7.6, color=RED, style="italic")
    axR.set_yscale("log")
    axR.set_xticks(range(len(ps)))
    axR.set_xticklabels(["$p=%d$" % p for p in ps], fontsize=9)
    axR.set_xlim(-0.6, len(ps) - 0.4)
    axR.set_ylabel(T("$\\zeta(p)\\,/\\,\\pi^{\\,p}$", "$\\zeta(p)\\,/\\,\\pi^{\\,p}$"))
    # NOT "rational exactly when p is even": for odd p the rationality of zeta(p)/pi^p is OPEN, and
    # the caption right below this panel says so.  The title was claiming more than the caption
    # admitted, which is the worst place for an overclaim to sit.  Carles caught it.
    axR.set_title(T("Provably rational for even $p$; odd $p$ unresolved",
                    "Demostradamente racional en $p$ par; $p$ impar, abierto"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    axR.text(2, 4e-4, T("integer-relation search at 40 digits", "búsqueda de relación entera a 40 dígitos"),
             ha="center", fontsize=7.6, color=MUTED, style="italic")
    tidy(axR)

    fig.savefig(os.path.join(OUT, "fig_parity%s.pdf" % SUF), bbox_inches="tight")
    plt.close(fig)


# =========================================================== 4. the 3D surface
def fig_surface():
    from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
    C = J("ceiling_ilp.json")
    PR = J("prediction.json")
    pk = C["per_k"]
    MW, Z3 = 80.4, 1.2020569031595943
    KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * C["g4"]
    mu = (C["mh_window"][1] / (KK * math.pi ** 2)) ** 2

    # 8D is plotted on a log axis (as u = log_3 of it): the whole structure lives at small D,
    # and a linear axis drowns it in the 1/sqrt(D) tail.
    kmax = 27
    tt = np.linspace(0, 1500, 220)
    kk = np.logspace(0, math.log10(kmax), 180)
    TT, KKg = np.meshgrid(tt, kk)
    U = np.log(KKg) / math.log(3.0)
    x = np.sqrt(12 * Z3 * (KKg / 8.0) / (6 * mu + TT))
    Z = 2 * math.pi * MW / x / 1000.0                # TeV

    # the feasibility boundary: largest A_4 at each 8D, from the certificate
    bk = np.array([r["k8D"] for r in pk], float)
    bt = np.array([r["A4"] for r in pk], float)
    keep = bk <= kmax
    bk, bt = bk[keep], bt[keep]
    bz = np.array([2 * math.pi * MW / math.sqrt(12 * Z3 * (k / 8.0) / (6 * mu + t)) / 1000.0
                   for k, t in zip(bk, bt)])

    mask = np.zeros_like(Z, dtype=bool)
    for j, kv in enumerate(kk):
        tlim = np.interp(kv, bk, bt)
        mask[j, :] = TT[j, :] > tlim
    Zm = np.ma.array(Z, mask=mask)

    fig = plt.figure(figsize=(9.6, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(TT, U, Zm, cmap=SEQ, rstride=2, cstride=3, linewidth=0.15,
                    edgecolors=(1, 1, 1, 0.16), antialiased=True, vmin=2.6, vmax=10.2, alpha=0.985)
    ax.contour(TT, U, Zm, levels=[3, 3.5, 4, 5, 6, 7, 8, 9, 10], colors=[STROKE], linewidths=0.55,
               alpha=0.55, offset=None)
    ax.plot(bt, np.log(bk) / math.log(3.0), bz, color=RED, lw=2.4, zorder=10)
    ax.plot(bt, np.log(bk) / math.log(3.0), np.zeros_like(bz), color=RED, lw=0.9, alpha=0.35, zorder=2)

    ct, ck = C["ceiling_A4"], C["ceiling_8D"]
    cz = C["ceiling_GeV"] / 1000.0
    ax.plot([ct, ct], [0, 0], [0, cz], color=AMBER, lw=1.6, zorder=11)
    ax.scatter([ct], [0], [cz], s=70, c=AMBER, edgecolors=STROKE, linewidths=1.3, zorder=12,
               depthshade=False)
    ax.text(ct + 120, 0.10, cz + 0.9,
            T("%.2f TeV  at  $A_4=%d$, $D=1/8$" % (cz, ct), "%.2f TeV  en  $A_4=%d$, $D=1/8$" % (cz, ct)),
            color=STROKE, fontsize=9, fontweight="bold", zorder=13)

    ax.set_xlabel(T("$A_4$, the fourth moment", "$A_4$, el cuarto momento"), labelpad=8)
    # corto A PROPOSITO: la condicion de semilla va en la nota al pie de la figura, donde
    # cabe.  Con la condicion dentro, la etiqueta del eje y se salia 5.4 pt del bloque de
    # texto en la edicion castellana.  [[the-pdf-text-layer-is-not-the-page]]
    ax.set_ylabel(T("$8D$  (odd [2], log)", "$8D$  (impar [2], log)"), labelpad=6)
    ax.set_zlabel(T("$1/R_5$  (TeV)", "$1/R_5$  (TeV)"), labelpad=4)
    ax.set_zlim(0, 11.5)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["1", "3", "9", "27"])
    ax.view_init(elev=20, azim=-58)
    ax.set_box_aspect((1.30, 0.95, 0.85))
    ax.xaxis.pane.set_facecolor(SURF); ax.yaxis.pane.set_facecolor(SURF)
    ax.zaxis.pane.set_facecolor("#F4F6F8")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_edgecolor(GRID)
        a._axinfo["grid"].update(color=GRID, linewidth=0.5)
    ax.set_title(T("Once $m_h$ is pinned the hierarchy is algebraic --- and the window ends it",
                   "Fijado $m_h$ la jerarquía es algebraica --- y la ventana la termina"),
                 loc="left", fontsize=11, color=INK, fontweight="bold", pad=2)
    handles = [Line2D([], [], color=RED, lw=2.4,
                      # NO es el borde de la ventana: la superficie ya esta dibujada a
                      # m_h = 127 GeV.  Es donde la identidad (I) se queda sin contenidos,
                      # que es como la nombran el pie y la nota al pie de esta figura.
                      label=T("where identity (I) runs out of contents",
                              "donde la identidad (I) se queda sin contenidos")),
               Line2D([], [], marker="o", ls="", mfc=AMBER, mec=STROKE, ms=8,
                      label=T("certified maximum", "máximo certificado"))]
    leg = ax.legend(handles=handles, frameon=False, fontsize=8.6, loc="upper right",
                    bbox_to_anchor=(1.02, 0.93))
    for t in leg.get_texts():
        t.set_color(INK)
    fig.text(0.055, 0.045,
             T("Surface: $1/R_5 = 2\\pi m_W\\sqrt{(6\\mu+A_4)/(12\\zeta(3)D)}$, identity (II), at "
               "$m_h = 127$ GeV.  Beyond the red edge no content satisfies identity (I).\n"
               "$8D$ is odd on the seed [2] print; on the candidate seed of the open-questions "
               "section it is even, and the staircase starts one step higher.",
               "Superficie: $1/R_5 = 2\\pi m_W\\sqrt{(6\\mu+A_4)/(12\\zeta(3)D)}$, identidad (II), a "
               "$m_h = 127$ GeV.  Más allá del borde rojo ningún contenido cumple la (I).\n"
               "$8D$ es impar con la semilla que imprime [2]; con la candidata de la sección de "
               "preguntas abiertas es par, y la escalera arranca un peldaño más arriba."),
             fontsize=7.8, color=MUTED, style="italic")
    fig.savefig(os.path.join(OUT, "fig_surface%s.pdf" % SUF), bbox_inches="tight")
    plt.close(fig)



# =========================================================== 5. the moment cone (Sage-certified)
def fig_cone():
    """The lattice the integer program runs on, with its geometry computed exactly.

    Everything structural here is READ from outputs/cone_moments.json, produced by
    cone_moments.sage: the extreme rays, the dual description, the index of the sublattice the
    eight multiplets reach, the Hilbert basis, and which of its elements the model cannot realise.
    Sage certifies the geometry; this only draws it.
    """
    K = J("cone_moments.json")
    C = J("ceiling_ilp.json")
    off = K["offset"]
    pts = np.array(K["points"], float)
    TMAX, KMAX = 620, 46
    sel = (pts[:, 0] <= TMAX) & (pts[:, 1] <= KMAX)
    pts = pts[sel]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.6),
                                   gridspec_kw=dict(width_ratios=[1.55, 1], wspace=0.24))

    # ---- left: the lattice AT THE SCALE WHERE THE ANSWER IS DECIDED.  Drawn over the whole cone
    # the points merge into lines and the cone constraint is not the binding one -- the window is.
    T0, T1, K1 = 150, 275, 13
    ct, ck = C["ceiling_A4"], C["ceiling_8D"]
    pk = C["per_k"]
    bk = np.array([r["k8D"] for r in pk], float)
    bt = np.array([r["A4"] for r in pk], float)
    inb = np.array([np.interp(k, bk, bt) for k in np.arange(1, K1 + 1, 2)])
    yy = np.arange(1, K1 + 1, 2)
    axL.fill_betweenx(yy, inb, T1, color=RED, alpha=0.07, lw=0, zorder=1)
    axL.plot(inb, yy, "-", color=RED, lw=2.2, zorder=6)

    sel2 = (pts[:, 0] >= T0) & (pts[:, 0] <= T1) & (pts[:, 1] <= K1)
    P2 = pts[sel2]
    feas = np.array([p[0] <= np.interp(p[1], bk, bt) for p in P2])
    axL.plot(P2[feas, 0], P2[feas, 1], "o", ms=4.2, mfc=BLUE, mec=SURF, mew=0.6,
             zorder=4, ls="")
    axL.plot(P2[~feas, 0], P2[~feas, 1], "o", ms=3.4, mfc="#FFFFFF", mec=MUTED, mew=0.8,
             zorder=3, ls="")
    axL.plot([ct], [ck], "o", ms=11, mfc=AMBER, mec=STROKE, mew=1.6, zorder=8)
    axL.annotate(T("the ceiling: $A_4=%d$, $8D=%d$\n%.2f TeV" % (ct, ck, C["ceiling_GeV"] / 1000),
                   "el techo: $A_4=%d$, $8D=%d$\n%.2f TeV" % (ct, ck, C["ceiling_GeV"] / 1000)),
                 xy=(ct, ck), xytext=(-14, 62), textcoords="offset points",
                 fontsize=8.6, color=STROKE, fontweight="bold", ha="right",
                 # the label sits at the height of the 8D = 5 row, and without a background the
                 # lattice dots run straight through the letters -- "the ceiling" came out
                 # unreadable at print size.  No gate sees this: check_layout measures bounding
                 # boxes, not legibility, so it took looking at the figure.
                 bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.88),
                 zorder=9,
                 arrowprops=dict(arrowstyle="-", color=STROKE, lw=0.8, zorder=9))
    axL.text(T1 - 4, K1 - 1.1, T("no content here:\nidentity (I) has no solution",
                                 "aquí no hay contenido:\nla identidad (I) no tiene solución"),
             fontsize=8, color=RED, ha="right", va="top", style="italic", bbox=CAPBOX, zorder=CAPZ)
    axL.text(T0 + 3, K1 - 1.1, T("admissible $(A_4,8D)$ --- one\ninteger point in six",
                                 "$(A_4,8D)$ admisibles --- un\npunto entero de cada seis"),
             fontsize=8, color=STROKE, ha="left", va="top", style="italic", bbox=CAPBOX, zorder=CAPZ)
    axL.set_xlim(T0, T1)
    axL.set_ylim(0, K1 + 0.7)
    axL.set_yticks(range(1, K1 + 1, 2))
    axL.set_xlabel(T("$A_4$, the fourth moment", "$A_4$, el cuarto momento"))
    axL.set_ylabel(T("$8D$  (odd on the seed of [2])", "$8D$  (impar con la semilla de [2])"))
    axL.set_title(T("The lattice, at the scale where the answer is decided",
                    "El retículo, a la escala en que se decide la respuesta"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    tidy(axL)

    # ---- right: why only one point in six is there
    axR.set_xlim(-0.6, 5.6)
    axR.set_ylim(-0.6, 5.6)
    for a in range(6):
        for b in range(6):
            ok = (b % 2 == 1) and ((a + b) % 3 == 0)
            axR.plot([a], [b], "o", ms=12 if ok else 6,
                     mfc=BLUE if ok else "#FFFFFF", mec=STROKE if ok else MUTED,
                     mew=1.4 if ok else 0.9, zorder=3)
    axR.set_xticks(range(6)); axR.set_yticks(range(6))
    axR.set_xlabel(T("$A_4$ mod 6", "$A_4$ mod 6"))
    axR.set_ylabel(T("$8D$ mod 6", "$8D$ mod 6"))
    axR.set_title(T("One residue class in six", "Una clase de residuos entre seis"),
                  loc="left", fontsize=10, color=INK, fontweight="bold", pad=9)
    # raw strings: "\e" is not a valid Python escape.  It happened to work because Python left
    # the backslash alone, but it emitted a SyntaxWarning on first compile and NOT on reruns
    # from cached bytecode -- so the archived output and a fresh run disagreed, and the
    # reproducibility gate reported this file as broken for a reason that was not about figures.
    # las dos mitades del "uno de cada seis" NO son del mismo rango: la imparidad depende de
    # la semilla y la congruencia mod 3 vale con las dos.  El pie lo distingue desde el
    # 24-ago-2026; el texto dibujado las juntaba, y es donde el lector mira primero.
    axR.text(2.5, -0.42, T(r"$8D$ odd (seed of [2])  and  $8D+A_4\equiv0$ (mod 3, either seed)",
                           r"$8D$ impar (semilla de [2])  y  $8D+A_4\equiv0$ (mod 3, las dos)"),
             ha="center", fontsize=7.6, color=STROKE)
    axR.set_aspect("equal")
    tidy(axR)

    fig.savefig(os.path.join(OUT, "fig_cone%s.pdf" % SUF), bbox_inches="tight")
    plt.close(fig)


# BOTH EDITIONS IN ONE INVOCATION.  The archive was built by running this twice, once with
# --es, so a single run reproduced only half of it and check_reproduces.py reported the file as
# differing.  A script's archived output has to be what the script prints, once.
for _es in (False, True):
    ES = _es
    T = (lambda en, es: es) if _es else (lambda en, es: en)
    SUF = "_es" if _es else ""
    for f in (fig_ladder, fig_ceiling, fig_parity, fig_surface, fig_cone):
        f()
        print("drew %s%s.pdf" % (f.__name__, SUF), flush=True)
    print("figures written to %s" % OUT)
