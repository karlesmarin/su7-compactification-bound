#!/usr/bin/env python3
"""make_fig_lambert.py -- the fold. The two Lambert branches as one surface, and what sits on it.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS FIGURE EXISTS.  Sections 5 and 6 of the paper carry its central result -- the hierarchy
in closed form, eq. (lambert), and the theorem that says which Lambert branch is the minimum,
eq. (branchmu) -- and they were the only sections with no picture at all.  The words "branch
point" and "fold" appear there and nothing draws them.

THE SIMPLIFICATION THAT MAKES IT DRAWABLE, and it is exact.  The closed form is

    x^2 = -d / W(z),      z = -d e^{-b},      b = (4G - A_4)/(2 A_4),   d = 12 zeta(3) D / A_4,

and two real branches exist exactly when L := b - 1 - ln d >= 0.  Change coordinates from (b, d)
to (L, d).  Then b = L + 1 + ln d, so

    z = -d e^{-b} = -d e^{-(L+1)}/d = -e^{-(1+L)} ,

which depends on L ALONE.  The whole two-parameter family collapses to a product:

    x^2 = d * f(L),        f(L) = -1 / W(-e^{-(1+L)}) ,

one factor per axis.  The fold, where the two branches meet, is W = -1, hence f = 1, hence the
STRAIGHT edge L = 0 -- and by eq. (branchmu), mu = -(A_4/6)(W+1), the whole of that edge is
mu = 0, that is m_h = 0.  So the picture is honest about the two things the section claims: the
physical branch is one sheet of a fold, and the Higgs mass is the height above its edge.

WHAT IS PLOTTED.  1/R_5 = 2 pi m_W / x, the quantity the paper is about, over (L, d) with d on a
log axis so that the five published rows and the certificate's own witness -- which differ by a
factor of twenty in d -- fit on one surface.  The upper sheet is W_{-1}, where mu >= 0 and the
stationary point is a minimum; the lower sheet is W_0, where mu < 0 and it is a maximum.  The
lower sheet is drawn muted because nothing physical lives on it.

CONTROLS.  C1 reads alpha_min back off the surface at each published row and compares it with the
archived closed form; C2 does the same for the certificate witness against ceiling_ilp.json; C3
checks that the fold really is mu = 0 by evaluating f at L = 0.  A figure that cannot fail its own
data is decoration.  [[a-figure-is-an-unaudited-document]]

Run:  python make_fig_lambert.py        (add --es for the Spanish edition)
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
from scipy.special import lambertw

HERE = os.path.dirname(os.path.abspath(__file__))
ES = "--es" in sys.argv
T = lambda en, es: es if ES else en
SUF = "_es" if ES else ""
OUT = os.path.join(HERE, "paper")
J = lambda n: json.load(open(os.path.join(HERE, "outputs", n)))

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
STROKE, INK, MUTED = "#1F4E79", "#1F2933", "#6B7280"
GRID, SURF = "#E5E7EB", "#FCFCFB"
SEQ = LinearSegmentedColormap.from_list("seq", ["#EAF2FA", "#9CC4E4", BLUE, STROKE, "#12314C"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.facecolor": SURF, "figure.facecolor": SURF,
})

MW = 80.4
Z3 = 1.2020569031595942854


def f_of(L, branch):
    """f(L) = -1/W(-e^{-(1+L)}).  branch -1 is the minimum, 0 the maximum.

    At L = 0 exactly, z = -1/e and both branches are W = -1, so f = 1 -- but scipy returns nan
    there rather than the limit, which is what a branch point looks like to a library.  The value
    is put in by hand; C3 checks that lambertw itself walks up to it from both sides, so the
    hand-placed point is verified rather than assumed."""
    L = np.asarray(L, float)
    z = -np.exp(-(1.0 + L))
    out = np.ones_like(L)
    m = L > 0
    if np.any(m):
        out[m] = -1.0 / lambertw(z[m], branch).real
    return out if out.ndim else float(out)


def invR(L, d, branch):
    """1/R_5 in TeV on the given sheet."""
    x = np.sqrt(np.asarray(d, float) * f_of(L, branch))
    return 2.0 * math.pi * MW / x / 1000.0


def main():
    LC = J("lambert_criticality.json")
    rows = LC["rows"]
    CE = J("ceiling_ilp.json")
    CF = J("certify_212_215.json")

    # the certificate's own attained vertex, in the same coordinates
    A4c, kc = 212.0, 1.0
    Gc = CF["vertex212"]["G_min"]
    dc = 12.0 * Z3 * (kc / 8.0) / A4c
    bc = (4.0 * Gc - A4c) / (2.0 * A4c)
    Lc = bc - 1.0 - math.log(dc)

    fails = []

    # ---- controls, before anything is drawn -------------------------------------------------
    worst = 0.0
    for r in rows:
        x = math.sqrt(r["d"] * float(f_of(r["L"], -1)))
        worst = max(worst, abs(x / math.pi - r["amin"]) / r["amin"])
    ok1 = worst < 1e-9
    print("  C1  the surface returns each published row's alpha_min ... %s (worst %.1e)"
          % ("PASS" if ok1 else "FAIL", worst))
    if not ok1:
        fails.append("C1")

    # C2 is NOT "does the surface return ceiling_GeV".  It does not, and it should not: the
    # certificate evaluates (II) at the TOP of the Higgs window, while the surface returns the
    # content's own stationary point at its own Higgs mass.  The two differ by 0.7 % and comparing
    # them would have been a control that fails for the right answer.  What is checkable without
    # circularity is that the x the surface gives at that vertex solves the closed form itself:
    #     x^2 [ 4G - A_4 (4 ln x + 1) ] = 24 zeta(3) D .
    xc = math.sqrt(dc * float(f_of(Lc, -1)))
    lhs = xc ** 2 * (4.0 * Gc - A4c * (4.0 * math.log(xc) + 1.0))
    rhs = 24.0 * Z3 * (kc / 8.0)
    res = abs(lhs - rhs) / abs(rhs)
    ok2 = res < 1e-12
    print("  C2  at the certificate's vertex the surface solves eq.(amin) %s (residual %.1e)"
          % ("PASS" if ok2 else "FAIL", res))
    if not ok2:
        fails.append("C2")
    print("      (its height there is %.0f GeV -- its own stationary point, not the %.0f GeV the"
          % (float(invR(Lc, dc, -1)) * 1000.0, CE["ceiling_GeV"]))
    print("       certificate quotes, which is read at the top of the Higgs window instead)")

    # C3: the two sheets really do meet at L = 0, approached from lambertw and from opposite
    # sides.  A control that just read back the hand-placed f(0) = 1 could not fail.
    eps = 1e-8
    fm, fp = float(f_of(eps, -1)), float(f_of(eps, 0))
    ok3 = abs(fm - 1) < 1e-3 and abs(fp - 1) < 1e-3 and fm < 1.0 < fp
    print("  C3  the two sheets meet at L=0, from below and above ..... %s (%.8f, %.8f)"
          % ("PASS" if ok3 else "FAIL", fm, fp))
    if not ok3:
        fails.append("C3")

    if fails:
        print("  REFUSING TO DRAW: %s" % ", ".join(fails))
        return 1

    # ---- the surface -------------------------------------------------------------------------
    Ls = np.array([r["L"] for r in rows] + [Lc])
    ds = np.array([r["d"] for r in rows] + [dc])
    LO, LH = 0.0, max(2.35, Ls.max() * 1.12)
    # the lower edge in d is set by the certificate's vertex rather than by taste: going much
    # below it sends the surface past 16 TeV, which is off the top of anything the paper bounds
    # and makes the whole picture about a corner nothing lives in.
    DO, DH = dc * 0.80, ds.max() * 1.30

    gl = np.linspace(LO + 1e-6, LH, 190)
    gd = np.exp(np.linspace(math.log(DO), math.log(DH), 190))
    GL, GD = np.meshgrid(gl, gd)
    Zlo = invR(GL, GD, 0)          # W_0 : mu < 0, a maximum
    Zhi = invR(GL, GD, -1)         # W_-1: mu >= 0, the vacuum
    LOGD = np.log10(GD)

    # zoom, not a bigger figure: a 3-D axes keeps a wide margin of its own, and at the size the
    # paper prints this at that margin was most of the panel.
    fig = plt.figure(figsize=(10.6, 5.7))
    ax = fig.add_axes([-0.020, 0.050, 0.765, 0.895], projection="3d", computed_zorder=False)
    ax.set_box_aspect((1.42, 1.0, 0.80), zoom=1.10)
    ax.view_init(elev=21, azim=-124)

    # THE UNPHYSICAL SHEET.  It is drawn as a wireframe, not as a filled surface: filled and grey
    # it read as a floor -- a shadow of the blue sheet rather than the other side of the fold.
    # Open, the reader can see the blue sheet through it and see the two leave one edge.
    ax.plot_wireframe(GL, LOGD, Zlo, rstride=14, cstride=14, color="#9AA3AD", lw=0.55,
                      alpha=0.75, zorder=1)
    norm = plt.Normalize(np.nanmin(Zhi), np.nanpercentile(Zhi, 99))
    ax.plot_surface(GL, LOGD, Zhi, facecolors=SEQ(norm(np.clip(Zhi, None, np.nanpercentile(Zhi, 99)))),
                    rstride=2, cstride=2, linewidth=0, antialiased=True, shade=False,
                    alpha=0.985, zorder=3)
    # contours ON the sheet, so the height is readable and not only shaded
    for lev in (3, 4, 5, 6, 8, 10):
        cs = plt.contour(GL, LOGD, Zhi, levels=[lev])
        for pth in cs.allsegs[0]:
            if len(pth) > 3:
                ax.plot(pth[:, 0], pth[:, 1], np.full(len(pth), lev), color="white",
                        lw=0.7, alpha=0.55, zorder=4)
        plt.close(cs.axes.figure) if cs.axes.figure is not fig else None

    # the fold itself: the straight edge L = 0, where the two sheets meet and m_h = 0
    fold = invR(np.zeros_like(gd), gd, -1)
    ax.plot(np.zeros_like(gd), np.log10(gd), fold, color=RED, lw=3.0, zorder=8,
            solid_capstyle="round")

    # the published rows, with stems dropped onto the unphysical sheet below them
    for r in rows:
        v = float(invR(r["L"], r["d"], -1))
        v0 = float(invR(r["L"], r["d"], 0))
        ax.plot([r["L"]] * 2, [math.log10(r["d"])] * 2, [v0, v],
                color=STROKE, lw=0.8, alpha=0.5, zorder=7)
        ax.scatter([r["L"]], [math.log10(r["d"])], [v], s=46, c=AMBER,
                   edgecolors=STROKE, linewidths=0.9, depthshade=False, zorder=9)
    vc, vc0 = float(invR(Lc, dc, -1)), float(invR(Lc, dc, 0))
    ax.plot([Lc] * 2, [math.log10(dc)] * 2, [vc0, vc], color=RED, lw=1.2, alpha=0.85, zorder=7)
    ax.scatter([Lc], [math.log10(dc)], [vc], s=150, marker="*", c=RED,
               edgecolors="white", linewidths=1.0, depthshade=False, zorder=10)

    ax.set_xlim(LO, LH)
    ax.set_ylim(math.log10(DO), math.log10(DH))
    ax.set_zlim(0, float(np.nanpercentile(Zhi, 99)) * 1.06)

    # labels in FIGURE coordinates, not data coordinates: a 3-D text anchor moves with the camera
    # and the previous version's captions landed on top of the markers they named.
    # the ratio is COMPUTED, not asserted: an earlier draft wrote "twelve times" from memory and
    # the rows give ten.  A number in a label is a claim like any other.
    gap = float(min(r["d"] for r in rows) / dc)
    box = dict(fc=SURF, ec="none", alpha=0.88, pad=1.8)
    fig.text(0.505, 0.855, T("the star: $(A_4,8D)=(212,1)$, where the ceiling is attained",
                             "la estrella: $(A_4,8D)=(212,1)$, donde se alcanza el techo"),
             color=RED, fontsize=8.6, fontweight="bold", ha="left", bbox=box, zorder=20)
    fig.text(0.505, 0.822,
             T("and $d$ there is %.0f times smaller than on their smallest row" % gap,
               "y ahí $d$ es %.0f veces menor que en su fila más pequeña" % gap),
             color=MUTED, fontsize=7.8, ha="left", bbox=box, zorder=20)
    fig.text(0.255, 0.375, T("their five published rows", "sus cinco filas publicadas"),
             color=STROKE, fontsize=8.6, fontweight="bold", ha="left", bbox=box, zorder=20)
    fig.text(0.028, 0.400, T("the fold, $L=0$:\n$W=-1$,  $m_h=0$",
                             "el pliegue, $L=0$:\n$W=-1$,  $m_h=0$"),
             color=RED, fontsize=9.0, fontweight="bold", style="italic", ha="left",
             bbox=box, zorder=20)

    ax.set_xlabel(T("$L=b-1-\\ln d$   (height above the fold)",
                    "$L=b-1-\\ln d$   (altura sobre el pliegue)"), labelpad=4, fontsize=8.6)
    ax.set_ylabel(T("$\\log_{10} d$,   $d=12\\zeta(3)D/A_4$", "$\\log_{10} d$,   $d=12\\zeta(3)D/A_4$"),
                  labelpad=4, fontsize=8.6)
    ax.set_zlabel(T("$1/R_5$  (TeV)", "$1/R_5$  (TeV)"), labelpad=-2, fontsize=8.6)
    ax.tick_params(labelsize=7.4, pad=0.5)
    fig.text(0.012, 0.955, T("One surface, and the physical branch is the upper sheet of a fold",
                             "Una superficie, y la rama física es la hoja de arriba de un pliegue"),
             fontsize=10.8, color=INK, fontweight="bold")

    handles = [Line2D([], [], color=STROKE, lw=6, alpha=0.85,
                      label=T("$W_{-1}$: $\\mu\\geq0$, the vacuum",
                              "$W_{-1}$: $\\mu\\geq0$, el vacío")),
               Line2D([], [], color="#9AA3AD", lw=1.2,
                      label=T("$W_{0}$: $\\mu<0$, a maximum", "$W_{0}$: $\\mu<0$, un máximo")),
               Line2D([], [], color=RED, lw=2.6, label=T("the fold, $m_h=0$", "el pliegue, $m_h=0$"))]
    fig.legend(handles=handles, frameon=False, fontsize=8.6, loc="upper left",
               bbox_to_anchor=(0.055, 0.925))

    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_facecolor("#F4F6F8")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a._axinfo["grid"]["color"] = GRID
        a._axinfo["grid"]["linewidth"] = 0.5

    # ---- the panel at the right: the theorem, as the straight line it is ----------------------
    ins = fig.add_axes([0.795, 0.155, 0.185, 0.42])
    ww = np.linspace(-5.2, -1.0, 200)
    for r in rows:
        ins.plot(ww, -(r["A4"] / 6.0) * (ww + 1.0), color=BLUE, lw=0.9, alpha=0.5)
        ins.plot([r["W"]], [r["mu"]], "o", ms=5, mfc=AMBER, mec=STROKE, mew=0.9, zorder=5)
    ins.axvline(-1.0, color=RED, lw=1.4)
    ins.axhline(0.0, color=MUTED, lw=0.7)
    ins.text(-1.06, ins.get_ylim()[1] * 0.62, T("$m_h=0$", "$m_h=0$"), color=RED,
             fontsize=7.6, ha="right", fontweight="bold")
    ins.set_xlabel("$W_{-1}$", fontsize=7.6, labelpad=1)
    ins.set_ylabel("$\\mu=m_h^2/(K\\pi^2)^2$", fontsize=7.2, labelpad=1)
    ins.set_title(T("$\\mu=-\\frac{A_4}{6}(W+1)$, exactly",
                    "$\\mu=-\\frac{A_4}{6}(W+1)$, exacta"), fontsize=8.2, color=INK, pad=3)
    ins.tick_params(labelsize=6.8, length=2.2, width=0.6)
    for s in ins.spines.values():
        s.set_linewidth(0.7)
    ins.spines["top"].set_visible(False)
    ins.spines["right"].set_visible(False)
    ins.grid(color=GRID, lw=0.5)
    ins.set_axisbelow(True)

    # no provenance strip along the bottom: the caption carries it, and here it only crowded the
    # x-axis label.  One statement, one place.

    path = os.path.join(OUT, "fig_lambert%s.pdf" % SUF)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("  drew %s" % os.path.basename(path))
    return 0


def both():
    """Draw both editions in ONE invocation, as make_fig_tower.py does.

    Archiving two separate runs with '>' then '>>' left a file the reproducibility gate could
    never match: one run produces half of it.  A script's archived output has to be what the
    script prints, once."""
    global ES, T, SUF
    rc = 0
    for es in (False, True):
        ES = es
        T = (lambda en, esp: esp) if es else (lambda en, esp: en)
        SUF = "_es" if es else ""
        rc |= main()
    return rc


if __name__ == "__main__":
    sys.exit(both())
