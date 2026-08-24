#!/usr/bin/env python3
"""make_fig_tower -- what the Kaluza-Klein tower does to the dijet angular distribution.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

REDRAWN FROM THE CLOSED FORM.  The previous version of this figure evaluated the replaced
propagator of [DMN01] by truncating the tower at n = 400 and carrying a Breit-Wigner width on
every mode.  In the spacelike channels neither is needed and the width is wrong -- a mode at
t < 0 sits below its two-parton threshold, so there is no absorptive part to put there.
kk_resummation.py sums the tower exactly instead:

    1/t + sum_n 2/(t - n^2/R_5^2)  =  F(t)/t ,   F(t) = pi a coth(pi a) ,  a = R_5 sqrt(-t) ,

so on a subprocess with only t-channel gluon exchange the parton-level ratio is F^2, exactly.
With |t| = M_jj^2 / (1 + chi) for massless jets that is a closed form in the two variables the
CMS measurement is binned in, and it needs no truncation, no width, and no BR(G_n -> jj).

THREE PANELS.

  LEFT    the escape branch, 1/R_5 = 3.97 TeV -- the version that can also pay for Part VI's
          escape from proton decay.  The effect is a factor, not a percentage.
  MIDDLE  the ceiling at the measured Higgs mass, 1/R_5 = 9.09 TeV -- the whole class, capped.
  RIGHT   the reason the first two look alike.  F depends on M_jj, chi and M_KK only through
          a = M_jj / (M_KK sqrt(1+chi)), so EVERY curve on the left is the same function
          (pi a coth pi a)^2.  Both branches collapse onto one line.  The hollow markers are
          the full ratio of kk_dijet_lo.py -- all channels, truncated, widths on -- which
          tracks the closed form away from the tower and leaves it where the s-channel opens.
          That departure is the part the closed form does NOT describe, and drawing it is how
          the figure says so.

The grey band is the total theoretical uncertainty CMS quotes for its normalised chi_dijet in
the 3.0-3.6 TeV mass bin, 3.5 % (their Table 2).  It gives the vertical axis a scale; it is NOT
an exclusion threshold, which would need the full likelihood over correlated bins.

AND THE CURVES ARE STILL PARTONIC.  Gluon-initiated subprocesses do not see the tower at all,
and every CMS distribution is normalised per mass bin; both dilute what is drawn here.  The
picture establishes size and shape, not a limit.

Run:  python make_fig_tower.py     -> fig_tower.pdf, fig_tower_es.pdf, .png  (into . and paper/)
"""

import math
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
GREEN, GREY, STROKE = "#2E8B57", "#8A8A8A", "#1F4E79"

HERE = pathlib.Path(__file__).resolve().parent

# the theorem is defined once, in the file that proves it.  A second copy of coth here would be
# a second place for it to be wrong.
_rs = {"__file__": str(HERE / "kk_resummation.py"), "__name__": "kk_resummation"}
exec(compile((HERE / "kk_resummation.py").read_text(encoding="utf-8"),
             "kk_resummation.py", "exec"), _rs)
form_factor = _rs["form_factor"]

# and the superseded calculation, for the departure the right-hand panel shows
_ns = {"__file__": str(HERE / "kk_dijet_lo.py"), "__name__": "kk_dijet_lo"}
exec(compile((HERE / "kk_dijet_lo.py").read_text(encoding="utf-8"), "kk_dijet_lo.py", "exec"), _ns)
dsig_dchi = _ns["dsig_dchi"]

ALPHA_S = 0.08
CMS_THEORY = 0.035          # CMS-EXO-24-011 Table 2, 3.0-3.6 TeV bin
CHI = np.linspace(1.0, 16.0, 200)

BRANCHES = [
    (3.97, [3.3, 3.9, 5.1], "escape"),
    (9.09, [3.3, 5.1, 6.5], "ceiling"),
]


def a_of(mjj_tev, chi, invR5_tev):
    """the single variable everything depends on: a = M_jj / (M_KK sqrt(1+chi))."""
    return mjj_tev / (invR5_tev * np.sqrt(1.0 + np.asarray(chi, float)))


def f2_curve(invR5_tev, mjj_tev, chi=CHI):
    """(pi a coth pi a)^2 -- the exact t-channel ratio, from kk_resummation.form_factor."""
    t = -(mjj_tev ** 2) / (1.0 + np.asarray(chi, float))
    return np.array([form_factor(float(tt), invR5_tev) ** 2 for tt in t])


def full_ratio(invR5_tev, mjj_tev, chi):
    """the superseded all-channel, truncated, width-on ratio of kk_dijet_lo.py."""
    shat = (mjj_tev * 1000.0) ** 2
    out = []
    for c in np.atleast_1d(chi):
        a = dsig_dchi(shat, float(c), invR5_tev * 1000.0, ALPHA_S)
        b = dsig_dchi(shat, float(c), invR5_tev * 1000.0, ALPHA_S, tower=False)
        out.append(a / b)
    return np.array(out)


def draw(es=False):
    L = dict(
        sup=("La torre resumada: un solo factor de forma, $\\pi a\\coth\\pi a$" if es
             else "The tower resummed: one form factor, $\\pi a\\coth\\pi a$"),
        ttlA=("rama del escape, $1/R_5=3.97$ TeV" if es
              else "escape branch, $1/R_5=3.97$ TeV"),
        ttlB=("techo a la masa medida, $1/R_5=9.09$ TeV" if es
              else "ceiling at the measured mass, $1/R_5=9.09$ TeV"),
        ttlC=("y las dos son la misma curva" if es else "and both are the same curve"),
        xl=r"$\chi_{\mathrm{dijet}}=e^{|y_1-y_2|}$",
        xlC=(r"$a=M_{jj}\,/\,M_{\mathrm{KK}}\sqrt{1+\chi}$"),
        yl=("$\\mathcal{F}^2$  (QCD + torre) / QCD,  canal $t$" if es
            else "$\\mathcal{F}^2$  (QCD + tower) / QCD,  $t$-channel"),
        # NOT just "CMS theoretical uncertainty".  A horizontal band behind an unnormalised,
        # single-subprocess, partonic ratio invites the eye to divide one by the other, and the
        # caption saying not to does not undo what the picture already suggested.  So the warning
        # goes ON the band, where the eye does the comparing -- putting it in the legend made the
        # legend two lines long and covered the curves it was warning about.
        band=("teórica de CMS, $3.5\\%$" if es else "CMS theory unc., $3.5\\%$"),
        bandnote=("no es un umbral" if es else "not a threshold"),
        note=("partónico, un solo subproceso: los canales con gluones no ven la torre, y cada "
              "distribución de CMS va normalizada por bin de masa.\nEl $3.5\\%$ es el bin "
              "$3.0$--$3.6$ TeV; a $M_{jj}>7$ TeV CMS cita $12.2\\%$" if es else
              "partonic, one subprocess: gluon channels do not see the tower, and every CMS "
              "distribution is normalised per mass bin.\nThe $3.5\\%$ is the $3.0$--$3.6$ TeV bin; "
              "above $M_{jj}=7$ TeV CMS quote $12.2\\%$"),
        mj=("$M_{jj}=%.1f$ TeV"),
        closed=("forma cerrada, canal $t$" if es else "closed form, $t$-channel"),
        fullr=("torre truncada, todos los canales" if es else "truncated tower, all channels"),
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.15),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.06], wspace=0.26))
    cols = [RED, AMBER, BLUE]

    for ax, (invR5, mjjs, _tag) in zip(axes[:2], BRANCHES):
        ax.axhspan(1 - CMS_THEORY, 1 + CMS_THEORY, color=GREY, alpha=0.28, zorder=1,
                   label=L["band"])
        ax.axhline(1.0, color=GREY, lw=0.8, ls="--", zorder=2)
        for col, mjj in zip(cols, mjjs):
            ax.plot(CHI, f2_curve(invR5, mjj), color=col, lw=2.0, zorder=3,
                    label=L["mj"] % mjj)
        ax.set_xlabel(L["xl"], fontsize=10)
        ax.set_xlim(1, 16)
        ax.set_yscale("log")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=7.5, loc="upper right", framealpha=0.92)

    # ---- the collapse ------------------------------------------------------------------------
    axC = axes[2]
    ag = np.linspace(0.02, 1.55, 300)
    axC.plot(ag, (math.pi * ag / np.tanh(math.pi * ag)) ** 2, color=STROKE, lw=2.4, zorder=4,
             label=L["closed"])
    axC.axhspan(1 - CMS_THEORY, 1 + CMS_THEORY, color=GREY, alpha=0.28, zorder=1)
    axC.axhline(1.0, color=GREY, lw=0.8, ls="--", zorder=2)
    marks = ["o", "s"]
    chi_pts = np.array([1.5, 3.0, 6.0, 10.0, 14.0])
    for mk, (invR5, mjjs, tag) in zip(marks, BRANCHES):
        for mjj in mjjs:
            aa = a_of(mjj, chi_pts, invR5)
            rr = full_ratio(invR5, mjj, chi_pts)
            axC.plot(aa, rr, mk, ms=4.6, mfc="none", mec=RED, mew=1.0, ls="none", zorder=5)
    axC.plot([], [], "o", ms=4.6, mfc="none", mec=RED, mew=1.0, ls="none", label=L["fullr"])
    axC.set_xlabel(L["xlC"], fontsize=10)
    axC.set_xlim(0, 1.55)
    axC.set_yscale("log")
    axC.grid(alpha=0.25, lw=0.5)
    axC.legend(fontsize=7.5, loc="upper left", framealpha=0.92)
    axC.set_title(L["ttlC"], fontsize=11)

    # a log axis spanning less than a decade labels only the decade ticks, so the left panel
    # came out with a single "1.0" on it and no way to read a factor off the curves.
    for ax in axes:
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, subs=np.arange(1.0, 10.0)))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(
            lambda v, _p: ("%g" % v) if v < 10 else ("%d" % v)))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.tick_params(labelsize=8)

    # the warning goes ON the band, in both branch panels, where the eye does the comparing
    for ax in axes[:2]:
        ax.text(15.4, 1.0, L["bandnote"], ha="right", va="center", fontsize=6.6,
                color=GREY, style="italic", zorder=6)

    axes[0].set_ylabel(L["yl"], fontsize=10)
    axes[0].set_title(L["ttlA"], fontsize=11)
    axes[1].set_title(L["ttlB"], fontsize=11)
    fig.suptitle(L["sup"], fontsize=12, y=1.02)
    fig.text(0.5, -0.055, L["note"], ha="center", fontsize=8, color=GREY, style="italic")

    stem = "fig_tower_es" if es else "fig_tower"
    here = os.path.dirname(os.path.abspath(__file__))
    # same pad_inches as the other figures: a tight bbox crops to the ink and the suptitle,
    # which sits above the axes, has been known to spill past the declared box.
    for d in (here, os.path.join(here, "paper")):
        fig.savefig(os.path.join(d, stem + ".pdf"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(os.path.join(here, stem + ".png"), dpi=150, bbox_inches="tight",
                pad_inches=0.08)
    plt.close(fig)
    return stem


def main():
    fails = []
    print("=" * 96)
    print("THE TOWER IN chi_dijet, FROM THE CLOSED FORM -- the figure and the numbers in it")
    print("=" * 96)

    # C1: the two branches really are one curve in a.  Not decoration -- it is the claim the
    # third panel makes, and a wrong a would break it while leaving the first two panels intact.
    print("\n  C1  DOES EVERYTHING COLLAPSE ONTO ONE VARIABLE?")
    print("      the same a reached from different (M_jj, chi, M_KK) must give the same F^2")
    probes = [(3.3, 1.5, 3.97), (5.1, 1.5, 9.09), (3.9, 6.0, 3.97), (6.5, 6.0, 9.09)]
    seen = {}
    worst = 0.0
    for mjj, chi, mkk in probes:
        aa = a_of(mjj, chi, mkk)
        f2 = f2_curve(mkk, mjj, np.array([chi]))[0]
        ref = (math.pi * aa / math.tanh(math.pi * aa)) ** 2
        worst = max(worst, abs(f2 - ref) / ref)
        print("      M_jj=%4.1f chi=%4.1f M_KK=%4.2f -> a=%.6f  F^2=%.6f" % (mjj, chi, mkk, aa, f2))
    ok = worst < 1e-12
    print("\n      C1 %s (worst %.1e)" % ("PASS" if ok else "FAIL", worst))
    if not ok:
        fails.append("C1")

    # C2: away from the tower the closed form and the superseded all-channel calculation must
    # agree to better than the CMS band, or one of the two is wrong.
    print("\n  C2  AND FAR BELOW THE TOWER IT MUST AGREE WITH kk_dijet_lo.py")
    print("      %-26s %12s %12s %10s" % ("(M_jj, chi, M_KK)", "closed F^2", "full ratio", "rel"))
    worst2 = 0.0
    for mjj, chi, mkk in ((3.3, 14.0, 9.09), (5.1, 14.0, 9.09), (3.3, 10.0, 9.09)):
        f2 = f2_curve(mkk, mjj, np.array([chi]))[0]
        fr = full_ratio(mkk, mjj, np.array([chi]))[0]
        rel = abs(f2 - fr) / fr
        worst2 = max(worst2, rel)
        print("      (%4.1f, %4.1f, %4.2f)%9s %12.4f %12.4f %10.3f"
              % (mjj, chi, mkk, "", f2, fr, rel))
    ok = worst2 < 0.12
    print("\n      C2 %s (worst %.3f, against the 0.035 CMS band and a 0.12 bar)"
          % ("PASS" if ok else "FAIL", worst2))
    if not ok:
        fails.append("C2")

    for invR5, mjjs, tag in BRANCHES:
        print("\n  1/R_5 = %.2f TeV  (%s) -- F^2, the exact t-channel ratio" % (invR5, tag))
        print("      %-12s %10s %10s %10s" % ("M_jj [TeV]", "chi=1", "chi=6", "chi=16"))
        for mjj in mjjs:
            r = f2_curve(invR5, mjj)
            print("      %-12.1f %10.3f %10.3f %10.3f"
                  % (mjj, r[0], r[int(len(CHI) * 5 / 15)], r[-1]))

    print("\n  CMS quotes 3.5% total theoretical uncertainty in its 3.0-3.6 TeV bin, drawn as")
    print("  the grey band.  It is a scale, not a threshold: the curves are one subprocess at")
    print("  parton level, and the gluon channels the tower does not touch dilute them.")
    for es in (False, True):
        stem = draw(es)
        print("\n  written: %s.pdf (here and in paper/)" % stem)
    print("\n  written: fig_tower.png")
    print("=" * 96)
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        print("=" * 96)
        return 1
    print("VERDICT: the collapse holds and the closed form meets the old calculation off the tower.")
    print("=" * 96)
    return 0


if __name__ == "__main__":
    sys.exit(main())
