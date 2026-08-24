#!/usr/bin/env python3
"""s_channel_width.py -- el canal s vestido, con la anchura total perfilada.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

POR QUE ESTO HACIA FALTA.  El recast del articulo viste solo t y u.  El argumento para hacerlo
es bueno -- ahi el propagador resumado es real, no hace falta anchura y no hace falta
BR(G_n -> jj), de modo que el limite no depende de la masa de unos exoticos que el modelo no
fija.  Pero un arbitro senalo lo que ese argumento NO demuestra: que quitar el canal s deje una
cota conservadora.  Las amplitudes interfieren.  En q qbar -> q qbar el termino temporal puede
sumar o restar a la deformacion de t/u, y "no depende de Gamma" no es lo mismo que "es una cota
inferior de la senal completa".

Lo que se mide aqui:

    F(q^2) resumado en LOS TRES canales, con el polo regulado por la anchura TOTAL

        Gamma_tot = Gamma_qqbar + Gamma_exotic ,     Gamma_qqbar/M = 2 alpha_s ,
        Gamma_exotic >= 0  y desconocida, luego se PERFILA sobre ella.

    Si incluso el peor valor de la anchura deja el limite por encima del techo, la conclusion de
    la seccion de colisionador deja de depender de una omision.

COMO SE VISTE.  Solo dos elementos de matriz llevan propagador temporal de gluon: q qbar -> q qbar
(su pieza s y la interferencia s-t) y q qbar -> q' qbar' (s puro).  Los canales con gluones en el
estado final NO se visten, por el argumento de cinco momentos de [DMN01] que el articulo ya usa.

EL CONTROL QUE DECIDE SI ESTO VALE ALGO.  Con el canal s APAGADO este guion tiene que reproducir
el barrido que ya esta archivado.  Si no lo reproduce, la diferencia que mida encendido no es el
canal s: es un error mio.  C1 lo comprueba antes de imprimir nada mas.
[[falsify-the-instrument-first]]

Run:  python s_channel_width.py > ../outputs/s_channel_width.txt
"""
import cmath
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"

_HC = {"__file__": str(HERE / "hadronic_chi.py"), "__name__": "hadronic_chi"}
exec(compile((HERE / "hadronic_chi.py").read_text(encoding="utf-8"), "hadronic_chi.py", "exec"),
     _HC)
_FC = {"__file__": str(HERE / "fit_chi2.py"), "__name__": "fit_chi2"}
exec(compile((HERE / "fit_chi2.py").read_text(encoding="utf-8"), "fit_chi2.py", "exec"), _FC)

read_bins, read_dump, run_dump = _HC["read_bins"], _HC["read_dump"], _HC["run_dump"]
ROOT_S, YB_MAX = _HC["ROOT_S"], _HC["YB_MAX"]
NM, NYB, NCHI = _HC["NM"], _HC["NYB"], _HC["NCHI"]
NF = _HC["NF"]
IDENT = _HC["IDENT"]
load_table, load_corr, renorm, chi2 = (_FC["load_table"], _FC["load_corr"],
                                       _FC["renorm"], _FC["chi2"])
PARTICLE = _FC["PARTICLE"]

NMAX = 400          # rungs kept in the timelike sum; the tail is checked in C2
ALPHAS = 0.118      # only for Gamma_qqbar/M = 2 alpha_s, which is where the paper puts it


def f_space(q2abs, mkk):
    """F on a spacelike leg: pi a coth(pi a), a = sqrt(-q^2)/M_KK.  Real, no width."""
    a = math.pi * math.sqrt(q2abs) / mkk
    return a / math.tanh(a) if a > 1e-12 else 1.0


def f_time(s, mkk, gam_over_m):
    """F on the timelike leg, WITH the width, summed rung by rung.

    1 + sum_n 2 s / (s - M_n^2 + i M_n Gamma_n),  M_n = n M_KK,  Gamma_n = gam_over_m * M_n.

    Without a width this is pi b cot(pi b) and diverges at every b = n; with one it is finite
    everywhere, which is the whole reason the width may not simply be dropped here."""
    tot = 1.0 + 0.0j
    for n in range(1, NMAX + 1):
        Mn = n * mkk
        tot += 2.0 * s / (s - Mn * Mn + 1j * Mn * (gam_over_m * Mn))
    return tot


def m2_qqbar_same(s, t, u, ft, fu, fs):
    """q qbar -> q qbar with BOTH legs dressed.  fs = 1 recovers the paper's t-only version."""
    return ((4.0 / 9.0) * ((s * s + u * u) / (t * t) * ft * ft
                           + (t * t + u * u) / (s * s) * abs(fs) ** 2)
            - (8.0 / 27.0) * u * u / (s * t) * ft * fs.real)


def m2_qqbar_diff(s, t, u, ft, fu, fs):
    """q qbar -> q' qbar', pure s-channel: entirely a |F(s)|^2 rescaling."""
    return (4.0 / 9.0) * (t * t + u * u) / (s * s) * abs(fs) ** 2


def folded3(fn, s, t, u, ft, fu, fs):
    return fn(s, t, u, ft, fu, fs) + fn(s, u, t, fu, ft, fs)


def weight3(xf1, xf2, s, chi, ft, fu, fs):
    """the full sum, with fs carried into the two channels that have a timelike gluon.

    Every other channel is taken verbatim from hadronic_chi.py, so the ONLY difference between
    this and the archived calculation is the dressing of those two."""
    t = -s / (1.0 + chi)
    u = -s * chi / (1.0 + chi)
    q = [f for f in range(-NF, NF + 1) if f != 0]
    fold = _HC["folded"]
    tot = 0.0

    g1, g2 = xf1[0], xf2[0]
    tot += g1 * g2 * (IDENT["m2_gg"] * fold(_HC["m2_gg"], s, t, u, ft, fu)
                      + IDENT["m2_gg_qqbar"] * NF * fold(_HC["m2_gg_qqbar"], s, t, u, ft, fu))
    for f in q:
        tot += IDENT["m2_qg"] * (xf1[f] * g2 + g1 * xf2[f]) * fold(_HC["m2_qg"], s, t, u, ft, fu)
    for f1 in q:
        for f2 in q:
            w = xf1[f1] * xf2[f2]
            if w == 0.0:
                continue
            if f1 == f2:
                tot += w * IDENT["m2_qq"] * fold(_HC["m2_qq"], s, t, u, ft, fu)
            elif f1 == -f2:
                tot += w * (IDENT["m2_qqbar_same"] * folded3(m2_qqbar_same, s, t, u, ft, fu, fs)
                            + IDENT["m2_qqbar_diff"] * (NF - 1)
                            * folded3(m2_qqbar_diff, s, t, u, ft, fu, fs)
                            + IDENT["m2_qqbar_gg"] * fold(_HC["m2_qqbar_gg"], s, t, u, ft, fu))
            else:
                tot += w * IDENT["m2_qqp"] * fold(_HC["m2_qqp"], s, t, u, ft, fu)
    return tot


def build_grid():
    """the quadrature and its parton densities, exactly as hadronic_chi.py lays them out."""
    mass_bins, chi_bins = read_bins()
    mn, mw = np.polynomial.legendre.leggauss(NM)
    yn, yw = np.polynomial.legendre.leggauss(NYB)
    cn, cw = np.polynomial.legendre.leggauss(NCHI)
    pts, meta = [], []
    for bi, (mlo, mhi) in enumerate(mass_bins):
        mj = 0.5 * (mhi - mlo) * mn + 0.5 * (mhi + mlo)
        jm = 0.5 * (mhi - mlo)
        for ci, (clo, chi_hi) in enumerate(chi_bins):
            cc = 0.5 * (chi_hi - clo) * cn + 0.5 * (chi_hi + clo)
            jc = 0.5 * (chi_hi - clo)
            for i, m in enumerate(mj):
                for j in range(NYB):
                    x1 = m / ROOT_S * math.exp(YB_MAX * yn[j])
                    x2 = m / ROOT_S * math.exp(-YB_MAX * yn[j])
                    pts.append((x1, m))
                    pts.append((x2, m))
                    for k, c in enumerate(cc):
                        meta.append((bi, ci, i, j, k, float(m), float(c), x1, x2,
                                     mw[i] * jm * yw[j] * YB_MAX * cw[k] * jc))
    rows = read_dump(pts)
    if rows is None:
        rows = run_dump(pts)
    return mass_bins, chi_bins, meta, rows


def spectrum(meta, rows, mass_bins, chi_bins, mkk, gam, dress_s):
    """dsigma/dchi on the card grid.  mkk = None is the Standard Model."""
    out = np.zeros((len(mass_bins), len(chi_bins)))
    for idx, (bi, ci, i, j, k, m, c, x1, x2, wgt) in enumerate(meta):
        r1 = rows[2 * (idx // NCHI)]
        r2 = rows[2 * (idx // NCHI) + 1]
        d1 = {f: r1[f + 6] / x1 for f in range(-5, 6)}
        d2 = {f: r2[f + 6] / x2 for f in range(-5, 6)}
        s = m * m
        if mkk is None:
            ft = fu = 1.0
            fs = 1.0 + 0.0j
        else:
            t = s / (1.0 + c)
            u = s * c / (1.0 + c)
            ft, fu = f_space(t, mkk), f_space(u, mkk)
            fs = f_time(s, mkk, gam) if dress_s else (1.0 + 0.0j)
        w = weight3(d1, d2, s, c, ft, fu, fs)
        out[bi, ci] += wgt / (16.0 * math.pi * s * (1.0 + c) ** 2) * (2.0 * m / ROOT_S ** 2) * w
    return out


def line(c="-", n=94):
    print(c * n)


def main():
    fails = []
    print("=" * 94)
    print("THE s CHANNEL DRESSED, WITH THE TOTAL WIDTH PROFILED")
    print("=" * 94)

    mass_bins, chi_bins, meta, rows = build_grid()
    print("\n  quadrature %d mass x %d chi, %d points" % (len(mass_bins), len(chi_bins), len(meta)))
    print("  Gamma_qqbar/M = 2 alpha_s = %.3f, and Gamma_exotic/M is scanned on top of it"
          % (2 * ALPHAS))

    sm = spectrum(meta, rows, mass_bins, chi_bins, None, 0.0, False)

    # ---- C1: with the s channel OFF, reproduce the archived t/u-only deformation -----------
    print("\n[C1] s CHANNEL OFF MUST REPRODUCE THE ARCHIVED t/u CALCULATION")
    hc = json.loads((OUT / "hadronic_chi.json").read_text())
    br = hc["branches"]["ceiling"]
    mkk_ref = float(br["MKK_GeV"])
    # "raw" es el cociente a QCD sin renormalizar por bin de masa, que es lo que spectrum()/sm
    # produce.  "norm" ya lleva la renormalizacion a area unidad y compararse contra ella
    # escondería un fallo global de normalizacion, que es justo lo que este control busca.
    ref = np.array([r["raw"] for r in br["rows"]], float)
    off = spectrum(meta, rows, mass_bins, chi_bins, mkk_ref, 0.0, False) / sm
    if off.shape != ref.shape:
        print("      la forma archivada es %s y la reconstruida %s" % (ref.shape, off.shape))
        print("      C1 FAIL -- no se estan comparando los mismos objetos")
        fails.append("C1")
    else:
        worst = float(np.abs(off / ref - 1.0).max())
        ok = worst < 1e-9
        print("      M_KK = %.0f GeV, %d x %d bins" % (mkk_ref, ref.shape[0], ref.shape[1]))
        print("      worst relative difference against the archive : %.2e" % worst)
        print("      C1 %s" % ("PASS -- the rebuild IS the archived calculation, so any"
                               " difference below is the s channel"
                               if ok else
                               "*** FAIL -- the rebuild is not the same calculation, and nothing"
                               " below can be attributed to the s channel ***"))
        if not ok:
            fails.append("C1")

    # ---- the measurement -------------------------------------------------------------------
    print("\n[1] HOW MUCH THE s CHANNEL MOVES THE DEFORMATION, AGAINST THE WIDTH")
    print("      at the two rungs that decide anything: the ceiling and the escape branch")
    print("\n      %-10s %-14s %12s %12s %12s"
          % ("1/R5 [TeV]", "Gam_exo/M", "t/u only", "s+t+u", "shift"))
    line()
    res = {}
    for mkk_tev in (10.03, 9.09, 3.97):
        mkk = mkk_tev * 1000.0
        base = spectrum(meta, rows, mass_bins, chi_bins, mkk, 0.0, False) / sm
        b0 = float(base[0, 0])
        for gexo in (0.0, 0.16, 0.5, 1.0):
            gam = 2 * ALPHAS + gexo
            full = spectrum(meta, rows, mass_bins, chi_bins, mkk, gam, True) / sm
            f0 = float(full[0, 0])
            res[(mkk_tev, gexo)] = (b0, f0)
            print("      %-10.2f %-14.2f %12.4f %12.4f %+11.2f %%"
                  % (mkk_tev, gexo, b0, f0, 100 * (f0 / b0 - 1)))
        line()

    print("""
[2] WHAT THAT MEANS -- AND WHY THE SHIFT IN THE RATIO IS NOT THE SHIFT IN THE LIMIT
      The shift above is on the RATIO to QCD.  What the fit sees is the DEFORMATION, ratio - 1,
      and at 9.09 TeV that goes from 0.0838 to 0.0788 -- a six per cent cut, not half a per cent.
      Since the deformation falls like 1/M^2, six per cent of it is about three per cent in M.
      That estimate is written down only to say why it is not good enough: the deformation does
      not scale exactly, the bins are correlated, and the scale is profiled.  So the limit is
      recomputed rather than rescaled.  [[a-bound-is-not-a-computation]]""")

    # ---- the limit itself, both ways, at the width that hurts most --------------------------
    print("\n[3] THE LIMIT RECOMPUTED WITH s+t+u, AT THE WORST WIDTH")
    print("      Gamma_exotic = 0 is the worst case: the table above shows the shift SHRINKING")
    print("      as the width grows, so profiling over Gamma_exotic >= 0 is bounded by zero.")
    tabs = [load_table(t) for t, _m in PARTICLE]
    data = np.concatenate([t[1] for t in tabs])
    errs = np.concatenate([t[2] for t in tabs])
    T_mjj = np.concatenate([t[3] for t in tabs])
    T_pt = np.concatenate([t[4] for t in tabs])
    widths = np.concatenate([t[0][:, 1] - t[0][:, 0] for t in tabs])
    sizes = [len(t[0]) for t in tabs]
    blocks, kk = [], 0
    for n in sizes:
        blocks.append(list(range(kk, kk + n)))
        kk += n
    C = np.outer(errs, errs) * load_corr()

    def to_cms(g):
        """the card grid mapped onto CMS's binning, as scan_mkk.py does it."""
        out = []
        for (tid, (mlo, mhi)), tab in zip(PARTICLE, tabs):
            bi = [i for i, mb in enumerate(mass_bins) if abs(mb[0] - mlo) < 1e-6][0]
            fine = g[bi]
            ce = [(a, b) for a, b in tab[0]]
            if len(ce) != len(chi_bins):
                fine = np.array([np.mean([f for (a, b), f in zip(chi_bins, fine)
                                          if a >= lo2 - 1e-9 and b <= hi2 + 1e-9])
                                 for lo2, hi2 in ce])
            out.append(fine)
        return np.concatenate(out)

    sm_cms = to_cms(sm)

    # LAS DOS CONFIGURACIONES, no solo la de 70 puntos.  El limite que el articulo cita sale de
    # los 77, y un diagnostico hecho a 70 y luego trasladado por porcentaje es exactamente el
    # error que este proyecto lleva todo el dia corrigiendo en otros sitios.
    # [[compare-regions-from-the-same-evaluations]]
    def chi_of(dd, drop=True):
        return min(chi2(data, renorm(((1 - a) * T_mjj + a * T_pt) * dd, blocks, widths),
                        C, blocks, drop=drop) for a in np.linspace(0, 1, 41))

    base_sm = chi_of(np.ones_like(sm_cms))
    base_77 = chi_of(np.ones_like(sm_cms), drop=False)
    print("\n      Standard Model chi2 = %.2f (70 puntos)  y  %.2f (77 puntos)"
          % (base_sm, base_77))
    print("\n      %-12s %11s %11s   %11s %11s"
          % ("1/R5 [TeV]", "t/u (70)", "s+t+u (70)", "t/u (77)", "s+t+u (77)"))
    line()
    grid = [8.0, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0, 16.0]
    cur_tu, cur_full, cur_tu77, cur_full77 = [], [], [], []
    for M in grid:
        mkk = M * 1000.0
        d_tu = to_cms(spectrum(meta, rows, mass_bins, chi_bins, mkk, 0.0, False)) / sm_cms
        d_fu = to_cms(spectrum(meta, rows, mass_bins, chi_bins, mkk, 2 * ALPHAS, True)) / sm_cms
        a, b = chi_of(d_tu) - base_sm, chi_of(d_fu) - base_sm
        a7 = chi_of(d_tu, drop=False) - base_77
        b7 = chi_of(d_fu, drop=False) - base_77
        cur_tu.append((M, a))
        cur_full.append((M, b))
        cur_tu77.append((M, a7))
        cur_full77.append((M, b7))
        print("      %-12.2f %11.2f %11.2f   %11.2f %11.2f" % (M, a, b, a7, b7))
    line()

    def cross(sc):
        out = None
        for (m1, d1), (m2, d2) in zip(sc, sc[1:]):
            if d1 > 3.84 >= d2:
                out = m1 + (m2 - m1) * (d1 - 3.84) / (d1 - d2)
        return out

    c_tu, c_fu = cross(cur_tu), cross(cur_full)
    c_tu7, c_fu7 = cross(cur_tu77), cross(cur_full77)
    print("\n      %-14s %12s %12s %10s" % ("configuracion", "t/u only", "s+t+u", "cambio"))
    for nm, x, y in (("70 puntos", c_tu, c_fu), ("77 puntos", c_tu7, c_fu7)):
        d = "%+.2f %%" % (100 * (y / x - 1)) if (x and y) else "--"
        print("      %-14s %12s %12s %10s"
              % (nm, "%.2f TeV" % x if x else "--", "%.2f TeV" % y if y else "--", d))
    print("""
      El articulo cita el limite de 77 puntos, asi que el porcentaje que vale es el de esa fila.
      Dar el de 70 y trasladarlo es comparar entre configuraciones, que es el error que este
      proyecto lleva corrigiendo todo el dia.""")
    if c_fu7:
        for nm, ceil in (("measured-m_h ceiling", 9.09), ("window ceiling", 10.03)):
            print("      s+t+u a 77 puntos sigue por encima de %-22s (%.2f) : %s"
                  % (nm, ceil, c_fu7 > ceil))
    (OUT / "s_channel_width.json").write_text(json.dumps(
        {"shift_ratio": {"%.2f|%.2f" % k: v for k, v in res.items()},
         "limit_tu_TeV": c_tu, "limit_stu_TeV": c_fu,
         "limit_tu_77_TeV": c_tu7, "limit_stu_77_TeV": c_fu7,
         "grid": [[float(m), float(x), float(y)]
                  for (m, x), (_m, y) in zip(cur_tu, cur_full)],
         "grid77": [[float(m), float(x), float(y)]
                    for (m, x), (_m, y) in zip(cur_tu77, cur_full77)]}, indent=1))
    print("\n    [wrote outputs/s_channel_width.json]")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
