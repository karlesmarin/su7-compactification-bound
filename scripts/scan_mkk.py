#!/usr/bin/env python3
"""scan_mkk -- an SM-referenced Delta chi^2 = 3.84 sensitivity on 1/R_5, spacelike channels only.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

fit_chi2.py scored two archived values of M_KK.  A limit needs a scan, and a naive scan is
hopeless: 46080 quadrature points times forty masses times a hundred-term flavour sum is a
hundred million Python operations.

IT DOES NOT HAVE TO BE.  M_KK enters ONLY through the two form factors, and the weight is
quadratic in them:

    W(F_t, F_u) = c0 + c_t F_t^2 + c_u F_u^2 + c_tu F_t F_u + c_t1 F_t + c_u1 F_u .

So six coefficients per quadrature point carry the whole mass dependence.  They are obtained by
evaluating the exact weight at six (F_t, F_u) pairs and solving, rather than by re-deriving the
matrix elements by hand -- and C1 then evaluates a SEVENTH pair, which the fit never saw.  If
the quadratic ansatz were incomplete that check would fail; it is the reason to trust the speed.

WHAT IS SCANNED.  The spacelike deformation only: no width, no BR(G_n -> jj).  The limit is
therefore the half of the question that survives knowing nothing about the exotic masses.

Run:  python scan_mkk.py > ../outputs/scan_mkk.txt
"""
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"

_FC = {"__file__": str(HERE / "fit_chi2.py"), "__name__": "fit_chi2"}
exec(compile((HERE / "fit_chi2.py").read_text(encoding="utf-8"), "fit_chi2.py", "exec"), _FC)
load_table, load_corr, rebin = _FC["load_table"], _FC["load_corr"], _FC["rebin"]
renorm, chi2, PARTICLE = _FC["renorm"], _FC["chi2"], _FC["PARTICLE"]

_HC = {"__file__": str(HERE / "hadronic_chi.py"), "__name__": "hadronic_chi"}
exec(compile((HERE / "hadronic_chi.py").read_text(encoding="utf-8"),
             "hadronic_chi.py", "exec"), _HC)
read_bins, weight_ff = _HC["read_bins"], _HC["weight_ff"]
ROOT_S, YB_MAX = _HC["ROOT_S"], _HC["YB_MAX"]
NM, NYB, NCHI = _HC["NM"], _HC["NYB"], _HC["NCHI"]

# six (F_t, F_u) pairs, well spread so the 6x6 solve is not ill-conditioned
PROBES = [(1.0, 1.0), (1.3, 1.0), (1.0, 1.3), (1.3, 1.3), (1.7, 1.1), (1.1, 1.7)]
CHECK = (1.45, 1.22)          # the seventh, never used in the solve


def alphas_1loop(mu, asmz=0.118, mz=91.1876, nf=5):
    b0 = (33.0 - 2.0 * nf) / (12.0 * math.pi)
    return asmz / (1.0 + asmz * b0 * 2.0 * math.log(mu / mz))


def basis(ft, fu):
    return np.array([1.0, ft * ft, fu * fu, ft * fu, ft, fu])


def line(c="-", n=94):
    print(c * n)


def main():
    fails = []
    print("=" * 94)
    print("AN SM-REFERENCED SENSITIVITY ON 1/R_5 FROM THE SPACELIKE CHANNELS ALONE")
    print("=" * 94)
    # el conjunto de partones va IMPRESO: un archivo que no dice con que PDF se hizo no se puede
    # comparar con otro.  [[a-number-i-hand-over-carries-my-instrument]]
    print("\n  parton densities : %s" % _HC["PDFSET"])

    mass_bins, chi_bins = read_bins()
    mn, mw = np.polynomial.legendre.leggauss(NM)
    yn, yw = np.polynomial.legendre.leggauss(NYB)
    cn, cw = np.polynomial.legendre.leggauss(NCHI)

    # ---- the parton densities, once, at mu = M_jj -----------------------------------------
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
    print("\n  quadrature points : %d   PDF requests : %d" % (len(meta), len(pts)))
    rows = _HC["read_dump"](pts)
    # a stderr por la misma razon que en hadronic_chi.py: de donde vinieron los numeros es estado
    # del disco y no una entrada, y una salida archivada que lo lleva no reproduce en una maquina
    # limpia aunque los numeros sean identicos.  [[stale-outputs-lie]]
    if rows is None:
        sys.stderr.write("  no cache for these points; asking LHAPDF\n")
        rows = _HC["run_dump"](pts)
    else:
        sys.stderr.write("  reusing the cache keyed to these points\n")

    # ---- the six coefficients per quadrature point ----------------------------------------
    print("\n[C1] THE QUADRATIC ANSATZ, CHECKED ON A PAIR IT NEVER SAW")
    B = np.array([basis(a, b) for a, b in PROBES])
    Binv = np.linalg.inv(B)
    print("      condition number of the 6x6 solve : %.1f" % np.linalg.cond(B))

    ncoef = np.zeros((len(mass_bins), len(chi_bins), 6))
    worst = 0.0
    kpt = 0
    for idx, (bi, ci, i, j, k, m, c, x1, x2, wgt) in enumerate(meta):
        r1 = rows[2 * (idx // NCHI)]
        r2 = rows[2 * (idx // NCHI) + 1]
        d1 = {f: r1[f + 6] / x1 for f in range(-5, 6)}
        d2 = {f: r2[f + 6] / x2 for f in range(-5, 6)}
        s = m * m
        vals = np.array([weight_ff(d1, d2, s, c, a, b) for a, b in PROBES])
        co = Binv @ vals
        pref = wgt / (16.0 * math.pi * s * (1.0 + c) ** 2) * (2.0 * m / ROOT_S ** 2)
        ncoef[bi, ci] += pref * co
        if idx % 2000 == 0:
            got = float(basis(*CHECK) @ co)
            want = weight_ff(d1, d2, s, c, *CHECK)
            worst = max(worst, abs(got - want) / abs(want))
    ok = worst < 1e-9
    print("      worst relative error at (F_t, F_u) = %s : %.1e" % (str(CHECK), worst))
    print("\n      C1 %s -- the weight really is quadratic in the two form factors, so the"
          % ("PASS" if ok else "FAIL"))
    print("         scan costs six numbers a point and not one evaluation per mass.")
    if not ok:
        fails.append("C1")
        return 1

    # ---- the measurement -------------------------------------------------------------------
    # CMS set their own contact-interaction limit on M_jj > 3.6 TeV, keeping 2.4-3.0 as a
    # control and leaving 3.0-3.6 out -- the bin where they also report a 2.0 sigma local
    # discrepancy.  --ci-selection asks what our limit becomes on the same region, which is the
    # honest robustness question: a limit that needs a bin its own closure test excludes is a
    # limit leaning on the one place the data are least quiet.
    CI_SEL = "--ci-selection" in sys.argv
    PART = [x for x in PARTICLE if not (CI_SEL and abs(x[1][0] - 3000.0) < 1e-6)]
    if CI_SEL:
        print("\n  CI SELECTION: dropping 3.0-3.6 TeV, as CMS do for their CI limit")
    tabs = [load_table(t) for t, _m in PART]
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
    R = load_corr()
    if CI_SEL:
        keep, off = [], 0
        for (tid, (mlo, mhi)) in PARTICLE:
            n = len(load_table(tid)[0])
            if abs(mlo - 3000.0) > 1e-6:
                keep += list(range(off, off + n))
            off += n
        assert len(keep) == len(errs)
        R = R[np.ix_(keep, keep)]
        print("  correlation matrix sub-blocked to %d x %d" % (len(keep), len(keep)))
    C = np.outer(errs, errs) * R

    def spectrum(mkk):
        """dsigma/dchi on the card grid, then mapped onto CMS's binning."""
        out = np.zeros((len(mass_bins), len(chi_bins)))
        for bi, (mlo, mhi) in enumerate(mass_bins):
            mmid = 0.5 * (mlo + mhi)
            for ci, (clo, chi_hi) in enumerate(chi_bins):
                cmid = 0.5 * (clo + chi_hi)
                s = mmid * mmid
                if mkk is None:
                    ft = fu = 1.0
                else:
                    t = -s / (1.0 + cmid)
                    u = -s * cmid / (1.0 + cmid)
                    ft = (math.pi * math.sqrt(-t) / mkk) / math.tanh(
                        math.pi * math.sqrt(-t) / mkk)
                    fu = (math.pi * math.sqrt(-u) / mkk) / math.tanh(
                        math.pi * math.sqrt(-u) / mkk)
                out[bi, ci] = float(basis(ft, fu) @ ncoef[bi, ci])
        return out

    def to_cms(sp):
        v = []
        for (tid, (mlo, mhi)), tab in zip(PART, tabs):
            bi = [i for i, mb in enumerate(mass_bins) if abs(mb[0] - mlo) < 1e-6][0]
            f = sp[bi]
            ce = [(a, b) for a, b in tab[0]]
            if len(ce) != len(chi_bins):
                f = rebin(chi_bins, f, ce)
            v.append(f)
        return np.concatenate(v)

    sm_cms = to_cms(spectrum(None))

    # ---- the NLO K-factors, bin by bin and Lambda by Lambda -------------------------------
    NLO = "--nlo" in sys.argv
    kint_blocks = ksm_cms = None
    if NLO:
        _NK = {"__file__": str(HERE / "nlo_kfactor.py"), "__name__": "nlo_kfactor"}
        exec(compile((HERE / "nlo_kfactor.py").read_text(encoding="utf-8"),
                     "nlo_kfactor.py", "exec"), _NK)
        blocks_all, LAMS = _NK["blocks_all"], _NK["LAM"]
        raw = {}
        for mb, cb, mu0, b, a, is_lo in blocks_all():
            raw.setdefault((mb, cb), {"lo": [], "nlo": [], "mu0": mu0})
            raw[(mb, cb)]["lo" if is_lo else "nlo"].append((b, a))
        kint_blocks = raw
        # K^SM in shape: CMS's NNLO against our LO, both normalised as densities
        ks = []
        for (tid, (mlo, mhi)), tab, tn in zip(PARTICLE, tabs, [t[3] for t in tabs]):
            bi = [i for i, mb in enumerate(mass_bins) if abs(mb[0] - mlo) < 1e-6][0]
            f = spectrum(None)[bi]
            ce = [(a, b) for a, b in tab[0]]
            if len(ce) != len(chi_bins):
                f = rebin(chi_bins, f, ce)
            w = np.array([b - a for a, b in ce])
            dens = f / w
            ours = dens / np.sum(dens * w)
            k = tn / ours
            ks.append(k / np.mean(k))
        ksm_cms = np.concatenate(ks)
        print("\n  NLO mode: K^SM shape in %.3f-%.3f, K^int recomputed at every Lambda"
              % (ksm_cms.min(), ksm_cms.max()))

    def kint_cms(mkk):
        """K^int on the CMS binning, at Lambda_8 corresponding to this M_KK."""
        als = alphas_1loop(mkk)
        lam8 = math.sqrt(3.0) / (math.pi * math.sqrt(als)) * mkk
        out = []
        for (tid, (mlo, mhi)), tab in zip(PART, tabs):
            mb = [m for m in mass_bins if abs(m[0] - mlo) < 1e-6][0]
            fine = []
            for cb in chi_bins:
                v = kint_blocks[(mb, cb)]
                blo, _ = v["lo"][1]
                bn, an = v["nlo"][1]
                r = math.log(lam8 / v["mu0"])
                lo = sum(LAMS[i] * blo[i] for i in LAMS)
                nl = sum(LAMS[i] * (bn[i] + an[i] * r) for i in LAMS)
                fine.append(nl / lo if lo != 0 else 1.0)
            fine = np.array(fine)
            ce = [(a, b) for a, b in tab[0]]
            if len(ce) != len(chi_bins):
                # a K-factor is an intensive ratio, so rebinning it means averaging, not summing
                fine = np.array([np.mean([f for (a, b), f in zip(chi_bins, fine)
                                          if a >= lo2 - 1e-9 and b <= hi2 + 1e-9])
                                 for lo2, hi2 in ce])
            out.append(fine)
        return np.concatenate(out)

    def at(dd, a, drop=True, project=False):
        return chi2(data, renorm(((1 - a) * T_mjj + a * T_pt) * dd, blocks, widths),
                    C, blocks, drop=drop, project=project, widths=widths)

    def best(dd, drop=True, project=False):
        return min(at(dd, a, drop, project) for a in np.linspace(0, 1, 41))

    def best_discrete(dd, drop=True):
        """DISCRETE profiling: the better of CMS's two published scale choices, and nothing
        between them.  The continuous nuisance interpolates linearly between two columns that are
        separately-computed predictions; the interpolants are not themselves calculations, so a
        limit that leans on one of them has bought freedom rather than measured it.  If the two
        agree, the interpolation is validated empirically -- which is the only way it can be."""
        return min(at(dd, 0.0, drop), at(dd, 1.0, drop))

    base = best(np.ones_like(sm_cms))
    print("\n[1] THE SCAN, SCALE PROFILED AT EVERY POINT")
    print("      Standard Model chi2 = %.1f for %d points" % (base, len(data)))
    print("\n      %-12s %12s %12s" % ("1/R_5 [TeV]", "chi2", "Delta chi2"))
    line()
    grid = list(np.arange(3.0, 40.1, 0.5))
    scan, scan_full, absol, disc, proj, per_scale = [], [], [], [], [], []
    base_full = best(np.ones_like(sm_cms), drop=False)
    base_disc = best_discrete(np.ones_like(sm_cms))
    base_proj = best(np.ones_like(sm_cms), project=True)
    for M in grid:
        dd = to_cms(spectrum(M * 1000.0)) / sm_cms
        if NLO:
            # delta -> delta * K^int / K^SM.  The deformation is 1 + delta, so the K acts on
            # delta and not on the whole ratio: scaling (1+delta) would rescale the Standard
            # Model too.
            dd = 1.0 + (dd - 1.0) * kint_cms(M * 1000.0) / ksm_cms
        c = best(dd)
        absol.append((M, c))
        # las DOS escalas por separado, no solo la perfilada: la figura tiene que poder ensenar
        # que el minimo finito no es un artefacto de interpolar entre columnas.
        per_scale.append((M, at(dd, 0.0), at(dd, 1.0)))
        dc = c - base
        scan.append((M, dc))
        scan_full.append((M, best(dd, drop=False) - base_full))
        disc.append((M, best_discrete(dd) - base_disc))
        proj.append((M, best(dd, project=True) - base_proj))
        if abs(M - round(M)) < 1e-9 and (M <= 12 or M % 4 == 0):
            print("      %-12.1f %12.1f %12.1f" % (M, base + dc, dc))
    line()

    # ---- the crossing ----------------------------------------------------------------------
    def crossing(sc, thr=3.84):
        out = None
        for (M1, d1), (M2, d2) in zip(sc, sc[1:]):
            if d1 > thr >= d2:
                out = M1 + (M2 - M1) * (d1 - thr) / (d1 - d2)
        return out

    # ---- C3: is the Standard Model actually the best fit? -----------------------------------
    # Delta chi^2 = 3.84 is Wilks' 95 % point for one parameter measured from the MINIMUM, not
    # from an arbitrary reference.  This scan measures it from the SM, which is only the same
    # thing if the SM IS the minimum.  It need not be: CMS report a 2.0 sigma local excess in the
    # 3.0-3.6 TeV bin in the same angular direction the tower pushes, so a finite M_KK could fit
    # better than no tower at all.  That is a question with an answer, and here it is.
    # [[a-control-that-cannot-fail]]
    print("\n[C3] IS THE STANDARD MODEL THE BEST FIT, OR IS THERE A FINITE PREFERRED M_KK?")
    Mmin, cmin = min(absol, key=lambda r: r[1])
    print("      chi2 at the Standard Model (M -> infinity)   : %.2f" % base)
    print("      lowest chi2 anywhere on the scanned grid     : %.2f at 1/R_5 = %.1f TeV"
          % (cmin, Mmin))
    prefers = base - cmin
    if prefers > 0.01 and Mmin < grid[-1] - 1e-9:
        print("      the data prefer a FINITE tower by Delta chi2 = %.2f  (%.1f sigma, 1 dof)"
              % (prefers, math.sqrt(max(prefers, 0.0))))
        ref, refname = cmin, "the minimum"
    else:
        print("      the Standard Model IS the best fit on this grid: no finite preference")
        ref, refname = base, "the Standard Model, which is also the minimum"
    print("      so Wilks is referenced to %s" % refname)

    # the three numbers a reader needs in order to see how much the criterion matters
    absd = [(M, c - ref) for M, c in absol]
    lim_sm = crossing(scan, 3.84)
    lim_min384 = crossing(absd, 3.84)
    lim_min271 = crossing(absd, 2.71)
    lim_disc = crossing(disc, 3.84)
    lim_proj = crossing(proj, 3.84)
    print("\n      %-52s %10s" % ("criterion", "1/R_5 [TeV]"))
    print("      %-52s %10s" % ("Delta chi2 = 3.84 against the SM (what we quote)",
                                "%.2f" % lim_sm if lim_sm else "--"))
    print("      %-52s %10s" % ("Delta chi2 = 3.84 against the minimum (Wilks, 2-sided)",
                                "%.2f" % lim_min384 if lim_min384 else "--"))
    print("      %-52s %10s" % ("Delta chi2 = 2.71 against the minimum (1-sided, theta >= 0)",
                                "%.2f" % lim_min271 if lim_min271 else "--"))
    print("      %-52s %10s" % ("scale profiled DISCRETELY over CMS's two choices",
                                "%.2f" % lim_disc if lim_disc else "--"))
    print("      %-52s %10s" % ("covariance PROJECTED through the normalisation",
                                "%.2f" % lim_proj if lim_proj else "--"))
    print("""
      The last row is the answer to the sharpest objection there is to this fit.  HEPData calls
      the matrix "the correlation matrix of the maximum likelihood estimators of the SIGNAL
      STRENGTH MODIFIERS ... after the fit to the data", not the covariance of the normalised
      unfolded points.  If the unit-area normalisation happens AFTER that fit, the right object is
      J C J^T and not C.  That row computes it that way.  The other rows assume the parametrisations
      already agree; this one assumes they do not.""")
    cands = [x for x in (lim_sm, lim_min384, lim_min271, lim_disc, lim_proj) if x]
    if cands:
        # el rotulo se escribe con len(cands) y no con un numero a mano: decia "four" cuando ya
        # eran cinco filas, que es como un rotulo se convierte en mentira sin que nadie lo toque.
        print("\n      spread across the %d criteria : %.2f to %.2f TeV (%.1f %%)"
              % (len(cands), min(cands), max(cands), 100 * (max(cands) / min(cands) - 1)))
        print("      all %d above the measured-m_h ceiling of 9.09 TeV : %s"
              % (len(cands), all(x > 9.09 for x in cands)))
        print("      all %d above the window ceiling of 10.03 TeV      : %s"
              % (len(cands), all(x > 10.03 for x in cands)))
    crit = list(cands)      # el bloque C2, mas abajo, compara el numero citado contra estos

    cross_drop, cross_full = crossing(scan), crossing(scan_full)

    # ---- C2: does the choice of which chi^2 to use move the answer? -------------------------
    # C1b of fit_chi2.py showed the published covariance is NOT singular and does NOT project out
    # the unit-area constraint, so dropping the last bin of each mass block is a choice.  A choice
    # that changes the quoted limit has to be reported; one that does not, has to be shown not to.
    print("\n[C2] THE CHI^2 BOTH WAYS -- 70 points with one bin per block dropped, and all 77")
    print("      70 points (one dropped per mass block) : SM chi2 = %.1f, limit = %s"
          % (base, "%.2f TeV" % cross_drop if cross_drop else "no crossing"))
    print("      77 points (the full published matrix)  : SM chi2 = %.1f, limit = %s"
          % (base_full, "%.2f TeV" % cross_full if cross_full else "no crossing"))
    if cross_drop and cross_full:
        spread = abs(cross_full - cross_drop) / min(cross_full, cross_drop)
        print("      they differ by %.1f %%, and what is quoted below is the WEAKER of the two"
              % (100 * spread))
    # EL MINIMO SOBRE TODAS LAS VARIACIONES, no sobre dos.  Cuando esto solo comparaba
    # descartar-o-no, en la seleccion de CMS citaba 13.95 mientras la covarianza proyectada daba
    # 13.75 -- mas debil, y por tanto la que hay que citar.  Un "citamos el mas debil" que solo
    # mira dos de las cinco variaciones no es la regla que dice ser.
    # [[a-control-that-cannot-fail]]
    allcross = [x for x in (cross_drop, cross_full, lim_sm, lim_disc, lim_proj) if x]
    cross = min(allcross) if allcross else None
    if cross and crit:
        print("      and the quoted %.2f TeV is the weakest of ALL %d variations : %s"
              % (cross, len(allcross), all(x >= cross for x in allcross)))
        print("      (the Wilks-referenced readings of C3 are stronger still, so they do not bind)")

    print("\n[2] THE LIMIT")
    if cross:
        # NO "the 95 % point": eso es Wilks medido desde el MINIMO, y C3 mide que el minimo no
        # esta en el Modelo Estandar.  Esta cifra es un criterio referenciado al SM, y asi se
        # nombra.  Las otras cuatro lecturas estan arriba y todas caen por encima de esta.
        print("""
      chi2 - chi2_SM falls through 3.84 -- an SM-REFERENCED criterion, not Wilks' 95 %% point,
      which C3 above measures and varies -- at

          1/R_5  =  %.2f TeV .

      Below that the tower deforms chi_dijet more than the data allow.""" % cross)
    else:
        print("      no crossing inside the scanned range")
        fails.append("crossing")

    ceil_meas, ceil_all, escape = 9.09, 10.03, 3.97
    print("\n      %-46s %10s" % ("", "TeV"))
    print("      %-46s %10.2f" % ("this limit, spacelike channels only, SM-referenced", cross or 0))
    print("      %-46s %10.2f" % ("the paper's ceiling at the measured m_h", ceil_meas))
    print("      %-46s %10.2f" % ("the paper's ceiling over all contents", ceil_all))
    print("      %-46s %10.2f" % ("the branch that can pay for the escape", escape))

    # ---- the band the closure test earns, carried rather than asserted -------------------
    print("""
[2b] AND THE BAND THIS INSTRUMENT IS ENTITLED TO
      cms_ci_closure.py returned CMS's own contact-interaction limits as 20.5 against 17 TeV
      and 25.5 against 37 TeV -- ratios 1.21 and 0.69.  That is what this chi^2 is worth
      against a full likelihood, and it must be carried onto our own number rather than left
      in the other file.  [[a-number-i-hand-over-carries-my-instrument]]""")
    lo, hi = cross * 0.69, cross * 1.21
    print("\n      central                       %6.2f TeV" % cross)
    print("      scaled by the closure ratios  %6.2f  to %6.2f TeV" % (lo, hi))

    verdicts = []
    if cross and lo > escape:
        verdicts.append("the proton-decay escape branch (%.2f TeV) is excluded across the whole "
                        "band" % escape)
    elif cross and cross > escape:
        verdicts.append("the escape branch is excluded centrally but not across the band")
    for name, val in (("the ceiling at the measured m_h", ceil_meas),
                      ("the ceiling over all contents", ceil_all)):
        if cross and lo > val:
            verdicts.append("the limit clears %s (%.2f TeV) across the whole band" % (name, val))
        elif cross and cross > val:
            verdicts.append("the limit clears %s (%.2f TeV) centrally, but the low end of the "
                            "band does not" % (name, val))
    print("\n      " + "\n      ".join(verdicts) if verdicts else "\n      nothing is excluded")

    print("""
[3] WHAT THIS IS AND IS NOT
      It IS an SM-referenced Delta chi^2 = 3.84 sensitivity from the published covariance
      with the QCD scale profiled -- NOT a 95 %% CL, because C3 measures that the Standard
      Model is not the best fit and because theta = 1/M^2 >= 0 puts the null on a boundary --
      on a
      machinery that reproduces CMS's own contact-interaction limit to 21 %% and 31 %%
      (cms_ci_closure.py).  It is NOT the full statement: only the spacelike channels are
      dressed, so it is the part that survives knowing nothing about BR(G_n -> jj).

      AND THE SIGN OF THAT OMISSION IS NOT UNIFORM, which an earlier version of this paragraph
      got wrong by asserting the s-channel "would only strengthen it".  channel_breakdown.py
      measures the pure s-channel piece, q qbar -> q' qbar', at -2.4 % of the interference in
      the contact limit: dressing it there would WEAKEN this limit by about a per cent in M.
      Near the tower it is the other way and strongly so, the cotangent branch of eq. (resum)
      having poles the spacelike form factor cannot see.  So the omission is conservative where
      the escape branch sits and mildly anti-conservative at the ceiling, and neither is worth
      more than a per cent or two against the closure band.

      And it is a chi^2, not CMS's likelihood: the closure test says how much that is worth,
      and no more.""")

    # UNA CLAVE POR VARIANTE, y no un fichero que la ultima corrida pisa.  Las tres variantes
    # escribian todas en scan_mkk.json, de modo que nlo_kfactor.py leia la que hubiese corrido la
    # ultima -- la de M_jj > 3.6 TeV -- y la imprimia bajo la etiqueta "LO limit".  El propio
    # docstring de both() dice que un archivo que depende de con que bandera se corrio no se puede
    # comprobar; el JSON lo incumplia.  Lo destapo check_reproduces.py el dia que aprendio a mirar
    # en tu_limit/.  [[an-overloaded-symbol-becomes-a-false-claim]]
    OUT.mkdir(exist_ok=True)
    key = ("nlo" if NLO else "lo") + ("_ci" if CI_SEL else "_all")
    if _HC["PDFSET"] != "CT10nlo":
        key += "_" + _HC["PDFSET"]
    path = OUT / "scan_mkk.json"
    allv = json.loads(path.read_text()) if path.exists() else {}
    if "chi2_SM" in allv:                     # el formato viejo, de una sola variante
        allv = {}
    allv[key] = {"chi2_SM": base, "n_points": int(len(data)), "limit_TeV": cross,
                 "limit_drop_TeV": cross_drop, "limit_full_TeV": cross_full,
                 "chi2_min": cmin, "M_at_min_TeV": Mmin,
                 "limit_wilks_384_TeV": lim_min384, "limit_wilks_271_TeV": lim_min271,
                 "limit_discrete_scale_TeV": lim_disc, "limit_projected_cov_TeV": lim_proj,
                 "scan": [[float(a), float(b)] for a, b in scan],
                 # LA CURVA DE 77 PUNTOS TAMBIEN, y no solo su cruce.  Sin ella no se puede leer
                 # el peine sobre el chi^2 DEBIL, que es el que produce el limite citado; y una
                 # sonda que la busca y no la encuentra se queda con su centinela y concluye lo
                 # que a uno le conviene.  Paso exactamente eso.
                 # [[a-helpers-sentinel-has-a-meaning]]
                 "scan_full": [[float(a), float(b)] for a, b in scan_full],
                 "chi2_abs": [[float(a), float(b)] for a, b in absol],
                 "chi2_per_scale": [[float(a), float(b), float(c)] for a, b, c in per_scale],
                 "chi2_SM_per_scale": [float(at(np.ones_like(sm_cms), 0.0)),
                                       float(at(np.ones_like(sm_cms), 1.0))]}
    path.write_text(json.dumps(allv, indent=1, sort_keys=True))
    print("\n    [wrote outputs/scan_mkk.json, key '%s']" % key)

    print("\n" + "=" * 94)
    if fails:
        print("VERDICT: %d FAILED: %s" % (len(fails), ", ".join(fails)))
        print("=" * 94)
        return 1
    print("VERDICT: 1/R_5 > %.2f TeV, SM-referenced Delta chi^2 = 3.84, spacelike only." % cross)
    print("=" * 94)
    return 0


def both():
    """Every variant the paper quotes, in ONE invocation.

    A script whose archive depends on which flag it was last run with is a script whose archive
    cannot be checked.  [[a-patch-on-disk-is-not-a-patch-that-ran]]"""
    rc = 0
    for nlo, ci, label in (
            (False, False, "LEADING ORDER, all seven mass bins"),
            (True, False, "NLO, all seven mass bins"),
            (True, True, "NLO, on CMS's own contact-interaction range M_jj > 3.6 TeV")):
        for flag, want in (("--nlo", nlo), ("--ci-selection", ci)):
            if want and flag not in sys.argv:
                sys.argv.append(flag)
            if not want and flag in sys.argv:
                sys.argv.remove(flag)
        print("\n" + "#" * 94)
        print("# %s" % label)
        print("#" * 94)
        rc |= main()
    return rc


if __name__ == "__main__":
    sys.exit(both())
