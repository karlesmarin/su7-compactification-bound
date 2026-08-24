#!/usr/bin/env python3
"""cms_ci_closure -- can this machinery reproduce CMS's OWN contact-interaction limit?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY.  fit_chi2.py returns Delta chi2 = 51 against the tower at the paper's own ceiling, which
would exclude the whole model class.  A number that strong has to be earned, and the way to earn
it is to point the same machinery at a limit somebody else has already published and see whether
it comes back.  If it does, the 51 can be believed.  If it does not, the error is ours and this
is where it shows.

WHAT IS BEING REPRODUCED.  CMS-EXO-24-011 quotes 95 % CL limits on the quark compositeness scale
of the standard colour-SINGLET left-handed contact interaction, 17 TeV for destructive and 37 TeV
for constructive interference.  Those numbers are in their text, not in the HEPData record --
the three limit tables there are a DM mediator, an ALP gluon coupling and an anomalous triple
gluon coupling, none of which is a four-quark contact operator.

HOW, AND WHY NOT WITH OUR OWN COLOUR ALGEBRA.  A colour-singlet operator interferes with the
QCD gluon only through the crossed terms of identical quarks -- Tr T^a = 0 kills the direct
contraction -- so writing that interference by hand is exactly the kind of colour bookkeeping
that goes wrong silently.  It does not have to be written: CIJET already computed it, and its
coefficients decode with the manual's own parity identities (b_1 = b_5, b_11 = b_55, b_12 = b_56)
as a check that the slots are being read right.  So the CI side is theirs and only the Standard
Model denominator and the chi^2 are ours.

    sigma_CI,i(Lambda) = eta b_1,i (5 TeV)^2 / Lambda^2  +  b_11,i (5 TeV)^4 / Lambda^4

WHAT A PASS LOOKS LIKE.  Not three digits.  CMS's number comes from a full likelihood with their
systematics, their unfolding and their nuisances; this is a chi^2 over the published covariance
with one nuisance.  Landing within a few tens of per cent means the machinery is sound; landing
at 60 TeV, or at 5, means it is not.

Run:  python cms_ci_closure.py > ../outputs/cms_ci_closure.txt
"""
import json
import math
import pathlib
import re
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"
GRID = pathlib.Path("E:/proyectos/Curiosity/tools/cijet/CIJET/data/grid_fitALresults.dat")

# CMS's CI limit is set on M_jj > 3.6 TeV, with 2.4-3.0 TeV as a control region: the 3.0-3.6
# bin does NOT enter the 17/37 TeV numbers.  A closure that keeps it is not closing on their
# analysis.  --ci-selection drops it.
import sys as _s
CI_SEL = "--ci-selection" in _s.argv

_FC = {"__file__": str(HERE / "fit_chi2.py"), "__name__": "fit_chi2"}
exec(compile((HERE / "fit_chi2.py").read_text(encoding="utf-8"), "fit_chi2.py", "exec"), _FC)
load_table, load_corr, rebin = _FC["load_table"], _FC["load_corr"], _FC["rebin"]
renorm, chi2, PARTICLE = _FC["renorm"], _FC["chi2"], _FC["PARTICLE"]

_CC = {"__file__": str(HERE / "cijet_control.py"), "__name__": "cijet_control"}
exec(compile((HERE / "cijet_control.py").read_text(encoding="utf-8"),
             "cijet_control.py", "exec"), _CC)
read_bins, cijet_scale = _CC["read_bins"], _CC["cijet_scale"]
_HC = _CC["_HC"]
ROOT_S, YB_MAX, GEV2_TO_PB = _CC["ROOT_S"], _CC["YB_MAX"], _CC["GEV2_TO_PB"]
NM, NYB, NCHI = _CC["NM"], _CC["NYB"], _CC["NCHI"]
FIVE2 = 5000.0 ** 2
FIVE4 = 5000.0 ** 4

NUM = r"[-+]?[0-9]*\.?[0-9]+(?:[EeDd][-+]?[0-9]+)?"


def cijet_lo_central():
    """b_1 (linear, slot 0) and b_11 (quadratic, slot 12) per (mass bin, chi bin), LO, mu_f = 1."""
    lines = GRID.read_text().splitlines()
    blocks, mass_now, i = [], None, 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("mass=["):
            mass_now = tuple(float(v) for v in re.findall(NUM, s))
            i += 1
            continue
        p = s.split()
        if mass_now and len(p) == 3 and all(_num(v) for v in p):
            nums = []
            for k in (1, 2, 3):
                nums += [float(v.replace("D", "E")) for v in lines[i + k].split()]
            if len(nums) == 38:
                blocks.append((mass_now, (float(p[0]), float(p[1])), nums))
                i += 4
                continue
        i += 1
    seen = {}
    for mb, cb, nums in blocks:
        if max(abs(x) for x in nums[1:12:2]) == 0.0:          # LO: every a vanishes
            seen.setdefault((mb, cb), []).append(nums)
    out = {}
    for k, lst in seen.items():
        if len(lst) == 3:
            n = lst[1]                                        # mu_f = 1.0, the middle of three
            # b3 y b33 hacen falta para el benchmark VECTORIAL, que es el que se parece a la
            # torre: (q gamma_mu T^a q)^2 ve las dos quiralidades, y VV es lambda_1 = lambda_3 =
            # lambda_5 donde el cierre LL es lambda_1 solo.  Sus ranuras NO se pueden verificar
            # con las identidades de paridad del manual -- b(a)1,2 = b(a)5,6 y
            # b(a)11,22,12 = b(a)55,66,56 no tocan el 3 ni el 4, porque O_3 y O_4 son LR y la
            # paridad los manda a si mismos.  Se INFIEREN del patron anclado (b11=12, b22=14,
            # b12=16 | b33=18, b44=20, b34=22 | b55=24, b66=26, b56=28), y quien las valida es el
            # propio cierre contra los 50 TeV publicados de Lambda^-_VV: con la ranura mal, no
            # sale ese numero.  [[a-number-i-hand-over-carries-my-instrument]]
            out[k] = {"b1": n[0], "b11": n[12], "b5": n[8],
                      "b55": n[24], "b22": n[14], "b66": n[26],
                      "b12": n[16], "b56": n[28],
                      "b3": n[4], "b33": n[18]}
    return out


def _num(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


def sm_absolute(mass_bins, chi_bins):
    """the Standard-Model dsigma/dchi per bin, in pb, at CIJET's scale, from the cached dump."""
    weight = _HC["weight"]
    mn, mw = np.polynomial.legendre.leggauss(NM)
    yn, yw = np.polynomial.legendre.leggauss(NYB)
    cn, cw = np.polynomial.legendre.leggauss(NCHI)

    # build the request ourselves and go through hadronic_chi's KEYED cache, so this cannot
    # pick up a dump made for a different scale.  That collision already happened once.
    pts = []
    for (mlo, mhi) in mass_bins:
        mj = 0.5 * (mhi - mlo) * mn + 0.5 * (mhi + mlo)
        for (clo, chi_hi) in chi_bins:
            cc = 0.5 * (chi_hi - clo) * cn + 0.5 * (chi_hi + clo)
            for m in mj:
                for j in range(NYB):
                    x1 = m / ROOT_S * math.exp(YB_MAX * yn[j])
                    x2 = m / ROOT_S * math.exp(-YB_MAX * yn[j])
                    for c in cc:
                        mu = cijet_scale(m, float(c))
                        pts.append((x1, mu))
                        pts.append((x2, mu))
    rows = _HC["read_dump"](pts)
    # a stderr, igual que en hadronic_chi.py y scan_mkk.py: si los numeros salieron del cache o
    # de LHAPDF es estado del disco, no una entrada.  [[stale-outputs-lie]]
    if rows is None:
        sys.stderr.write("      no cache for these %d points; asking LHAPDF\n" % len(pts))
        rows = _HC["run_dump"](pts)
    else:
        sys.stderr.write("      reusing the cache keyed to these %d points\n" % len(pts))
    k = 0
    out = np.zeros((len(mass_bins), len(chi_bins)))
    for bi, (mlo, mhi) in enumerate(mass_bins):
        mj = 0.5 * (mhi - mlo) * mn + 0.5 * (mhi + mlo)
        jm = 0.5 * (mhi - mlo)
        for ci, (clo, chi_hi) in enumerate(chi_bins):
            cc = 0.5 * (chi_hi - clo) * cn + 0.5 * (chi_hi + clo)
            jc = 0.5 * (chi_hi - clo)
            acc = 0.0
            for i, m in enumerate(mj):
                s = m * m
                for j in range(NYB):
                    x1 = m / ROOT_S * math.exp(YB_MAX * yn[j])
                    x2 = m / ROOT_S * math.exp(-YB_MAX * yn[j])
                    for kk, c in enumerate(cc):
                        r1, r2 = rows[k], rows[k + 1]
                        k += 2
                        d1 = {f: r1[f + 6] / x1 for f in range(-5, 6)}
                        d2 = {f: r2[f + 6] / x2 for f in range(-5, 6)}
                        als = r1[13]
                        w = weight(d1, d2, s, float(c), None)
                        acc += (mw[i] * jm * yw[j] * YB_MAX * cw[kk] * jc
                                * (2.0 * m / ROOT_S ** 2) * (4.0 * math.pi * als) ** 2
                                * w / (16.0 * math.pi * s * (1.0 + c) ** 2))
            out[bi, ci] = acc * GEV2_TO_PB
    return out


def line(c="-", n=94):
    print(c * n)


def main():
    fails = []
    print("=" * 94)
    print("CLOSURE: CMS's OWN CONTACT-INTERACTION LIMIT, THROUGH OUR MACHINERY")
    print("=" * 94)

    card_mass, card_chi = read_bins()
    coef = cijet_lo_central()
    print("\n  CIJET LO central-scale coefficients for %d bins" % len(coef))

    # ---- C1: the manual's parity identities, on the slots we are about to use -------------
    print("\n[C1] ARE WE READING THE RIGHT SLOTS?")
    print("""      The manual states b_1 = b_5, b_11 = b_55, b_22 = b_66 and b_12 = b_56, forced
      by QCD parity.  If the slot mapping were wrong none of them would hold.""")
    w1 = max(abs(v["b1"] - v["b5"]) for v in coef.values())
    w2 = max(abs(v["b11"] - v["b55"]) for v in coef.values())
    w3 = max(abs(v["b22"] - v["b66"]) for v in coef.values())
    w4 = max(abs(v["b12"] - v["b56"]) for v in coef.values())
    print("\n      worst |b_1-b_5|   %.1e" % w1)
    print("      worst |b_11-b_55| %.1e" % w2)
    print("      worst |b_22-b_66| %.1e" % w3)
    print("      worst |b_12-b_56| %.1e" % w4)
    ok = max(w1, w2, w3, w4) < 1e-12
    print("\n      C1 %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C1")
        return 1

    # ---- the Standard Model denominator ---------------------------------------------------
    print("\n[1] THE STANDARD-MODEL DENOMINATOR, ABSOLUTE, AT CIJET's SCALE")
    sm = sm_absolute(card_mass, card_chi)
    print("      dsigma/dchi in the first bin of each mass bin, pb:")
    for bi, mb in enumerate(card_mass):
        print("        %5.1f-%-6.1f TeV   %12.4e" % (mb[0] / 1000, mb[1] / 1000, sm[bi, 0]))

    # ---- map both onto the CMS particle-level binning -------------------------------------
    PART = [x for x in PARTICLE if not (CI_SEL and abs(x[1][0] - 3000.0) < 1e-6)]
    if CI_SEL:
        print("\n  CI SELECTION: dropping 3.0-3.6 TeV, which CMS excludes from its CI limit")
    tabs = [load_table(t) for t, _m in PART]
    data = np.concatenate([t[1] for t in tabs])
    errs = np.concatenate([t[2] for t in tabs])
    T_mjj = np.concatenate([t[3] for t in tabs])
    T_pt = np.concatenate([t[4] for t in tabs])
    widths = np.concatenate([t[0][:, 1] - t[0][:, 0] for t in tabs])
    sizes = [len(t[0]) for t in tabs]
    blocks, k = [], 0
    for n in sizes:
        blocks.append(list(range(k, k + n)))
        k += n
    # The correlation matrix is indexed by the FULL 77-point binning, so dropping a mass bin
    # means taking the corresponding sub-block -- not just a smaller matrix.  Getting this wrong
    # is silent: numpy would broadcast a 65-vector against a 77x77 only by accident.
    R = load_corr()
    if CI_SEL:
        keep, off = [], 0
        for (tid, (mlo, mhi)) in PARTICLE:
            n = len(load_table(tid)[0])
            if not (abs(mlo - 3000.0) < 1e-6):
                keep += list(range(off, off + n))
            off += n
        assert len(keep) == len(errs), "sub-block %d vs data %d" % (len(keep), len(errs))
        R = R[np.ix_(keep, keep)]
        print("      correlation matrix sub-blocked to %d x %d" % (len(keep), len(keep)))
    C = np.outer(errs, errs) * R

    # LL es lambda_1 solo; VV es lambda_1 = lambda_3 = lambda_5.  Con lambda_2 = lambda_4 =
    # lambda_6 = 0 NINGUN termino cruzado de su ec. (3) sobrevive -- los que existen son (1,2),
    # (3,4), (5,6) y los (i,4), y todos tocan un indice par -- de modo que
    #     sigma_VV = lambda (b1+b3+b5)/L^2 + lambda^2 (b11+b33+b55)/L^4 .
    VV = "--vv" in sys.argv
    LIN = ("b1", "b3", "b5") if VV else ("b1",)
    QUA = ("b11", "b33", "b55") if VV else ("b11",)
    if VV:
        print("\n  VV: lambda_1 = lambda_3 = lambda_5, sin terminos cruzados (su ec. 3)")

    b1, b11, smv = [], [], []
    for (tid, (mlo, mhi)), tab in zip(PART, tabs):
        bi = [i for i, mb in enumerate(card_mass) if abs(mb[0] - mlo) < 1e-6][0]
        f1 = np.array([sum(coef[(card_mass[bi], cb)][k] for k in LIN) for cb in card_chi])
        f11 = np.array([sum(coef[(card_mass[bi], cb)][k] for k in QUA) for cb in card_chi])
        fsm = sm[bi]
        ce = [(a, b) for a, b in tab[0]]
        if len(ce) != len(card_chi):
            f1 = rebin(card_chi, f1, ce)
            f11 = rebin(card_chi, f11, ce)
            fsm = rebin(card_chi, fsm, ce)
        b1.append(f1)
        b11.append(f11)
        smv.append(fsm)
    b1 = np.concatenate(b1)
    b11 = np.concatenate(b11)
    smv = np.concatenate(smv)
    print("\n      mapped onto %d CMS points (%s)" % (len(b1), sizes))

    # ---- the scan --------------------------------------------------------------------------
    # EXPECTED (Asimov) o OBSERVADO.  El cierre solo contra el limite observado no puede separar
    # un defecto del metodo de una fluctuacion de los datos: si nuestro numero sale bajo, puede
    # ser que la maquinaria pierda sensibilidad o que estos datos concretos fueran generosos.  El
    # Asimov ---sustituir los datos por la propia prediccion del Modelo Estandar--- quita la
    # fluctuacion y deja solo el metodo, y CMS publica sus esperados, 18 y 34 TeV, contra los que
    # comparar.  [[a-control-that-cannot-fail]]
    ASIMOV = "--asimov" in sys.argv
    obs = data
    if ASIMOV:
        a_hat = min(np.linspace(0, 1, 41),
                    key=lambda a: chi2(data, renorm((1 - a) * T_mjj + a * T_pt, blocks, widths),
                                       C, blocks))
        obs = renorm((1 - a_hat) * T_mjj + a_hat * T_pt, blocks, widths)
        print("\n  ASIMOV: los datos se sustituyen por la prediccion del Modelo Estandar de")
        print("  CMS con la escala en su mejor ajuste, a = %.2f.  La covarianza no se toca."
              % a_hat)

    def best_chi2(dd):
        return min(chi2(obs, renorm(((1 - a) * T_mjj + a * T_pt) * dd, blocks, widths),
                        C, blocks)
                   for a in np.linspace(0, 1, 41))

    base = best_chi2(np.ones_like(b1))
    print("\n[2] THE SCAN OVER Lambda, WITH THE SCALE PROFILED")
    print("      Standard Model chi2 = %.1f for %d points" % (base, len(data)))
    print("\n      %-8s %12s %12s   %-8s %12s %12s"
          % ("L [TeV]", "eta=+1 chi2", "Dchi2", "L [TeV]", "eta=-1 chi2", "Dchi2"))
    line()
    lim = {}
    for eta in (+1.0, -1.0):
        rows_out = []
        for L in np.arange(5.0, 60.1, 0.5):
            L2 = (L * 1000.0) ** 2
            d = 1.0 + (eta * b1 * FIVE2 / L2 + b11 * FIVE4 / L2 ** 2) / smv
            rows_out.append((L, best_chi2(d) - base))
        # the 95% point: the largest L still excluded, scanning down from high L
        cross = None
        for L, dc in rows_out:
            if dc > 3.84:
                cross = L
        lim[eta] = (cross, rows_out)
    for i, L in enumerate(np.arange(6.0, 42.1, 4.0)):
        row = []
        for eta in (+1.0, -1.0):
            dc = dict(lim[eta][1]).get(round(L * 2) / 2, float("nan"))
            row += ["%.1f" % L, "%.1f" % (base + dc), "%.1f" % dc]
        print("      %-8s %12s %12s   %-8s %12s %12s" % tuple(row))
    line()

    print("""
[3] THE COMPARISON THAT DECIDES WHETHER ANYTHING ELSE HERE CAN BE BELIEVED""")
    # los publicados, de su Tabla 5.  LL/RR: 17 y 37 observados, 18 y 34 esperados.
    # VV: Lambda^+ observado sale DISJUNTO ("<19 and 27-31 TeV"), asi que como objetivo limpio
    # solo sirve el Lambda^- = 50 TeV; los esperados son 41 y 45.
    if VV:
        THEIRS = (41.0, 45.0) if ASIMOV else (float("nan"), 50.0)
    else:
        THEIRS = (18.0, 34.0) if ASIMOV else (17.0, 37.0)
    what = ("VV, " if VV else "LL/RR, ") + ("expected (Asimov)" if ASIMOV else "observed")
    print("\n      %-34s %14s %14s" % ("(%s)" % what, "eta = +1", "eta = -1"))
    print("      %-34s %14s %14s" % ("CMS, published (their text)",
                                     "%.0f TeV" % THEIRS[0], "%.0f TeV" % THEIRS[1]))
    print("      %-34s %14s %14s" % ("this machinery, Dchi2 = 3.84",
                                     "%.1f TeV" % lim[+1.0][0] if lim[+1.0][0] else "none",
                                     "%.1f TeV" % lim[-1.0][0] if lim[-1.0][0] else "none"))
    r1 = (lim[+1.0][0] / THEIRS[0]) if lim[+1.0][0] else float("nan")
    r2 = (lim[-1.0][0] / THEIRS[1]) if lim[-1.0][0] else float("nan")
    print("      %-34s %14.2f %14.2f" % ("ratio ours / theirs", r1, r2))
    ok2 = all(0.6 < r < 1.7 for r in (r1, r2) if r == r)
    print("""
      A pass is not three digits: CMS's numbers come from a full likelihood with their
      systematics, their unfolding and their nuisances, and this is a chi^2 over the published
      covariance with one.  Landing within a few tens of per cent says the machinery is sound.""")
    print("\n   C2  the published contact-interaction limit is reproduced  %s"
          % ("PASS" if ok2 else "FAIL"))
    if not ok2:
        fails.append("C2")

    # una clave por modo, no un fichero que la ultima corrida pisa -- la misma leccion que
    # scan_mkk.json costo hace unas horas.  [[an-overloaded-symbol-becomes-a-false-claim]]
    OUT.mkdir(exist_ok=True)
    path = OUT / "cms_ci_closure.json"
    allm = json.loads(path.read_text()) if path.exists() else {}
    if "cms_published" in allm:
        allm = {}
    allm["asimov" if ASIMOV else "observed"] = {
        "cms_published": {"eta=+1": THEIRS[0], "eta=-1": THEIRS[1]},
        "ours": {"eta=+1": lim[+1.0][0], "eta=-1": lim[-1.0][0]},
        "ratio": {"eta=+1": r1, "eta=-1": r2},
        "chi2_SM": base, "n_points": int(len(data))}
    path.write_text(json.dumps(allm, indent=1, sort_keys=True))
    print("\n    [wrote outputs/cms_ci_closure.json, key '%s']"
          % ("asimov" if ASIMOV else "observed"))

    print("\n" + "=" * 94)
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        print("         Until this closes, fit_chi2.py's Delta chi2 = 51 is not to be quoted.")
        print("=" * 94)
        return 1
    print("VERDICT: the machinery returns CMS's own published limit, so what it says about the")
    print("         tower is worth reading.")
    print("=" * 94)
    return 0


def both():
    """Both selections in one invocation: all seven mass bins, and CMS's own CI range.

    The second is the one that closes -- 1.03 and 0.80 against 1.21 and 0.69 -- and keeping the
    first in the archive is what makes that visible rather than asserted."""
    global CI_SEL
    rc = 0
    for ci, label in ((False, "ALL SEVEN MASS BINS (not CMS's CI selection)"),
                      (True, "CMS's OWN CI SELECTION, M_jj > 3.6 TeV")):
        CI_SEL = ci
        print("\n" + "#" * 94)
        print("# %s" % label)
        print("#" * 94)
        rc |= main()
    return rc


if __name__ == "__main__":
    sys.exit(both())
