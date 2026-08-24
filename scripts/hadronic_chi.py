#!/usr/bin/env python3
"""hadronic_chi -- the tower's deformation of chi_dijet at hadron level, spacelike channels only.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT THIS IS FOR.  kk_resummation.py gives the tower's effect on ONE subprocess at parton level.
A limit needs the inclusive sample: parton densities, every subprocess including the ones the
tower does not touch, and the per-mass-bin normalisation CMS applies.  All three dilute the
parton-level factor, and the point of this file is to measure by how much rather than to guess.

WHY THE SPACELIKE CHANNELS ONLY.  In t and u every mode is below its two-parton threshold, so
the propagator is real and the tower enters as the exact form factor F = pi a coth(pi a) of
kk_resummation.py -- no truncation, no width, and NO dependence on BR(G_n -> jj), which is the
one number [KM25] does not supply.  The s-channel is left at its Standard-Model value here.  So
what comes out is a deformation that survives knowing nothing at all about the exotic masses.
That is the whole point; it is not an approximation to the full answer, it is a different and
weaker question, asked because it is the one that can be answered.

WHICH SUBPROCESSES THE TOWER TOUCHES, AND WHY ONLY THOSE.  A Kaluza-Klein gluon couples to the
brane quarks with sqrt(2) g_s and cannot couple to two zero-mode gluons -- five-momentum in x5
forbids g0 g0 G(n).  So it dresses a gluon propagator only where that propagator joins two QUARK
lines.  That leaves the four-quark processes; qg -> qg, gg -> gg, q qbar -> gg and gg -> q qbar
are untouched, their internal gluon ending on a triple-gluon vertex or being a quark line.

  q q' -> q q'   (distinct)      t only            -> F_t^2
  q qbar' -> q qbar' (distinct)  t only            -> F_t^2
  q q -> q q     (identical)     t, u, interf.     -> F_t^2, F_u^2, F_t F_u
  q qbar -> q qbar (same flav.)  t, s, interf.     -> F_t^2, 1, F_t
  q qbar -> q' qbar' (distinct)  s only            -> 1  (left SM)

PDFs.  CT10nlo through LHAPDF inside the cijet:1.0 image -- the same set the CIJET NLO grid was
run with, so the two are not on different parton densities.  The quadrature points are known in
advance, so LHAPDF is asked for exactly those and nothing is interpolated.

BINNING.  The eight mass bins and twelve chi bins of CMS-EXO-24-011, read from the CIJET card
that was checked against their selection, not retyped here.  At leading order in this phase
space the pT and rapidity cuts are inactive -- C1 proves it -- so the fiducial region is exactly
|y_boost| < 1.11 and 1 < chi < 16.

Run:  python hadronic_chi.py            (needs docker for the PDF dump; --reuse skips it)
"""
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent

# CT10nlo por defecto porque es el conjunto con el que se genero la rejilla NLO de CIJET, y el
# control contra CIJET tiene que correr con el suyo o no es un control.  Para el LIMITE, en
# cambio, la base es la de CMS -- NNPDF3.1 -- y el efecto de esa diferencia se mide en vez de
# suponerse: SU7_PDF=NNPDF31_nnlo_as_0118 o CT14nlo.  Los ficheros viven fuera del repositorio,
# en tools/lhapdf/, y se montan en el contenedor.
PDFSET = os.environ.get("SU7_PDF", "CT10nlo")
PDFDIR = pathlib.Path("E:/proyectos/Curiosity/tools/lhapdf")
PART7 = HERE.parent
OUT = PART7 / "outputs"
CARD = pathlib.Path("E:/proyectos/Curiosity/tools/cijet/CIJET/data/bininput.card")

ROOT_S = 13000.0          # GeV
YB_MAX = 1.11
NF = 5
MKK = {"ceiling": 9090.0, "escape": 3970.0}      # GeV

# quadrature: Gauss-Legendre in each direction
NM, NYB, NCHI = 8, 10, 6


def read_bins():
    """the CMS binning, from the card that was verified against their paper."""
    mass, chi, mode = [], [], None
    for ln in CARD.read_text().splitlines():
        s = ln.split()
        if not s:
            continue
        if s[0] == "massbin":
            mode = "m"
            continue
        if s[0] == "rapbin":
            mode = "c"
            continue
        if mode == "m":
            mass.append((float(s[0]), float(s[1])))
        elif mode == "c":
            chi.append((float(s[0]), float(s[1])))
    return mass, chi


# ---------------------------------------------------------------- matrix elements
# Ellis-Stirling-Webber Table 7.1, in units of g_s^4, summed over final and averaged over
# initial colours and spins.  Each is split by channel so the form factor can dress the
# propagator that actually carries it.

def m2_qqp(s, t, u, ft, fu):
    """q q' -> q q', distinct flavours (also q qbar' -> q qbar'): t-channel only."""
    return (4.0 / 9.0) * (s * s + u * u) / (t * t) * ft * ft


def m2_qq(s, t, u, ft, fu):
    """q q -> q q, identical: t, u and their interference."""
    return ((4.0 / 9.0) * ((s * s + u * u) / (t * t) * ft * ft
                           + (s * s + t * t) / (u * u) * fu * fu)
            - (8.0 / 27.0) * s * s / (u * t) * ft * fu)


def m2_qqbar_same(s, t, u, ft, fu):
    """q qbar -> q qbar: t (dressed), s (left SM) and their interference (one t leg)."""
    return ((4.0 / 9.0) * ((s * s + u * u) / (t * t) * ft * ft
                           + (t * t + u * u) / (s * s))
            - (8.0 / 27.0) * u * u / (s * t) * ft)


def m2_qqbar_diff(s, t, u, ft, fu):
    """q qbar -> q' qbar', distinct: s-channel only, so untouched here."""
    return (4.0 / 9.0) * (t * t + u * u) / (s * s)


def m2_qqbar_gg(s, t, u, ft, fu):
    return (32.0 / 27.0) * (t * t + u * u) / (t * u) - (8.0 / 3.0) * (t * t + u * u) / (s * s)


def m2_gg_qqbar(s, t, u, ft, fu):
    return (1.0 / 6.0) * (t * t + u * u) / (t * u) - (3.0 / 8.0) * (t * t + u * u) / (s * s)


def m2_qg(s, t, u, ft, fu):
    return -(4.0 / 9.0) * (s * s + u * u) / (s * u) + (u * u + s * s) / (t * t)


def m2_gg(s, t, u, ft, fu):
    return 4.5 * (3.0 - t * u / (s * s) - s * u / (t * t) - s * t / (u * u))


def form_factors(s, chi, mkk):
    """F on the two spacelike legs.  mkk = None gives the Standard Model, F = 1."""
    if mkk is None:
        return 1.0, 1.0
    t = -s / (1.0 + chi)
    u = -s * chi / (1.0 + chi)
    ft = (math.pi * math.sqrt(-t) / mkk) / math.tanh(math.pi * math.sqrt(-t) / mkk)
    fu = (math.pi * math.sqrt(-u) / mkk) / math.tanh(math.pi * math.sqrt(-u) / mkk)
    return ft, fu


def folded(fn, s, t, u, ft, fu):
    """SUM the two jet assignments.

    chi = exp|y1-y2| >= 1 covers only cos(theta*) >= 0, that is half the angular range.  The
    other half maps onto the same chi with t and u exchanged, and it is a DIFFERENT physical
    configuration -- so the two add.  The 1/2 that belongs here for some processes is the
    identical-particle symmetry factor, applied by IDENT below, and not a property of the fold.

    The first version of this file averaged unconditionally.  That is right for q q -> q q and
    wrong for every distinct final state, which is 42 % of the interference; measured against
    CIJET it came out a factor 1.435 low.  [[a-constant-difference-is-an-object-mismatch]]"""
    return fn(s, t, u, ft, fu) + fn(s, u, t, fu, ft)


# 1/2 exactly when the two final-state particles are identical.  A quark and an antiquark are
# not, so q qbar -> q qbar does NOT get it; two gluons are.
IDENT = {"m2_qqp": 1.0, "m2_qq": 0.5, "m2_qqbar_same": 1.0, "m2_qqbar_diff": 1.0,
         "m2_qqbar_gg": 0.5, "m2_gg_qqbar": 1.0, "m2_qg": 1.0, "m2_gg": 0.5}


def weight(xf1, xf2, s, chi, mkk):
    """sum |M|^2 over every initial state, weighted by the two parton densities.

    xf1, xf2 are dicts flavour -> x f(x), flavour in -5..5 with 0 the gluon."""
    ft, fu = form_factors(s, chi, mkk)
    return weight_ff(xf1, xf2, s, chi, ft, fu)


def weight_ff(xf1, xf2, s, chi, ft, fu):
    """the same sum, with the two form factors given directly.

    Split out so that a scan over M_KK can reuse one evaluation of the parton densities: M_KK
    enters ONLY through (ft, fu), and the sum is quadratic in them, so six coefficients per
    quadrature point carry the whole M_KK dependence."""
    t = -s / (1.0 + chi)
    u = -s * chi / (1.0 + chi)
    q = [f for f in range(-NF, NF + 1) if f != 0]
    tot = 0.0

    # gluon-gluon
    g1, g2 = xf1[0], xf2[0]
    tot += g1 * g2 * (IDENT["m2_gg"] * folded(m2_gg, s, t, u, ft, fu)
                      + IDENT["m2_gg_qqbar"] * NF * folded(m2_gg_qqbar, s, t, u, ft, fu))
    # quark-gluon, both orderings
    for f in q:
        tot += IDENT["m2_qg"] * (xf1[f] * g2 + g1 * xf2[f]) * folded(m2_qg, s, t, u, ft, fu)
    # quark-quark and quark-antiquark
    for f1 in q:
        for f2 in q:
            w = xf1[f1] * xf2[f2]
            if w == 0.0:
                continue
            if f1 == f2:                                   # q q  (or qbar qbar) -- identical
                tot += w * IDENT["m2_qq"] * folded(m2_qq, s, t, u, ft, fu)
            elif f1 == -f2:                                # q qbar, same flavour
                tot += w * (IDENT["m2_qqbar_same"] * folded(m2_qqbar_same, s, t, u, ft, fu)
                            + IDENT["m2_qqbar_diff"] * (NF - 1)
                            * folded(m2_qqbar_diff, s, t, u, ft, fu)
                            + IDENT["m2_qqbar_gg"] * folded(m2_qqbar_gg, s, t, u, ft, fu))
            else:                                          # distinct flavours, t-channel
                tot += w * IDENT["m2_qqp"] * folded(m2_qqp, s, t, u, ft, fu)
    return tot


# ---------------------------------------------------------------- the PDF dump
def pdf_points(mass_bins):
    """every (x, Q) the quadrature will need.  Q = M_jj, the CMS central scale."""
    pts, index = [], {}
    for bi, (mlo, mhi) in enumerate(mass_bins):
        mn, mw = np.polynomial.legendre.leggauss(NM)
        mj = 0.5 * (mhi - mlo) * mn + 0.5 * (mhi + mlo)
        yn, yw = np.polynomial.legendre.leggauss(NYB)
        yb = YB_MAX * yn
        for i, m in enumerate(mj):
            tau = (m / ROOT_S) ** 2
            for j, y in enumerate(yb):
                x1 = math.sqrt(tau) * math.exp(y)
                x2 = math.sqrt(tau) * math.exp(-y)
                for k, x in ((0, x1), (1, x2)):
                    index[(bi, i, j, k)] = len(pts)
                    pts.append((x, m))
    return pts, index


def cache_paths(pts):
    """Name the cache after the POINTS it holds, not after the script that wrote it.

    Both this file and cijet_control.py used _pdf_in/_pdf_out.  They ask LHAPDF for different
    point sets -- 1280 at mu = M_jj here, 36864 at CIJET's scale there -- and the second run
    silently overwrote the first.  A later `--reuse` here then read the first 1280 rows of the
    OTHER file: right shape, wrong scale, wrong order, no error anywhere.  Keying the file by a
    digest of the request makes a stale cache impossible to read by accident, and read_dump()
    below refuses one whose digest does not match.

    EL CONJUNTO ENTRA EN LA CLAVE por la misma razon: desde que pdfdump.f lo acepta por stdin,
    dos corridas con los mismos puntos y distinto PDF ya no piden lo mismo, y una clave que solo
    mira los puntos volveria a servir el fichero de otro.
    [[a-guard-must-not-stand-on-a-public-name]]"""
    h = hashlib.sha1(("%s|%d|" % (PDFSET, len(pts))
                      + "|".join("%.10e,%.10e" % p for p in pts)).encode()).hexdigest()[:16]
    return HERE / ("_pdf_in_%s.txt" % h), HERE / ("_pdf_out_%s.txt" % h), h


def read_dump(pts):
    """the cached dump for exactly these points, or None."""
    _inp, outp, _h = cache_paths(pts)
    if not outp.exists():
        return None
    rows = [[float(v) for v in ln.split()] for ln in outp.read_text().splitlines()
            if len(ln.split()) == 14]
    if len(rows) != len(pts):
        return None
    return np.array(rows)


def run_dump(pts):
    inp, outp, _h = cache_paths(pts)
    with open(inp, "w") as fh:
        fh.write("%s\n" % PDFSET)          # pdfdump.f lee el conjunto en la primera linea
        fh.write("%d\n" % len(pts))
        for x, q in pts:
            fh.write("%.16e %.16e\n" % (x, q))
    cmd = ("gfortran -O2 -o pdfdump pdfdump.f $(lhapdf-config --libs) -lstdc++ "
           "&& ./pdfdump < %s > %s" % (inp.name, outp.name))
    args = ["docker", "run", "--rm", "-v", "%s:/w" % HERE.as_posix()]
    if PDFDIR.exists():
        args += ["-v", "%s:/pdfs" % PDFDIR.as_posix(),
                 "-e", "LHAPDF_DATA_PATH=/pdfs:/usr/local/share/LHAPDF"]
    args += ["-w", "//w", "cijet:1.0", "sh", "-c", cmd]
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("PDF dump failed:\n%s\n%s" % (r.stdout[-2000:], r.stderr[-2000:]))
    rows = []
    for ln in outp.read_text().splitlines():
        s = ln.split()
        if len(s) != 14:
            continue                      # the LHAPDF banner
        rows.append([float(v) for v in s])
    if len(rows) != len(pts):
        raise SystemExit("dump returned %d rows for %d points" % (len(rows), len(pts)))
    return np.array(rows)


def main():
    mass_bins, chi_bins = read_bins()
    print("=" * 100)
    print("THE TOWER IN chi_dijet AT HADRON LEVEL -- SPACELIKE CHANNELS ONLY")
    print("=" * 100)
    print("\n  binning read from %s" % CARD)
    print("  %d mass bins, %d chi bins, sqrt(s) = %.0f GeV, |y_boost| < %.2f"
          % (len(mass_bins), len(chi_bins), ROOT_S, YB_MAX))

    fails = []

    # ---- C1: the pT and rapidity cuts really are inactive at LO in this region ------------
    print("\n[C1] ARE THE pT AND RAPIDITY CUTS ACTIVE AT LO HERE?")
    worst_pt, worst_y = 1e9, 0.0
    for mlo, mhi in mass_bins:
        for clo, chi in chi_bins:
            for m in (mlo, mhi if mhi < 13000 else 8000.0):
                for c in (clo, chi):
                    pt = m * math.sqrt(c) / (1.0 + c)
                    worst_pt = min(worst_pt, pt)
                    worst_y = max(worst_y, YB_MAX + 0.5 * math.log(c))
    ok = worst_pt > 500.0 and worst_y < 2.5
    print("      lowest LO jet pT anywhere in the fiducial region : %.0f GeV (cut 500)" % worst_pt)
    print("      largest |y| reachable                            : %.3f  (cut 2.5)" % worst_y)
    print("\n      C1 %s -- the region is exactly |y_b|<1.11 and 1<chi<16"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C1")

    # ---- PDFs ------------------------------------------------------------------------------
    pts, index = pdf_points(mass_bins)
    print("\n  asking LHAPDF for %d (x, Q) points, no interpolation" % len(pts))
    # SI VINO DE CACHE O DE LHAPDF VA A stderr, no aqui.  Esa linea no es funcion de las entradas
    # sino del estado del disco: el mismo guion, con los mismos datos, imprimia "reusing the
    # cache" en una maquina y "dump returned 1280 rows" en otra, de modo que su salida archivada
    # no podia reproducirse nunca.  Lo destapo check_reproduces.py el dia que el conjunto de
    # partones entro en la clave del cache y el fichero viejo dejo de casar.  Es el mismo defecto
    # que los relojes de pared, con otra cara.  [[stale-outputs-lie]]
    rows = read_dump(pts)
    if rows is not None:
        sys.stderr.write("  reusing the cache keyed to these exact points\n")
    else:
        rows = run_dump(pts)
        sys.stderr.write("  dump returned %d rows\n" % len(rows))

    def xf(idx):
        r = rows[idx]
        return {f: r[f + 6] for f in range(-NF, NF + 1)}

    # ---- the distributions ------------------------------------------------------------------
    mn, mw = np.polynomial.legendre.leggauss(NM)
    yn, yw = np.polynomial.legendre.leggauss(NYB)
    cn, cw = np.polynomial.legendre.leggauss(NCHI)

    results = {}
    for tag, mkk in (("SM", None), ("ceiling", MKK["ceiling"]), ("escape", MKK["escape"])):
        table = np.zeros((len(mass_bins), len(chi_bins)))
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
                        f1 = xf(index[(bi, i, j, 0)])
                        f2 = xf(index[(bi, i, j, 1)])
                        # x f(x) -> the 1/(x1 x2) of the flux is folded in below
                        x1 = math.sqrt(s) / ROOT_S * math.exp(YB_MAX * yn[j])
                        x2 = math.sqrt(s) / ROOT_S * math.exp(-YB_MAX * yn[j])
                        d1 = {k: v / x1 for k, v in f1.items()}
                        d2 = {k: v / x2 for k, v in f2.items()}
                        for k, c in enumerate(cc):
                            w = weight(d1, d2, s, float(c), mkk)
                            # dsigma/dchi = |M|^2 / (16 pi s (1+chi)^2), times 2 M / s_had
                            acc += (mw[i] * jm * yw[j] * YB_MAX * cw[k] * jc
                                    * (2.0 * m / ROOT_S ** 2)
                                    * w / (16.0 * math.pi * s * (1.0 + c) ** 2))
                table[bi, ci] = acc
        results[tag] = table
        print("  computed %s" % tag)

    # ---- C2: pure QCD must come out nearly flat in chi (Rutherford) --------------------------
    print("\n[C2] THE QCD-ONLY chi DISTRIBUTION MUST BE NEARLY FLAT")
    sm = results["SM"]
    spread = []
    for bi, (mlo, mhi) in enumerate(mass_bins):
        dens = sm[bi] / np.array([c[1] - c[0] for c in chi_bins])
        dens = dens / dens.mean()
        spread.append(dens.max() - dens.min())
    print("      peak-to-peak of the normalised density, per mass bin:")
    print("      " + "  ".join("%.2f" % s for s in spread))
    ok = max(spread) < 0.9
    print("\n      C2 %s (worst %.2f; a Rutherford 1/t^2 gives a mild rise, not a factor)"
          % ("PASS" if ok else "FAIL", max(spread)))
    if not ok:
        fails.append("C2")

    # ---- the answer --------------------------------------------------------------------------
    print("\n[1] THE DEFORMATION, AS CMS WOULD SEE IT")
    print("""
    CMS normalises the chi distribution to unit area in each mass bin, so the observable is a
    RATIO OF SHAPES and the overall enhancement divides straight out.  Both numbers are given
    below because the difference between them is the single largest dilution in this file.""")
    payload = {"mass_bins": mass_bins, "chi_bins": chi_bins, "branches": {}}
    for tag in ("ceiling", "escape"):
        tw = results[tag]
        print("\n    1/R_5 = %.2f TeV" % (MKK[tag] / 1000.0))
        print("      %-16s %10s %10s %10s %10s"
              % ("mass bin [TeV]", "raw 1st", "norm 1st", "norm last", "max dev"))
        rows_out = []
        for bi, (mlo, mhi) in enumerate(mass_bins):
            raw = tw[bi] / sm[bi]
            nn = (tw[bi] / tw[bi].sum()) / (sm[bi] / sm[bi].sum())
            rows_out.append({"mass": [mlo, mhi], "raw": raw.tolist(), "norm": nn.tolist(),
                             "sm": sm[bi].tolist(), "tw": tw[bi].tolist()})
            print("      %5.1f-%-10.1f %10.3f %10.4f %10.4f %10.4f"
                  % (mlo / 1000, mhi / 1000, raw[0], nn[0], nn[-1],
                     float(np.abs(nn - 1).max())))
        payload["branches"][tag] = {"MKK_GeV": MKK[tag], "rows": rows_out}

    print("""
    Read the two middle columns together.  The raw column is the enhancement of the first chi
    bin over the Standard Model; the normalised one is what survives CMS's per-bin
    normalisation, and it is much smaller because the tower lifts the whole mass bin and only
    the SHAPE is measured.  Anything quoted against the published distribution has to be the
    normalised column.""")

    OUT.mkdir(exist_ok=True)
    (OUT / "hadronic_chi.json").write_text(json.dumps(payload, indent=1))
    print("\n    [wrote outputs/hadronic_chi.json]")

    print("\n" + "=" * 100)
    print("""STATUS.  This is a LEADING-ORDER hadronic prediction with the spacelike channels only.
        What it is NOT: a limit.  A limit needs the CMS covariance between chi bins and mass
        bins, which is on HEPData and was unreachable when this was written, and it needs the
        QCD scale choice profiled -- CMS's own small tension below 4.8 TeV moves between their
        two scale choices, so a fit that does not carry that as a nuisance is measuring the
        scale and not the tower.""")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        print("=" * 100)
        return 1
    print("VERDICT: controls pass.  The numbers above are the input to a fit, not its output.")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    sys.exit(main())
