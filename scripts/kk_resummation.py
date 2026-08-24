#!/usr/bin/env python3
"""kk_resummation -- the whole coloured tower, in closed form, and what that settles.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

kk_dijet_lo.py evaluates the replaced propagator of [DMN01] by TRUNCATING the tower at n = 400
and carrying a Breit-Wigner width on every mode.  In the spacelike channels neither is necessary,
and one of them is wrong.

THE RESUMMATION.  For t < 0 the propagator never reaches a pole, so the width may be dropped and
the sum done exactly.  With Q = sqrt(-t) and a = Q R_5, Euler's

    sum_{n>=1} 1/(n^2 + a^2)  =  [ pi a coth(pi a) - 1 ] / (2 a^2)

turns the replaced propagator into a form factor on the gluon one:

    1/t + sum_{n>=1} 2/(t - n^2/R_5^2)  =  F(t)/t ,     F(t) = pi R_5 sqrt(-t) coth(pi R_5 sqrt(-t)).

For a subprocess with only t-channel gluon exchange (q q' -> q q', distinct flavours) every
QCD amplitude is multiplied by F, so the parton-level ratio is F^2 exactly -- no truncation, no
width, and NO dependence on BR(G_n -> jj), because nothing is being produced and decayed.

WHAT IS OURS AND WHAT IS NOT.  The identity is Euler's and resumming a Kaluza-Klein tower into a
coth is standard in extra-dimensional phenomenology; it is how a brane-to-brane propagator is
normally written.  What is specific here is that it APPLIES: every mode reaches the brane quarks
with the same sqrt(2), by the parity theorem and [KM25]'s localisation of the quarks at the fixed
point, so the coefficients c_n are all equal and the sum really is sum 1/(n^2 + a^2).  With a
mode-dependent coupling there would be no closed form.

WHAT IT SETTLES.  kk_dijet_lo.py's control C4b reports that restoring Gamma_n = 2 alpha_s m_n
changes the contact coefficient by +2.55 %, and the paper carries that as a correction to
Lambda_8.  It is not a correction, it is an artefact: a Kaluza-Klein gluon below its two-parton
threshold has NO absorptive part, so Im Pi = 0 and the spacelike propagator is real.  C4 below
makes that measurable rather than asserted -- a genuine self-energy correction must vanish as
q^2 -> 0^-, and the fixed-width prescription's does not; it tends to a constant.

Run:  python kk_resummation.py > outputs/kk_resummation.txt
"""
import cmath
import json
import math
import pathlib
import sys

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"

ALPHA_S = 0.0789          # at the escape branch, as kk_dijet_lo.py uses it


def tower_sum(t, invR5, nmax, width=False, alpha_s=ALPHA_S):
    """1/t + sum_{n=1..nmax} 2/(t - m_n^2 + i m_n Gamma_n), by brute force."""
    tot = complex(1.0 / t, 0.0)
    for n in range(1, nmax + 1):
        mn = n * invR5
        gam = 2.0 * alpha_s * mn if width else 0.0
        tot += 2.0 / complex(t - mn * mn, mn * gam)
    return tot


def form_factor(t, invR5):
    """F(t) = pi R_5 sqrt(-t) coth(pi R_5 sqrt(-t)), the exact resummation for t < 0."""
    a = math.pi * math.sqrt(-t) / invR5
    return a / math.tanh(a)


def closed(t, invR5):
    return form_factor(t, invR5) / t


def line(c="-", n=94):
    print(c * n)


def main():
    fails = []
    line("=")
    print("THE COLOURED TOWER RESUMMED: A FORM FACTOR, AND WHAT THAT SETTLES ABOUT THE WIDTH")
    line("=")

    invR5 = 3.97          # TeV
    print("\n[1] THE CLOSED FORM AGAINST THE SUM ITSELF,  1/R_5 = %.2f TeV" % invR5)
    print("""
    The comparison is against the INFINITE sum evaluated to thirty digits, not against a
    truncation: a truncated tower leaves a tail of size 2 R_5^2 / N, which at N = 20000 is
    6.3e-06 -- constant in t, and therefore a growing RELATIVE error exactly where the closed
    form is smallest.  The first draft of this control compared against the truncation and
    read that tail as a failure of the resummation.  The truncated column is kept below so the
    size of what truncation costs is visible.""")
    mp.mp.dps = 30
    print("\n    %10s %8s %18s %18s %12s %10s"
          % ("-t [TeV^2]", "a", "exact sum (30 dps)", "pi a coth(pi a)/t", "N=20000", "rel"))
    line()
    worst = 0.0
    for mt in (0.01, 0.1, 1.0, 4.0, 9.0, 25.0, 100.0, 400.0):
        t = -mt
        ex = mp.mpf(1) / t + mp.nsum(lambda n: 2 / (t - (n * mp.mpf(invR5)) ** 2), [1, mp.inf])
        cf = closed(t, invR5)
        brute = tower_sum(t, invR5, 20000).real
        rel = abs(float(ex) - cf) / abs(cf)
        worst = max(worst, rel)
        a = math.pi * math.sqrt(mt) / invR5
        print("    %10.2f %8.4f %18.12f %18.12f %12.6f %10.1e"
              % (mt, a, float(ex), cf, brute, rel))
    line()
    ok = worst < 1e-12
    print("\n   C1  the resummation IS the sum                      %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst))
    if not ok:
        fails.append("C1")

    print("\n    F1  THE CONTROL THAT MUST FAIL: the same test with cot in place of coth,")
    bad = 0.0
    for mt in (1.0, 9.0, 100.0):
        t = -mt
        a = math.pi * math.sqrt(mt) / invR5
        wrong = (a / math.tan(a)) / t
        bad = max(bad, abs(tower_sum(t, invR5, 20000).real - wrong) / abs(wrong))
    print("        worst relative difference : %.2e   (must be large)  %s"
          % (bad, "ok" if bad > 1e-2 else "THE TEST IS BLIND"))
    if bad <= 1e-2:
        fails.append("F1")

    # ---------------------------------------------------------------------------------------
    print("\n[2] THE EFT LIMIT IS THE PAPER'S OWN COEFFICIENT")
    print("    pi a coth(pi a) = 1 + (pi a)^2/3 - (pi a)^4/45 + ... , so D - 1/t -> -pi^2 R_5^2/3")
    want = -math.pi ** 2 / (3.0 * invR5 ** 2)
    print("\n    %12s %18s %18s %10s" % ("-t [TeV^2]", "D(t) - 1/t", "-pi^2 R_5^2/3", "rel"))
    line()
    worst2 = 0.0
    for mt in (1e-4, 1e-3, 1e-2):
        t = -mt
        got = closed(t, invR5) - 1.0 / t
        rel = abs(got - want) / abs(want)
        worst2 = max(worst2, rel)
        print("    %12.0e %18.10f %18.10f %10.1e" % (mt, got, want, rel))
    line()
    ok = worst2 < 1e-3
    print("\n   C2  and it is the pi^2/6 of eq. (tower)                %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst2))
    if not ok:
        fails.append("C2")

    # ---------------------------------------------------------------------------------------
    print("\n[3] THE WIDTH AT SPACELIKE MOMENTUM IS A PRESCRIPTION, NOT A CORRECTION")
    print("""
    A Kaluza-Klein gluon at t < 0 sits below its two-parton threshold, so its self-energy has
    no absorptive part and the propagator is REAL.  A genuine correction of that kind must
    therefore switch off as the momentum leaves the resonance.  A constant i m_n Gamma_n does
    not.  Measured, rather than argued:""")
    print("\n    %12s %20s %20s" % ("-t [TeV^2]", "with width / no width", "1/(1+(2 alpha_s)^2)"))
    line()
    pred = 1.0 / (1.0 + (2.0 * ALPHA_S) ** 2)
    shifts = []
    for mt in (1e-6, 1e-4, 1e-2, 1.0):
        t = -mt
        w = (tower_sum(t, invR5, 4000, width=True) - 1.0 / t).real
        z = (tower_sum(t, invR5, 4000, width=False) - 1.0 / t).real
        shifts.append(w / z)
        print("    %12.0e %20.6f %20.6f" % (mt, w / z, pred))
    line()
    # the shift must NOT die away as t -> 0: that is exactly what makes it unphysical there.
    drift = abs(shifts[0] - shifts[-1])
    persists = abs(shifts[0] - pred) < 5e-3
    ok = persists and drift < 5e-2
    print("\n    at t -> 0 the prescription still shifts the coefficient by %+.2f %%, and an"
          % (100.0 * (shifts[0] - 1.0)))
    print("    absorptive part that survives below threshold is not one.")
    print("\n   C3  the shift persists to t = 0, so it is an artefact   %s"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3")
    print("""
    CONSEQUENCE, and it is a correction to this paper.  The 'second correction' of
    \\S collider -- Lambda_8 low by about 1.3 %, from restoring Gamma_n -- is withdrawn.  In the
    channel Lambda_8 describes there is no width to restore, and the zero-width coefficient
    pi^2/6 is the right one, not an approximation to something else.""")

    # ---------------------------------------------------------------------------------------
    print("\n[3b] AND THE TIMELIKE SIDE IS THE SAME FUNCTION")
    print("""
    Nothing about the sum required t < 0; that was only where the width could be dropped.  For
    q^2 > 0 and away from the poles the same Euler expansion gives, by coth(i x) = -i cot x,

        1/s + sum_n 2/(s - n^2/R_5^2)  =  (1/s) pi b cot(pi b) ,   b = R_5 sqrt(s) ,

    so ONE analytic function covers both:  F(q^2) = pi a coth(pi a) for q^2 < 0 and
    pi b cot(pi b) for q^2 > 0.  The cotangent's poles at b = 1, 2, 3, ... ARE the Kaluza-Klein
    resonances -- they are not put in, they are where the closed form diverges -- and it is
    exactly there, and only there, that a width has to replace the pole.

    Which makes the three statements of the collider section one statement: the contact
    operator is the small-b expansion, the angular distortion is the spacelike branch, and the
    resonances are the poles of the timelike branch.""")
    print("\n    %12s %10s %20s %20s %10s"
          % ("s [TeV^2]", "b", "truncated (N=20000)", "pi b cot(pi b)/s", "rel"))
    line()
    worst3 = 0.0
    for ss in (0.25, 1.0, 4.0, 9.0, 30.0, 60.0):
        b = math.sqrt(ss) / invR5
        if abs(b - round(b)) < 0.02:          # keep clear of the poles themselves
            continue
        brute = tower_sum(ss, invR5, 20000).real
        cf = (math.pi * b / math.tan(math.pi * b)) / ss
        rel = abs(brute - cf) / abs(cf)
        worst3 = max(worst3, rel)
        print("    %12.2f %10.4f %20.9f %20.9f %10.1e" % (ss, b, brute, cf, rel))
    line()
    # the truncation tail is the same 2 R_5^2 / N as on the spacelike side, so the bar is the
    # absolute one, not a relative one.  [[a-constant-difference-is-an-object-mismatch]]
    ok = worst3 < 3e-3
    print("\n   C3b the cotangent branch reproduces the tower too  %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst3))
    if not ok:
        fails.append("C3b")

    print("\n    and the poles are the tower, not an input:")
    for n in (1, 2, 3):
        s_pole = (n * invR5) ** 2
        b = math.sqrt(s_pole * 0.9999) / invR5
        val = abs((math.pi * b / math.tan(math.pi * b)) / (s_pole * 0.9999))
        print("       n = %d :  sqrt(s) = %.2f TeV = %d/R_5, |F/s| just below it = %.3e"
              % (n, math.sqrt(s_pole), n, val))

    # ---------------------------------------------------------------------------------------
    print("\n[3c] THE SECOND CIRCLE IS UV-SENSITIVE, AND THE HONEST STATEMENT IS ITS COEFFICIENT")
    print("""
    An earlier version of the collider section said the logarithmic divergence a codimension-two
    brane would give "never switches on".  That is too strong and a reader who knows
    brane-to-brane propagators would object: a two-dimensional tower IS logarithmically
    sensitive to the cut-off, and no hierarchy makes a divergence go away.  What a hierarchy
    does is suppress its COEFFICIENT, and that is computable.

    Summing n at fixed m first, with the same Euler identity,

        sum_{n in Z} 1/(n^2/R_5^2 + m^2/R_6^2) = pi R_5 (R_6/m) coth(pi m R_5/R_6)
                                               -> pi R_5 R_6 / m      for R_5 >> R_6,

    so the m >= 1 modes contribute pi R_5 R_6 sum_{m<=Lambda R_6} 1/m = pi R_5 R_6 ln(Lambda R_6),
    against the m = 0 tower's pi^2 R_5^2/6.  Hence

        Delta C(m>=1) / C(m=0)  =  (6/pi) (R_6/R_5) ln(Lambda R_6) ,

    logarithmically UV-sensitive, with a coefficient carrying one power of R_6/R_5.""")
    print("\n    %8s %12s %18s %18s %10s"
          % ("R5/R6", "ln(L R_6)", "2D sum, direct", "(6/pi)(R6/R5)lnL", "rel"))
    line()
    worst4 = 0.0
    for ratio in (20.0, 50.0, 200.0, 1000.0):
        R5 = 1.0
        R6 = R5 / ratio
        mmax = 400                      # the cut-off, in units of m
        direct = 0.0
        for m in range(1, mmax + 1):
            am = math.pi * m * R5 / R6
            # sum over all integer n, including n = 0
            direct += math.pi * R5 * (R6 / m) / math.tanh(am) if am < 700 else \
                math.pi * R5 * (R6 / m)
        c0 = math.pi ** 2 * R5 ** 2 / 6.0
        pred = (6.0 / math.pi) * (R6 / R5) * (math.log(mmax) + 0.5772156649)
        got = direct / c0
        rel = abs(got - pred) / pred
        worst4 = max(worst4, rel)
        print("    %8.0f %12.4f %18.6e %18.6e %10.1e"
              % (ratio, math.log(mmax), got, pred, rel))
    line()
    # the harmonic sum is ln(M) + gamma + O(1/M), so the prediction above carries gamma; what is
    # being tested is the R_6/R_5 power and the 6/pi, not the constant under the log.
    ok = worst4 < 5e-3
    print("\n   C3c the suppression is (6/pi)(R_6/R_5) ln, not zero  %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst4))
    if not ok:
        fails.append("C3c")
    print("""
    At [KM25]'s own hierarchy the number is small -- the table above gives 1.3 % at
    R_5/R_6 = 10^3, and their p. 3 contemplates 1/R_6 at the Planck scale, where it is far
    smaller still -- but it is small BECAUSE of the hierarchy, and that is what the paper should
    say.  A divergence that is suppressed is not a divergence that is absent.""")

    # ---------------------------------------------------------------------------------------
    print("\n[4] THE OBSERVABLE, IN CLOSED FORM")
    print("""
    For massless jets |t| = M_jj^2/(1+chi), so on the pure t-channel subprocess q q' -> q q'

        dsigma(QCD+tower)/dsigma(QCD)  =  F^2 ,   F = pi M_jj / (M_KK sqrt(1+chi))  coth( same ).

    Independent of BR(G_n -> jj): nothing is produced and decayed.  Quoted at chi = 1.5, the
    centre of the first angular bin of the CMS measurement.""")
    rows = {}
    print("\n    %10s %14s %14s" % ("M_jj [TeV]", "M_KK=9.09", "M_KK=3.97"))
    line()
    for mjj in (3.3, 3.9, 4.5, 5.1, 6.5):
        vals = []
        for mkk in (9.09, 3.97):
            t = -(mjj ** 2) / (1.0 + 1.5)
            vals.append(form_factor(t, mkk) ** 2)
        rows["%.1f" % mjj] = {"9.09": vals[0], "3.97": vals[1]}
        print("    %10.1f %14.3f %14.3f" % (mjj, vals[0], vals[1]))
    line()
    print("""
    These are ONE subprocess at parton level.  The inclusive CMS sample also carries q g and
    g g, which do not see the tower at all, and every distribution is normalised per mass bin.
    Both dilute these numbers.  They say where to look; they are not a prediction of the
    measured ratio, and nothing here is a limit.""")

    OUT.mkdir(exist_ok=True)
    (OUT / "kk_resummation.json").write_text(json.dumps(
        {"form_factor": "pi R_5 sqrt(-t) coth(pi R_5 sqrt(-t))",
         "eft_coefficient": want, "invR5_used": invR5,
         "width_shift_at_zero": shifts[0], "ratios_chi_1p5": rows}, indent=1))

    line("=")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line("=")
        return 1
    print("VERDICT: the tower resums exactly in the spacelike channels; the EFT coefficient is")
    print("         recovered; and the width correction the paper carried there is withdrawn.")
    line("=")
    return 0


if __name__ == "__main__":
    sys.exit(main())
