#!/usr/bin/env python3
"""kk_dijet_lo -- the Kaluza-Klein gluon tower in the dijet angular distribution, at LO.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  An independent leading-order calculation of what the tower of [KM25] does to chi_dijet, built
  to be checked against the LO coefficients that CIJET 1.0 produces.  Two instruments, one
  number: if they disagree we are feeding CIJET the wrong operator, and that is exactly the
  step where a recast goes wrong silently.

THE IDEA, AND WHY IT IS BETTER THAN A CONTACT OPERATOR.  Because the quarks of [KM25] sit at a
fixed point, every mode of the colour tower couples to them with the same sqrt(2) g_s and with
the colour structure of the gluon itself.  So the tower is not a new vertex: it is the gluon
propagator, replaced.  Following [DMN01] eq. (7),

      D(q^2) = 1/q^2 + sum_{n>=1} c_n / (q^2 - m_n^2 + i m_n Gamma_n),   c_n = 2,

with m_n = n/R_5 and Gamma_n = 2 alpha_s m_n.  Writing the QCD matrix elements with D in place
of 1/q^2 gives the exact tree-level answer at any M_jj -- on the resonance as well as below it.
The contact-operator form (our Lambda_8) is only the leading term of D as q^2 -> 0, and
collapses precisely where the LHC has data, which is why it is used here only as a control.

WHAT IS AND IS NOT INCLUDED.  Tree level; the subprocesses in which the tower actually appears.
By [DMN01], five-momentum conservation forbids internal KK modes in any tree-level dijet diagram
with external gluons -- q qbar -> g g and g g -> anything are untouched -- so those are pure QCD
and are carried unmodified.  No PDFs here: this file is the partonic layer, which is where the
physics of the angular shape lives and where a comparison against CIJET is cleanest.

Run:  python kk_dijet_lo.py > outputs/kk_dijet_lo.txt
"""

import cmath
import math
import sys

NC = 3.0
CF = 4.0 / 3.0
NMAX = 200          # tower modes kept; convergence is a control below


# ==============================================================================================
# the effective propagator
# ==============================================================================================
def D_eff(q2, invR5, alpha_s, nmax=NMAX, tower=True, width=True):
    """[DMN01] eq. (7).  q2 and invR5 in GeV^2 / GeV.  Returns a complex number.

    width=False switches the tower widths off.  That is not physics -- it is the approximation
    under which the contact coefficient Lambda_8 of the paper is derived, and control C4 needs
    it to compare like with like.  C4b then measures what the width actually costs.
    """
    d = complex(1.0 / q2, 0.0)
    if not tower:
        return d
    for n in range(1, nmax + 1):
        mn = n * invR5
        mn2 = mn * mn
        gam = 2.0 * alpha_s * mn if width else 0.0    # Gamma_n = 2 alpha_s m_n, [DMN01]
        d += 2.0 / complex(q2 - mn2, mn * gam)
    # analytic tail.  Every mode above nmax is far off shell, so 2/(q2 - m_n^2) -> -2 R^2/n^2
    # and the remainder of the sum is zeta(2) minus the part already taken.  Without this the
    # truncation leaves a ~1/nmax bias, which is 0.3% at nmax=200 -- small, and exactly the size
    # of the effect C4b is trying to measure, so it is not allowed to stay in.
    tail = math.pi ** 2 / 6.0 - sum(1.0 / n ** 2 for n in range(1, nmax + 1))
    d += complex(-2.0 * tail / invR5 ** 2, 0.0)
    return d


# ==============================================================================================
# LO squared matrix elements, summed/averaged over spin and colour, in units of g_s^4.
# Each channel's 1/q^2 is carried explicitly so the tower can be inserted.
# Reference forms: Ellis-Stirling-Webber table 7.1, rewritten with explicit propagators.
# ==============================================================================================
def m2_qqp(s, t, u, Dt):
    """q q' -> q q'  (distinct flavours): t-channel only."""
    return (4.0 / 9.0) * (s * s + u * u) * abs(Dt) ** 2


def m2_qqbarp(s, t, u, Dt):
    """q qbar' -> q qbar'  (distinct flavours): t-channel only."""
    return (4.0 / 9.0) * (s * s + u * u) * abs(Dt) ** 2


def m2_qq(s, t, u, Dt, Du):
    """q q -> q q  (identical flavours): t and u channels, with their interference."""
    return ((4.0 / 9.0) * ((s * s + u * u) * abs(Dt) ** 2 + (s * s + t * t) * abs(Du) ** 2)
            - (8.0 / 27.0) * s * s * (Dt * Du.conjugate()).real)


def m2_qqbar_same(s, t, u, Dt, Ds):
    """q qbar -> q qbar  (same flavour): t and s channels, with their interference."""
    return ((4.0 / 9.0) * ((s * s + u * u) * abs(Dt) ** 2 + (u * u + t * t) * abs(Ds) ** 2)
            - (8.0 / 27.0) * u * u * (Dt * Ds.conjugate()).real)


def m2_qqbar_ann(s, t, u, Ds):
    """q qbar -> q' qbar'  (annihilation into a different flavour): s-channel only."""
    return (4.0 / 9.0) * (t * t + u * u) * abs(Ds) ** 2


# ==============================================================================================
def kinematics(shat, chi):
    """Massless 2->2 in the partonic CM, parametrised by chi = exp(|y1-y2|)."""
    t = -shat / (1.0 + chi)
    u = -shat * chi / (1.0 + chi)
    return t, u


def line(c='-', n=94):
    print(c * n)


def dsig_dchi(shat, chi, invR5, alpha_s, tower=True):
    """Sum of the tower-sensitive channels, times the chi Jacobian.  Arbitrary common units."""
    t, u = kinematics(shat, chi)
    Dt = D_eff(t, invR5, alpha_s, tower=tower)
    Du = D_eff(u, invR5, alpha_s, tower=tower)
    Ds = D_eff(shat, invR5, alpha_s, tower=tower)
    tot = (m2_qqp(shat, t, u, Dt)
           + m2_qqbarp(shat, t, u, Dt)
           + m2_qq(shat, t, u, Dt, Du)
           + m2_qqbar_same(shat, t, u, Dt, Ds)
           + m2_qqbar_ann(shat, t, u, Ds))
    return tot * 2.0 / (1.0 + chi) ** 2      # dcos(theta*)/dchi


# ==============================================================================================
def main():
    fails = []
    line('=')
    print("THE KALUZA-KLEIN GLUON TOWER IN chi_dijet, AT LEADING ORDER")
    print("independent of CIJET, and built to be compared with it")
    line('=')

    ALPHA_S = 0.08
    CHIS = [1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]

    # ------------------------------------------------------------------ controls
    print("\n[0] CONTROLS")

    # C1: with the tower switched off, D must be exactly the gluon propagator
    d = D_eff(-1.0e6, 4000.0, ALPHA_S, tower=False)
    ok = abs(d - complex(-1e-6, 0)) < 1e-18
    print("   C1  tower off  =>  D(q2) = 1/q2 exactly ................... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C1")

    # C2: QCD alone must give a nearly FLAT chi distribution -- the textbook Rutherford result
    qcd = [dsig_dchi(3.0e6, c, 4000.0, ALPHA_S, tower=False) for c in CHIS]
    spread = (max(qcd) - min(qcd)) / (sum(qcd) / len(qcd))
    ok = spread < 0.45
    print("   C2  QCD-only chi distribution is flat to %4.1f%% ........... %s"
          % (100 * spread, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C2")
    print("       (the Rutherford 1/t^2 behaviour; a broken matrix element loses this at once)")

    # C3: the tower sum must have converged in NMAX
    a = D_eff(-1.0e6, 4000.0, ALPHA_S, nmax=NMAX)
    b = D_eff(-1.0e6, 4000.0, ALPHA_S, nmax=2 * NMAX)
    rel = abs(a - b) / abs(a)
    ok = rel < 1e-3
    print("   C3  tower converged: |D(200) - D(400)| / |D| = %.2e ...... %s"
          % (rel, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3")

    # C4: in the deep EFT region the tower must reproduce our Lambda_8 contact coefficient
    #     D(q2) -> 1/q2 - 2 R^2 zeta(2) = 1/q2 - (pi^2/3) R^2
    invR5 = 9090.0
    q2 = -(300.0 ** 2)                      # |q2| << (1/R5)^2
    R2 = 1.0 / invR5 ** 2
    want = -(math.pi ** 2 / 3.0) * R2
    got0 = (D_eff(q2, invR5, ALPHA_S, width=False) - 1.0 / q2).real
    rel = abs(got0 - want) / abs(want)
    ok = rel < 2e-3
    print("   C4  deep-EFT limit reproduces -(pi^2/3) R_5^2 ............ %s (rel %.1e)"
          % ("PASS" if ok else "FAIL", rel))
    if not ok:
        fails.append("C4")
    print("       this is the SAME pi^2/6 that colour_octet_recast.py and the paper carry,")
    print("       so the two files are now tied together by a number and not by a claim.")
    print("       NOTE the comparison is made at zero width, because that is the approximation")
    print("       Lambda_8 is derived under.  C4b is what that approximation costs.")

    # C4b: the finite width is not decoration -- it shifts the contact coefficient
    gotw = (D_eff(q2, invR5, ALPHA_S, width=True) - 1.0 / q2).real
    shift = (gotw - want) / abs(want)
    # a finite width REDUCES the magnitude of each off-shell mode by 1/(1+(Gamma/m)^2), and
    # the coefficient is negative, so the shift relative to |want| comes out POSITIVE.
    r2 = (2.0 * ALPHA_S) ** 2
    expected = r2 / (1.0 + r2)
    ok = abs(shift - expected) < 2e-3
    print("   C4b width shifts the contact coefficient by %+.2f%% ....... %s"
          % (100 * shift, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C4b")
    print("       (Gamma_n/m_n)^2 = (2 alpha_s)^2 = %.4f, so 1/(1+r^2) predicts %+.2f%%."
          % (r2, 100 * expected))
    print("       So the Lambda_8 quoted in the paper is a ZERO-WIDTH number: the tower's true")
    print("       contact coefficient is ~%.0f%% smaller, i.e. Lambda_8 is ~%.1f%% HIGHER."
          % (100 * shift, 50 * shift))

    # F1: a falsification -- with the tower on, the distribution must NOT stay flat
    kk = [dsig_dchi(3.0e6, c, 4000.0, ALPHA_S) for c in CHIS]
    r = [k / q for k, q in zip(kk, qcd)]
    varies = (max(r) - min(r)) / (sum(r) / len(r)) > 0.02
    print("   F1  tower ON does change the shape (else nothing is tested)  %s"
          % ("PASS" if varies else "FAIL"))
    if not varies:
        fails.append("F1")

    # ------------------------------------------------------------------ the physics
    line('=')
    print("[1] RATIO (QCD + tower) / QCD IN chi, AT FIXED DIJET MASS")
    line('=')
    print("\n    The tower is s-channel-like and therefore CENTRAL: it populates low chi, where")
    print("    t-channel QCD is weakest.  That is the discriminant the angular analysis uses.\n")

    MJJ = (2.7, 3.3, 3.9, 4.5, 5.1, 6.5, 9.0)
    record = {"what": "ratio (QCD+KK tower)/QCD in dchi at fixed M_jj, LO, partonic",
              "alpha_s": ALPHA_S, "nmax": NMAX, "chi": CHIS, "mjj_TeV": list(MJJ),
              "propagator": "DMN01 eq.(7): 1/q2 + sum_n 2/(q2 - m_n^2 + i m_n Gamma_n)",
              "Gamma_n": "2 alpha_s m_n", "m_n": "n / R_5", "ratios": {}}

    for invR5 in (3970.0, 9090.0):
        print("    1/R_5 = %.2f TeV" % (invR5 / 1000.0))
        print("      %-10s" % "M_jj[TeV]" + "".join("%8.1f" % c for c in CHIS))
        line('-', 94)
        rows = []
        for mjj in MJJ:
            shat = (mjj * 1000.0) ** 2
            row = []
            for c in CHIS:
                a = dsig_dchi(shat, c, invR5, ALPHA_S)
                b = dsig_dchi(shat, c, invR5, ALPHA_S, tower=False)
                row.append(a / b)
            print("      %-10.1f" % mjj + "".join("%8.3f" % v for v in row))
            rows.append([round(v, 6) for v in row])
        record["ratios"]["invR5_GeV_%d" % int(invR5)] = rows
        line('-', 94)
        print()

    # the numbers a paper would have to carry, written where a gate can find them
    import json
    import os
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
    with open(os.path.join(outdir, 'kk_dijet_lo.json'), 'w') as fh:
        json.dump(record, fh, indent=1, sort_keys=True)
    print("    [wrote outputs/kk_dijet_lo.json -- every number above, for the paper and the gates]")
    print()

    print("    Read the first column: at low chi the ratio departs from 1 by tens of percent")
    print("    well below the resonance, and the departure grows towards it.  CMS quote a total")
    print("    theoretical uncertainty of 3.5% in the 3.0-3.6 TeV bin and 12.2% above 7 TeV")
    print("    (their Table 2), so an effect of this size is not something the analysis could")
    print("    have absorbed -- which is the whole reason the recast is worth finishing.")

    line('=')
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line('=')
        return 1
    print("VERDICT: controls pass (C1-C4), falsification control behaves (F1).")
    print("STATUS: partonic layer only.  PDFs, NLO and the CMS covariance are still missing,")
    print("        and until CIJET's LO agrees with this, neither number should be quoted.")
    line('=')
    return 0


if __name__ == '__main__':
    sys.exit(main())
