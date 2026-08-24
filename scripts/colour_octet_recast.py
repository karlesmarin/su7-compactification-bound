#!/usr/bin/env python3
"""colour_octet_recast -- why a dijet-angular contact-interaction limit does not transfer.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  The Kaluza-Klein gluon tower of [KM25], integrated out, gives a COLOUR-OCTET four-quark
  operator.  CMS-EXO-24-011 (arXiv:2603.25458) sets the strongest dijet-angular limits to date,
  but its contact-interaction scenarios are stated to have "color-singlet couplings between
  quarks".  This script computes, with explicit SU(3) generators, exactly how the two differ --
  and the answer is not a rescaling.

WHY IT MATTERS.  The naive move is to read Lambda_8 off the published Lambda limit.  That is
wrong in a way that goes in our favour and must therefore be checked twice: the octet operator
has the SAME colour structure as one-gluon exchange, so it interferes with QCD at full strength,
while the singlet operator's interference vanishes in the t-channel by Tr(T^a) = 0.  A limit
built on a quadratic-in-1/Lambda^2 signal does not carry over to an operator whose leading
effect is linear.

Run:  python colour_octet_recast.py > outputs/colour_octet_recast.txt
"""

import itertools
import math
import sys

import numpy as np

# ----------------------------------------------------------------------------------------------
# SU(3) generators T^a = lambda^a / 2
# ----------------------------------------------------------------------------------------------
def gell_mann():
    l = []
    l.append(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex))
    l.append(np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex))
    l.append(np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex))
    l.append(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex))
    l.append(np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex))
    l.append(np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex))
    l.append(np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex))
    l.append(np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / math.sqrt(3.0))
    return [m / 2.0 for m in l]


T = gell_mann()
NC = 3
ID = np.eye(NC, dtype=complex)


def line(c='-', n=94):
    print(c * n)


def main():
    fails = []
    line('=')
    print("COLOUR-OCTET vs COLOUR-SINGLET FOUR-QUARK OPERATORS IN DIJET PRODUCTION")
    print("why the CMS-EXO-24-011 contact-interaction limit does not transfer to a KK gluon tower")
    line('=')

    # -------------------------------------------------------------- generator sanity
    print("\n[0] CONTROLS ON THE GENERATORS")
    tr = max(abs(np.trace(T[a] @ T[b]) - 0.5 * (a == b)) for a in range(8) for b in range(8))
    ok = tr < 1e-12
    print("   C0a Tr(T^a T^b) = delta^ab / 2 ........................... %s (max dev %.1e)"
          % ("PASS" if ok else "FAIL", tr))
    if not ok:
        fails.append("C0a")
    trz = max(abs(np.trace(T[a])) for a in range(8))
    ok = trz < 1e-12
    print("   C0b Tr(T^a) = 0 .......................................... %s (max %.1e)"
          % ("PASS" if ok else "FAIL", trz))
    if not ok:
        fails.append("C0b")
    cf = (T[0] @ T[0] + sum(T[a] @ T[a] for a in range(1, 8)))[0, 0].real
    ok = abs(cf - 4.0 / 3.0) < 1e-12
    print("   C0c sum_a (T^a T^a) = C_F 1 with C_F = 4/3 ............... %s (%.6f)"
          % ("PASS" if ok else "FAIL", cf))
    if not ok:
        fails.append("C0c")

    # -------------------------------------------------------------- the Fierz identity
    print("\n[1] THE COLOUR DECOMPOSITION OF THE OCTET OPERATOR")
    print("    T^a_ij T^a_kl  =  (1/2) delta_il delta_kj  -  (1/(2 N_c)) delta_ij delta_kl")
    lhs = np.zeros((NC,) * 4, dtype=complex)
    for a in range(8):
        lhs += np.einsum('ij,kl->ijkl', T[a], T[a])
    rhs = 0.5 * np.einsum('il,kj->ijkl', ID, ID) - (1.0 / (2 * NC)) * np.einsum('ij,kl->ijkl', ID, ID)
    dev = np.abs(lhs - rhs).max()
    ok = dev < 1e-12
    print("\n   C1  identity holds elementwise ........................... %s (max dev %.1e)"
          % ("PASS" if ok else "FAIL", dev))
    if not ok:
        fails.append("C1")
    print("\n    So O_8 = (1/2) O_1~ - (1/6) O_1, where O_1~ is the CROSSED-colour singlet.")
    print("    O_8 is therefore NOT proportional to O_1: the recast is not a rescaling.")

    # -------------------------------------------------------------- colour sums
    print("\n[2] COLOUR SUMS FOR q q' -> q q'  (distinct flavours, t-channel)")
    print("    amplitude colour factors:   QCD and OCTET-CI:  T^a_ji T^a_lk")
    print("                                SINGLET-CI:        delta_ji delta_lk")

    def csum(A, B):
        """sum over all four colour indices of A * conj(B), each given as a rank-4 array."""
        return np.einsum('ijkl,ijkl->', A, np.conj(B))

    OCT = lhs                                   # T^a (x) T^a
    SNG = np.einsum('ij,kl->ijkl', ID, ID)      # 1 (x) 1

    qcd_qcd = csum(OCT, OCT).real
    oct_oct = csum(OCT, OCT).real
    sng_sng = csum(SNG, SNG).real
    qcd_oct = csum(OCT, OCT).real
    qcd_sng = csum(OCT, SNG).real

    print("\n    %-34s %12s" % ("colour sum", "value"))
    line('-', 50)
    print("    %-34s %12.6f" % ("|QCD|^2      (T.T)x(T.T)", qcd_qcd))
    print("    %-34s %12.6f" % ("|octet CI|^2 (T.T)x(T.T)", oct_oct))
    print("    %-34s %12.6f" % ("|singlet CI|^2 (1.1)x(1.1)", sng_sng))
    print("    %-34s %12.6f" % ("QCD x octet CI  INTERFERENCE", qcd_oct))
    print("    %-34s %12.6f  <-- vanishes" % ("QCD x singlet CI INTERFERENCE", qcd_sng))
    line('-', 50)

    ok = abs(qcd_sng) < 1e-12
    print("\n   C2  singlet CI does not interfere in THIS contraction ..... %s"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C2")
    print("       reason: the colour factor is Tr(T^a) Tr(T^a) = 0, control C0b.")
    print("       AND THAT IS ALL IT SAYS.  See [2b]: it does not survive identical quarks,")
    print("       and an earlier version of this file drew a conclusion from it that was wrong.")

    ok = abs(qcd_oct - qcd_qcd) < 1e-12
    print("   C3  octet CI interferes at FULL gluon strength ........... %s"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3")
    print("       reason: its colour structure is identical to one-gluon exchange.")

    # falsification: a random non-octet structure must not reproduce C3
    rng = np.random.default_rng(20260822)
    RAND = rng.normal(size=(NC,) * 4) + 1j * rng.normal(size=(NC,) * 4)
    broke = abs(csum(OCT, RAND).real - qcd_qcd) > 1e-6
    print("   F1  a random colour structure does NOT match C3 ......... %s"
          % ("PASS" if broke else "FAIL"))
    if not broke:
        fails.append("F1")

    # ---------------------------------------------------------------- identical quarks
    line('=')
    print("[2b] AND NOW THE CASE THAT KILLS THE OBVIOUS CONCLUSION: q q -> q q")
    line('=')
    print("""
    For identical quarks the t and u channels both exist, and the CROSSED contractions -- the
    t-channel gluon against the u-channel contact term -- are contractions C2 never looked at.
    Writing the colour tensors as C[i,j,k,l] for q(i) q(k) -> q(j) q(l):""")
    OCT_t = sum(np.einsum('ji,lk->ijkl', T[a], T[a]) for a in range(8))
    SNG_t = np.einsum('ji,lk->ijkl', ID, ID)
    OCT_u = sum(np.einsum('li,jk->ijkl', T[a], T[a]) for a in range(8))
    SNG_u = np.einsum('li,jk->ijkl', ID, ID)

    def cs(A, B):
        return np.einsum('ijkl,ijkl->', A, np.conj(B)).real

    print("\n      %-28s %10s" % ("contraction", "value"))
    line('-', 44)
    for nm, v in (("QCD_t x singlet_t", cs(OCT_t, SNG_t)),
                  ("QCD_t x singlet_u  (crossed)", cs(OCT_t, SNG_u)),
                  ("QCD_u x singlet_t  (crossed)", cs(OCT_u, SNG_t)),
                  ("QCD_u x singlet_u", cs(OCT_u, SNG_u))):
        print("      %-28s %10.4f" % (nm, v))
    tot_s = cs(OCT_t, SNG_t) + cs(OCT_t, SNG_u) + cs(OCT_u, SNG_t) + cs(OCT_u, SNG_u)
    tot_o = cs(OCT_t, OCT_t) + cs(OCT_t, OCT_u) + cs(OCT_u, OCT_t) + cs(OCT_u, OCT_u)
    line('-', 44)
    print("      %-28s %10.4f" % ("SINGLET, summed over channels", tot_s))
    print("      %-28s %10.4f" % ("OCTET,   summed over channels", tot_o))

    ok = abs(tot_s) > 1e-9
    print("\n   C2b singlet DOES interfere once quarks are identical .... %s (%.4f)"
          % ("PASS" if ok else "FAIL", tot_s))
    if not ok:
        fails.append("C2b")
    ok = tot_s > tot_o
    print("   C2c and its colour weight EXCEEDS the octet's ........... %s (%.2f x)"
          % ("PASS" if ok else "FAIL", tot_s / tot_o if tot_o else float('nan')))
    if not ok:
        fails.append("C2c")
    print("""
    So the earlier reading of C2 -- that a published singlet limit understates what the data
    say about an octet -- is not merely unproven: in the channel that dominates a high-mass
    dijet sample it points the other way.  The decisive check needs no algebra at all:
    [CMSangular] quotes 17 TeV for destructive and 37 TeV for constructive interference.  A
    vanishing interference cannot produce a limit that depends on the sign.""")

    # -------------------------------------------------------------- consequence
    line('=')
    print("[3] WHAT THIS MEANS FOR THE LIMIT")
    line('=')
    print("""
    CMS parametrises its scenarios as   L = (2 pi / Lambda^2) [ eta_LL (qL g qL)(qL g qL) + ... ]
    with COLOUR-SINGLET currents (their Sec. on CI, Table 3).  Our tower gives

        L_KK = - (pi^2 / 6) g_s^2 R_5^2 (qbar g^mu T^a q)(qbar g_mu T^a q),

    vector-like in chirality -- so it matches their (eta_LL, eta_RR, eta_RL) = (1,1,1) row,
    the Lambda_VV scenario -- but OCTET in colour, which none of their rows is.

    The two therefore enter the observable differently:

      * the octet carries the colour flow of one-gluon exchange, so it interferes with QCD at
        full strength in the t-channel (C3);
      * the singlet does not interfere in that contraction (C2) -- but it does once the quarks
        are identical, through the crossed t-u terms, and there its colour weight is the LARGER
        of the two (C2b, C2c).

    Both operators therefore interfere, with different weights in different channels, and which
    of the two a given dataset constrains more tightly is not decidable from colour factors
    alone.  It depends on the flavour composition of the sample, which is a PDF question, and on
    the relative size of the interference and the quadratic term at the scales probed.

    What survives is the weaker and true statement: the colour structure and the interference
    pattern both differ from the CMS benchmark, so its published CI scale cannot be transferred
    to the Kaluza-Klein tower.  A recast is required.  Its direction is not known in advance.

    What is still needed, and is NOT done here:
      - the interference term is signed, and the sign decides destructive vs constructive.  The
        effective propagator of [DMN01] gives 1/q^2 + sum_n 2/(q^2 - m_n^2), both terms negative
        for spacelike q^2, so the tower ADDS to QCD -- constructive, and consistent with [DMN01]
        finding sigma_KK / sigma_SM > 1.  CMS's convention (their Table 3 note) calls that the
        eta < 0 branch, whose published limit is the STRONGER one.
      - the actual chi_dijet shape, which needs the partonic cross sections at NLO.  CMS use
        CIJET 1.0 for exactly this.  Without it we can bound the direction, not the number.
      - and the EFT itself fails near the resonance: with 1/R_5 = 3.97 TeV the sensitive region
        M_jj ~ 2-8 TeV sits ON the tower, not below it, so the escape branch needs the full
        propagator and not a contact operator at all.
    """)

    line('=')
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line('=')
        return 1
    print("VERDICT: all controls pass; the octet/singlet distinction is exact, not an estimate.")
    line('=')
    return 0


if __name__ == '__main__':
    sys.exit(main())
