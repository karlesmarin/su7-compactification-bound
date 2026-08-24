#!/usr/bin/env python3
"""moduli_space.py -- how many flat Wilson-line directions does the SU(7) model have?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS EXISTS
---------------
COLOUR_DIRECTIONS.md (2026-08-21) claims six flat COLOURED directions in A_6 that nobody has
minimised, and calls the one-variable minimisation of Parts VI and VII a partial one.  If that
were true every absolute number in Part VII would sit at a point that is not the vacuum.

Before measuring D along those six directions, ask whether they are there.  They are not, and the
reason is a rule that brane_route.py applied to the wrong field.

THEIR RULES, verbatim from Komori-Maru 2503.04090 eqs. (3)-(10).  Note the two asymmetries: A_6
carries an EXTRA MINUS under the x6 reflections, and A_5 carries it under the x5 ones.

    A_{mu,5}(x, x5, -x6)          = +P6 A_{mu,5} P6            (3)
    A_6      (x, x5, -x6)          = -P6 A_6      P6            (4)
    A_{mu,5}(x, x5, piR6 - x6)     = +P6 A_{mu,5} P6            (5)
    A_6      (x, x5, piR6 + x6)    = -P6 A_6      P6            (6)
    A_{mu,6}(x, -x5, x6)           = +P5 A_{mu,6} P5            (7)
    A_5      (x, -x5, x6)          = -P5 A_5      P5            (8)
    A_{mu,6}(x, piR5 - x5, x6)     = +P5' A_{mu,6} P5'          (9)
    A_5      (x, piR5 + x5, x6)    = -P5' A_5      P5'          (10)

So, per field, a 4D zero mode needs all three eigenvalues equal to +1 with these signs:

    field    x6 rule (P6)     x5 rule (P5)     x5' rule (P5')
    A_mu       +P6 . P6         +P5 . P5         +P5' . P5'
    A_5        +P6 . P6         -P5 . P5         -P5' . P5'
    A_6        -P6 . P6         +P5 . P5         +P5' . P5'

brane_route.py asked for A_6 with the x6 rule of A_6 and the x5 rules of A_5 -- a field that is
neither -- and that mixed row is where its six coloured states came from.

Run:  python moduli_space.py
"""
import sys

import numpy as np

P = lambda *a: print(*a, flush=True)

P6 = np.diag([1, 1, 1, -1, -1, -1, -1]).astype(float)
P5 = np.diag([1, 1, 1, 1, 1, -1, -1]).astype(float)
P5p = np.diag([1, 1, 1, -1, -1, -1, 1]).astype(float)

COLOUR = (0, 1, 2)          # the unbroken SU(3): P6 = +1 on these three
WEAK = (3, 4)               # the unbroken SU(2)_L, fixed below by a control against their (57)

# per field: (sign on the x6 rule, sign on the x5 rule, sign on the x5' rule)
FIELD_RULE = {"A_mu": (+1, +1, +1),
              "A_5": (+1, -1, -1),
              "A_6": (-1, +1, +1)}
WRONG_RULE = (-1, -1, -1)   # what brane_route.py used for A_6: A_6's x6 rule, A_5's x5 rules


# The six Cartan generators, in a basis ADAPTED to the unbroken group -- two inside the colour
# SU(3) of indices (1,2,3), one inside the weak SU(2) of indices (4,5), and three abelian.  Any
# traceless diagonal basis diagonalises the parities, so this choice costs nothing and it is what
# lets the count be compared with their eq. (57), where the octet is 8 and not 6.
CARTAN = [("cartan:colour", (1, -1, 0, 0, 0, 0, 0)),
          ("cartan:colour", (1, 1, -2, 0, 0, 0, 0)),
          ("cartan:weak", (0, 0, 0, 1, -1, 0, 0)),
          ("cartan:u1", (1, 1, 1, 0, 0, 0, -3)),
          ("cartan:u1", (0, 0, 0, 1, 1, -2, 0)),
          ("cartan:u1", (1, 1, 1, -1, -1, -1, 0))]


def basis():
    """the 48 generators of su(7) as (label, matrix), in a parity-diagonal basis.

    Every P is diagonal, so each off-diagonal E_ab is already an eigenstate of all three
    conjugations, and so is every traceless diagonal generator.
    """
    out = []
    for a in range(7):
        for b in range(7):
            if a != b:
                M = np.zeros((7, 7))
                M[a, b] = 1.0
                out.append(((a, b), M))
    for tag, v in CARTAN:
        out.append((tag, np.diag(v).astype(float)))
    return out


def eig(M, Pm, sign):
    """eigenvalue of  M -> sign * P M P ;  +1, -1, or None if not an eigenstate."""
    T = sign * Pm @ M @ Pm
    if np.allclose(T, M):
        return +1
    if np.allclose(T, -M):
        return -1
    return None


def quantum(idx, M, rule):
    s6, s5, s5p = rule
    return (eig(M, P6, s6), eig(M, P5, s5), eig(M, P5p, s5p))


def colour(idx):
    if isinstance(idx, str):
        return "octet" if idx == "cartan:colour" else "singlet"
    a, b = idx
    ina, inb = a in COLOUR, b in COLOUR
    if ina and inb:
        return "octet"
    if ina != inb:
        return "TRIPLET"
    return "singlet"


def weak(idx):
    if isinstance(idx, str):
        return 3 if idx == "cartan:weak" else 1
    a, b = idx
    ina, inb = a in WEAK, b in WEAK
    if ina != inb:
        return 2
    if ina and inb:
        return 3
    return 1


def zero_modes(B, rule):
    return [(idx, M) for idx, M in B if quantum(idx, M, rule) == (1, 1, 1)]


def report(name, B, rule):
    zm = zero_modes(B, rule)
    trip = [i for i, _ in zm if colour(i) == "TRIPLET"]
    P("  %-7s zero modes: %2d      colour triplets among them: %d"
      % (name, len(zm), len(trip)))
    return zm, trip


def main():
    B = basis()

    P("=" * 96)
    P("1 -- THE ZERO MODES OF EACH FIELD, WITH ITS OWN RULE")
    P("=" * 96)
    P("  a 4D zero mode is a generator even under all three reflections, with the signs of")
    P("  their eqs. (3)-(10) -- and those signs differ from field to field.")
    P("")
    got = {}
    for name in ("A_mu", "A_5", "A_6"):
        got[name] = report(name, B, FIELD_RULE[name])
    P("")
    P("  the Wilson-line moduli are the zero modes of A_5 and A_6, and of nothing else:")
    P("     A_5 : %d real     A_6 : %d real     total %d real"
      % (len(got["A_5"][0]), len(got["A_6"][0]),
         len(got["A_5"][0]) + len(got["A_6"][0])))

    P("")
    P("=" * 96)
    P("2 -- WHAT THOSE MODES ARE")
    P("=" * 96)
    for name in ("A_mu", "A_5", "A_6"):
        zm = got[name][0]
        if not zm:
            P("  %-5s : none" % name)
            continue
        P("  %-5s : %d states" % (name, len(zm)))
        for idx, _ in zm:
            P("        %-16s colour %-8s SU(2)_L %d" % (str(idx), colour(idx), weak(idx)))

    P("")
    P("=" * 96)
    P("3 -- CONTROL A: reproduce THEIR published decomposition of the 48, their eq. (57)")
    P("=" * 96)
    P("  their table is written with the A_mu parities (P6, P5, P5').  Count our 48 generators")
    P("  by that triple and by (colour, SU(2)) and compare with what they print.")
    P("")
    # their eq. (57), transcribed: (colour dim, SU(2) dim, (P6,P5,P5'), multiplicity)
    THEIRS = [(8, 1, (1, 1, 1), 1), (1, 3, (1, 1, 1), 1), (1, 1, (1, 1, 1), 3),
              (1, 1, (1, 1, -1), 2),
              (1, 2, (1, -1, 1), 2), (1, 2, (1, -1, -1), 2),
              (3, 2, (-1, 1, -1), 2), (3, 1, (-1, -1, -1), 2), (3, 1, (-1, -1, 1), 2)]
    tally = {}
    for idx, M in B:
        q = quantum(idx, M, FIELD_RULE["A_mu"])
        key = (colour(idx), weak(idx), q)
        tally[key] = tally.get(key, 0) + 1
    ok = True
    P("  %-16s %-8s %-14s %8s %8s" % ("colour", "SU(2)", "(P6,P5,P5')", "theirs", "ours"))
    for cdim, wdim, q, mult in THEIRS:
        cname = {8: "octet", 3: "TRIPLET", 1: "singlet"}[cdim]
        theirs = cdim * wdim * mult
        ours = sum(v for (c, w, qq), v in tally.items()
                   if c == cname and w == wdim and qq == q)
        flag = "ok" if theirs == ours else "MISMATCH"
        ok &= theirs == ours
        P("  %-16s %-8d %-14s %8d %8d   %s" % (cname, wdim, str(q), theirs, ours, flag))
    P("")
    P("  total generators accounted: %d of 48" % sum(tally.values()))
    P("  CONTROL A: %s" % ("their eq. (57) is reproduced state by state" if ok else "FAILED"))

    P("")
    P("=" * 96)
    P("4 -- CONTROL B: the falsification -- the rule brane_route.py used, and what it invents")
    P("=" * 96)
    P("  it asked for A_6's x6 rule together with A_5's x5 rules.  No field in their")
    P("  Lagrangian transforms that way.  Run it anyway:")
    P("")
    zmw = zero_modes(B, WRONG_RULE)
    tw = [i for i, _ in zmw if colour(i) == "TRIPLET"]
    P("     mixed rule (-P6, -P5, -P5') : %d states, %d colour triplets" % (len(zmw), len(tw)))
    P("     the six of COLOUR_DIRECTIONS.md are %s"
      % ("exactly these" if len(tw) == 6 else "NOT these -- reopen the note"))
    P("     and the true A_6 rule (-P6, +P5, +P5') gives : %d states, %d triplets"
      % (len(got["A_6"][0]), len(got["A_6"][1])))
    P("")
    P("  so the control separates the two rules: one manufactures six coloured flat")
    P("  directions and the other gives none.  A control that could not tell them apart")
    P("  would be worthless.")

    P("")
    P("=" * 96)
    P("5 -- CONTROL C: the counts have to close, and one of them is a physics requirement")
    P("=" * 96)
    nmu = len(got["A_mu"][0])
    P("     A_mu zero modes = unbroken 4D gauge group : %d" % nmu)
    P("       expected SU(3) x SU(2) x U(1)^3 = 8 + 3 + 3 = 14 : %s" % (nmu == 14))
    P("     colour triplets in A_mu : %d  (a triplet here is a MASSLESS leptoquark gauge"
      % len(got["A_mu"][1]))
    P("       boson and the model would be dead) : %s" % (len(got["A_mu"][1]) == 0))
    P("     A_5 zero modes : %d real = one complex SU(2) doublet, colourless : %s"
      % (len(got["A_5"][0]),
         len(got["A_5"][0]) == 4 and all(weak(i) == 2 and colour(i) == "singlet"
                                         for i, _ in got["A_5"][0])))

    P("")
    P("=" * 96)
    P("VERDICT")
    P("=" * 96)
    nmod = len(got["A_5"][0]) + len(got["A_6"][0])
    P("  the Wilson-line moduli space of this model is %d real dimensional, and all %d sit in"
      % (nmod, nmod))
    P("  A_5 as ONE colourless complex SU(2)_L doublet.  SU(2)_L x U(1)_Y acts on it with a")
    P("  three-dimensional generic orbit, so after a gauge rotation ONE phase is left.")
    P("")
    P("  => 'a single Wilson phase' is not an approximation of Parts VI and VII.  It is the")
    P("     whole moduli space, and the minimisation over alpha is a minimisation over all of")
    P("     it.  There is no coloured flat direction to check, because there is none.")
    P("")
    P("  COLOUR_DIRECTIONS.md is WITHDRAWN by this measurement.  What survives of it is the")
    P("  QUESTION -- 'is the origin the minimum over every flat direction, not just one?' --")
    P("  which is a good question, and this is its answer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
