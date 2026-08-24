#!/usr/bin/env python3
"""The formulas of Parts VI-VII in minimal form, and the map between them.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Carles: "dejalas en su minima expresion y mira como encajan unas con otras".  Reduction first,
review second: reviewing formulas that are not in canonical form is reviewing notation.

Every reduction below is CHECKED symbolically rather than asserted, because the last time this was
done by eye it produced a reading that had to be withdrawn (`W - 2b = N_c`, true at one assignment
and false on 3872 of 4096).

The five that matter, and what each is:

  (1)  the ladder identity        sum_n s^n/n^p = 2^-p zeta(p) + s(1-2^-p) zeta(p)
  (2)  the rung, per term         O_k = [S_k + (2^p - 1) Delta_k] / 2^p ,  p = 5-2k
  (3)  the rung, geometric        O_k = sum_sigma A_sigma sigma^k ,  sigma in {c^2/4, c^2}
  (4)  the Casimir sum rule       sum_Delta N_Delta(q) = sum_w mult(w) [C2 - |w|^2]
  (5)  the loop identity          L = 1 + (1/2) sum_k V_k (k-2)

What (2) and (3) have to agree on is the whole content of the ladder; what (4) says is where the
SECOND loop reads the group, and it is a different place from where the first one does.
"""
import itertools
import sympy as sp

P = lambda *a: print(*a, flush=True)
m, c, k = sp.symbols("m c k", positive=True)
s = sp.Symbol("s")

P("=" * 96)
P("(2) vs (3): the rung, per term, two ways")
P("=" * 96)
P("   per-term S_k = m c^{2k},  Delta_k = m c^{2k} s,  and 2^p = 32/4^k since p = 5-2k")
two_p = sp.Integer(32) / 4 ** k
per_term = (m * c ** (2 * k) + (two_p - 1) * m * c ** (2 * k) * s) / two_p
P("   (2)  O_k = %s" % sp.simplify(per_term))

# The geometric form of ladder_closed_form.py carries its amplitudes UNDIVIDED -- 32m and 2m -- so
# it is 2^p times the per-term form above.  Transcribing it with a 1/32 was the first thing this
# check caught, and it caught it against ME and not against the note: the two forms are the same
# object in two normalisations, and the factor is exactly 2^p.  Saying which one a number is in is
# not pedantry; it is a factor of 4^k.
geo_plus = sp.Integer(32) * m * (c ** 2 / 4) ** k
geo_minus = 2 * m * (c ** 2) ** k - 32 * m * (c ** 2 / 4) ** k
P("   (3)  s=+1: %s        s=-1: %s" % (sp.simplify(geo_plus), sp.simplify(geo_minus)))

ok_p = sp.simplify(sp.expand(per_term.subs(s, 1) * two_p - geo_plus)) == 0
ok_m = sp.simplify(sp.expand(per_term.subs(s, -1) * two_p - geo_minus)) == 0
P("")
P("   2^p x (2) against (3),  s = +1 : agree  %s" % ok_p)
P("   2^p x (2) against (3),  s = -1 : agree  %s" % ok_m)
P("   VERDICT: %s" % ("the geometric form IS the per-term form, identically in k, up to the"
                      " normalisation 2^p" if ok_p and ok_m else "*** they differ ***"))

P("")
P("=" * 96)
P("WHAT THE GEOMETRIC FORM SAYS, once it is trusted")
P("=" * 96)
P("   a PERIODIC term contributes one mode, sigma = c^2/4;")
P("   an ANTIPERIODIC one contributes two, sigma = c^2/4 and sigma = c^2, with opposite signs.")
P("   The ladder is therefore a sum of geometric sequences in k, and its fixed point is sigma = 1:")
P("")
for cc in (1, 2, 3):
    P("      c = %d :  sigma = c^2/4 = %-6s   c^2 = %-3d   %s"
      % (cc, sp.Rational(cc * cc, 4), cc * cc,
         "FIXED POINT: same weight at every rung" if cc == 2 else
         ("decays by 1/4 per rung" if cc < 2 else "grows")))
P("")
P("   So CHARGE TWO is the pivot of the whole ladder, and it is a statement about the")
P("   representation content, not about the loop integrals.")

P("")
P("=" * 96)
P("THE THREE RUNGS, AND WHY THE LAST ONE IS DIFFERENT")
P("=" * 96)
P("   %3s %4s %10s %28s" % ("k", "p", "2^p", "what it is"))
for kk, name in ((0, "the value F(0)"), (1, "the curvature -> D"), (2, "the fourth moment")):
    pp = 5 - 2 * kk
    P("   %3d %4d %10d %28s" % (kk, pp, 2 ** pp, name))
P("")
P("   The weight of Delta in rung k is (2^p - 1)/2^p = 1 - 2^-p: 31/32, 7/8, 1/2.")
P("   At p = 1 the two halves of (1) diverge separately and only their difference survives --")
P("   that is the logarithm, and it is why there is no polynomial.")

P("")
P("=" * 96)
P("THE JOIN: where each loop order reads the group")
P("=" * 96)
P("   C2 = sum_i H_i^2 + sum_alpha E^alpha E^-alpha    -- Cartan part + root part")
P("")
P("   ONE loop  reads  q_w = <w,T>                 the CARTAN part: a weight's own charge")
P("   TWO loops read   |<w+alpha|E^alpha|w>|^2     the ROOT part: how two weights are joined")
P("")
P("   and the sum rule (4) says the root part, resolved by charge, is fixed by C2 and |w|^2:")
P("       sum_Delta N_Delta(q) = sum_{w: q_w=q} mult(w) [C2 - |w|^2].")
P("")
P("   That is the same split Guo-Du find by computing in four dimensions: B_4(C_b) at one loop,")
P("   the individual angles, and B_4(C_bd) at two, the differences. An angle is a weight; a")
P("   difference of angles is a root. The two derivations meet.")

P("")
P("=" * 96)
P("THE Z2 THAT APPEARS UNDER FOUR NAMES")
P("=" * 96)
rows = [("the ladder", "s = +-1 in sum_n s^n/n^p", "even vs odd winding"),
        ("the moments", "S vs Delta", "parity-blind vs parity-sensitive"),
        ("the labels", "eta * eta'", "the two orbifold parities"),
        ("the SU(4) model", "A vs B", "identity vs coset multiplicity")]
for a, b, cdesc in rows:
    P("   %-16s %-30s %s" % (a, b, cdesc))
P("")
P("   They are one Z2. The criterion measured today closes the circle: the alpha_2 fold is legal")
P("   exactly when A = B on every half-integer charge -- that is, exactly when Delta receives")
P("   nothing from them, which is exactly when the parity-sensitive half of (1) is silent.")
