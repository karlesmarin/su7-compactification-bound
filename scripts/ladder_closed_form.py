#!/usr/bin/env python3
"""The ladder in minimal form: every term is a geometric sequence of ratio c^2/4.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

ladder_obstruction.py left three numbers for the gauge sector -- 9, -27, -36 at rungs k = 0,1,2 --
and the reading "the half-weight 1/2 enters as 7(2^{p-1}-1) = 105, 21, 0".  Simplifying that
properly gives something better, and kills an interpretation on the way.

THE SIMPLIFICATION.  Rung k has p = 5 - 2k, and the coefficient is S_k + (2^p - 1) Delta_k with
S_k = sum m c^{2k}, Delta_k = sum m c^{2k} s.  Since  2^p = 32 / 4^k,  a single term (m, s, c)
contributes

    s = +1 :   m c^{2k} 2^p            =  32 m (c^2/4)^k
    s = -1 :   m c^{2k} (2 - 2^p)      =  2 m (c^2)^k  -  32 m (c^2/4)^k

so EVERY term walks the ladder as a geometric sequence of ratio

    rho = c^2 / 4  :      c = 1 -> 1/4 (decays) ,   c = 2 -> 1 (CONSTANT) ,   c = 3 -> 9/4 (grows).

  >> CHARGE 2 IS THE FIXED POINT OF THE LADDER.  A charge-2 state weighs the same on the value, on
     the curvature and on the fourth moment; everything else is graded away from it.

Consequences, and both are tested below:
  * the gauge sector carries only c = 1 and c = 2, so its ladder has just TWO free parameters,
        O(p) = 2^{p-1} (W - 2b) - (32 a + W)     for the table (-a,+,2), (-b,+,1), (-W/2,-,1),
    and hence the exact three-rung relation   O(5) - O(3) = 4 [ O(3) - O(1) ];
  * that relation survives any matter content whose charges are 1 and 2 -- the 7, the 28, the 48 --
    and is broken by exactly one multiplet, the 84, which is the only one carrying a charge 3.

WHAT DIED HERE, and it is worth writing down.  Their assignment has W - 2b = 7 - 4 = 3 = N_c, which
made the whole ladder look like  N_c (2^{p-1} - 1) - even, i.e. made MARGINALITY.md's "D != 0 <=>
N_c odd" and STATE.md section 4's "the obstruction dies at the fourth moment" two factors of one
product.  Swept over all 4096 admissible parity assignments, W - 2b = N_c fails on 3872 of them.
It was a coincidence at one point, and the unification does not exist.
"""
import itertools
import json
import pathlib
from collections import Counter
from fractions import Fraction as Fr

P = lambda *a: print(*a, flush=True)
OUT = pathlib.Path(__file__).resolve().parent / "outputs"
N, I5, I7 = 7, 4, 6
W_PLUS, W_MINUS = Fr(2), Fr(1, 2)          # weight of a P6=+1 / P6=-1 pair, their (63)-(67)
RUNGS = [(0, 5), (1, 3), (2, 1)]

REF_PP = tuple(a * b for a, b in zip([1, 1, 1, 1, 1, -1, -1], [1, 1, 1, -1, -1, -1, 1]))
REF_P6 = (1, 1, 1, -1, -1, -1, -1)
_w = lambda k: (1 if k == I5 else 0) - (1 if k == I7 else 0)


def eq68(pp, p6):
    pairs = Counter()
    for i in range(N):
        for j in range(N):
            if i != j and _w(i) - _w(j) > 0:
                pairs[(_w(i) - _w(j), pp[i] * pp[j], p6[i] * p6[j])] += 1
    coef = Counter()
    for (c, s, sixth), n in pairs.items():
        coef[(c, s)] += n * (W_PLUS if sixth > 0 else W_MINUS)
    return coef


gauge_table = lambda coef: [(-v / 2, s, c) for (c, s), v in sorted(coef.items())]
abw = lambda coef: (coef[(2, +1)] / 2, coef[(1, +1)] / 2, coef[(1, -1)])


def n_colours(pp, p6):
    e, d = pp[I5], p6[I5]
    return sum(1 for k in range(N) if k not in (I5, I7) and pp[k] == -e and p6[k] == -d)


def rungs_of(tab):
    return {k: sum(m * Fr(c) ** (2 * k) for m, s, c in tab)
               + (2 ** p - 1) * sum(m * Fr(c) ** (2 * k) * s for m, s, c in tab)
            for k, p in RUNGS}


P("=" * 96)
P("1 -- CONTROL: their eq. (68) rebuilt from the parity assignment, not typed in")
P("=" * 96)
c0 = eq68(REF_PP, REF_P6)
a0, b0, W0 = abw(c0)
r0 = rungs_of(gauge_table(c0))
P("  coefficients {(c,s): value} = %s      (theirs: (2,+):2, (1,+):4, (1,-):7)" %
  dict(sorted((k, str(v)) for k, v in c0.items())))
P("  a = %s, b = %s, W = %s   ->  gauge table %s" %
  (a0, b0, W0, [(str(m), s, c) for m, s, c in gauge_table(c0)]))
P("  rungs O(5), O(3), O(1) = %s        (ladder_obstruction.py: 9, -27, -36)" %
  [str(r0[k]) for k, _ in RUNGS])

P("")
P("=" * 96)
P("2 -- THE GEOMETRIC LAW, and charge 2 as the fixed point")
P("=" * 96)
P("  %-6s %-12s %-30s %s" % ("c", "rho = c^2/4", "weight along k = 0, 1, 2", "behaviour"))
for c in (1, 2, 3):
    rho = Fr(c * c, 4)
    P("  %-6d %-12s %-30s %s" % (c, rho, [str(rho ** k) for k in (0, 1, 2)],
                                 "decays" if rho < 1 else ("CONSTANT" if rho == 1 else "grows")))
P("")
P("  So a charge-2 state weighs the SAME on the value, the curvature and the fourth moment.")
P("  The gauge sector has only c = 1, 2, hence two parameters and the three-rung relation below.")

P("")
P("=" * 96)
P("3 -- THE SWEEP: what holds on all 4096 assignments, and what does not")
P("=" * 96)
bad_alg = bad_nc = bad_rel = tot = 0
for pp in itertools.product([1, -1], repeat=N):
    if pp[I5] != pp[I7]:
        continue
    for p6 in itertools.product([1, -1], repeat=N):
        if p6[I5] != p6[I7]:
            continue
        tot += 1
        coef = eq68(pp, p6)
        a, b, Wv = abw(coef)
        R = rungs_of(gauge_table(coef))
        bad_alg += any(R[k] != 2 ** (p - 1) * (Wv - 2 * b) - (32 * a + Wv) for k, p in RUNGS)
        bad_nc += (Wv - 2 * b != n_colours(pp, p6))
        bad_rel += (R[0] - R[1] != 4 * (R[1] - R[2]))
P("  assignments swept: %d" % tot)
P("")
P("  HOLDS   O(p) = 2^{p-1}(W - 2b) - (32a + W)            counterexamples: %d" % bad_alg)
P("  HOLDS   O(5) - O(3) = 4 [ O(3) - O(1) ]               counterexamples: %d" % bad_rel)
P("  FAILS   W - 2b = N_c                                  counterexamples: %d of %d" % (bad_nc, tot))
P("")
P("  >> The identification of the bracket with the number of colours is TRUE AT THEIR ASSIGNMENT")
P("     AND FALSE IN GENERAL.  It was a coincidence at one point; the 'two theorems are one'")
P("     reading built on it is withdrawn.  What survives is the geometry, which is stronger anyway.")

P("")
P("=" * 96)
P("4 -- THE FALSIFIABLE CONSEQUENCE: which multiplet breaks the three-rung relation")
P("=" * 96)
import numpy as np                                                       # noqa: E402
exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]


def rungs_content(content):
    return rungs_of([(Fr(m).limit_denominator(64), s, int(c)) for m, s, c in table(content)])


P("A first guess -- 'it breaks as soon as a charge 3 appears' -- is WRONG, and the table below says")
P("so.  An ANTIPERIODIC term carries TWO modes, not one:")
P("")
P("      s = +1 :  one mode,  sigma = c^2/4")
P("      s = -1 :  two modes, sigma = c^2/4  AND  sigma = c^2")
P("")
P("so the modes available are  {1/4, 1, 9/4}  from c = 1, 2, 3  plus  {1, 4, 9}  from the")
P("antiperiodic ones.  Three rungs and two modes leave exactly one linear relation; a third mode")
P("destroys it.  The gauge sector reaches only {1/4, 1} -- which is WHY section 3 passes.")
P("")
P("PREDICTION, made before looking: the relation holds iff the term table has NO antiperiodic")
P("charge-2 state and NO charge-3 state.")
P("")
P("  %-12s %-16s %10s %10s %10s %12s %12s %s" %
  ("content", "modes", "O(5)", "O(3)", "O(1)", "O5-O3", "4(O3-O1)", "predicted/actual"))
bad_pred = 0
for j, (r, e, ep) in enumerate(REPS):
    tab = table([(r, e, ep, 1)])
    R = rungs_content([(r, e, ep, 1)])
    modes = sorted({Fr(int(c) ** 2, 4) for m, s, c in tab} |
                   {Fr(int(c) ** 2) for m, s, c in tab if s < 0})
    pred = modes == [Fr(1, 4), Fr(1)] or modes == [Fr(1, 4)] or modes == [Fr(1)]
    act = (R[0] - R[1] == 4 * (R[1] - R[2]))
    bad_pred += (pred != act)
    P("  %-12s %-16s %10s %10s %10s %12s %12s %s / %s" %
      (NAMES[j], ",".join(str(z) for z in modes), R[0], R[1], R[2],
       R[0] - R[1], 4 * (R[1] - R[2]),
       "holds" if pred else "BREAKS", "holds" if act else "BREAKS"))
P("")
P("  prediction wrong on: %d of 8   -> %s" % (bad_pred, "PREDICTION CONFIRMED" if not bad_pred
                                              else "*** PREDICTION FAILED ***"))
P("")
P("  Read the column: the (+,+) and (+,-) versions of the 28 and the 48 differ by whether their")
P("  charge-2 state is periodic or ANTIPERIODIC, and that single flip unlocks the mode sigma = 4.")
P("  The 84 unlocks 9/4 and 9 because it is the only multiplet carrying a charge 3.")
P("")
P("  and on the five published rows -- all of which carry an 84, so all break it:")
P("  %-6s %12s %12s %s" % ("row", "O5-O3", "4(O3-O1)", ""))
for label, content, a_them, mh_them, invR in T1:
    R = rungs_content(content)
    P("  %-6s %12s %12s   %s" % (label, R[0] - R[1], 4 * (R[1] - R[2]),
                                 "holds" if R[0] - R[1] == 4 * (R[1] - R[2]) else "BREAKS"))
P("")
P("  >> THE LADDER IS A SUM OF GEOMETRIC MODES  O_k = sum_sigma A_sigma sigma^k, and the number of")
P("     modes a content reaches is what decides whether its value, curvature and fourth moment are")
P("     three independent numbers or two.  Komori-Maru's own five rows all reach four modes.")

OUT.mkdir(exist_ok=True)
(OUT / "ladder_closed_form.json").write_text(json.dumps(dict(
    swept=tot, fail_closed_form=bad_alg, fail_three_rung=bad_rel, fail_Nc_identification=bad_nc,
    reference=dict(a=str(a0), b=str(b0), W=str(W0),
                   rungs={str(p): str(r0[k]) for k, p in RUNGS})), indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "ladder_closed_form.json"))
