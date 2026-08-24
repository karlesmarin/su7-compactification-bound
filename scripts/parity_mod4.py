#!/usr/bin/env python3
"""Is D != 0 an accident of their parity choice, or a rule?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Root of the question (socratic pass, 2026-08-08).  D = (4 A2 - 3 B2)/4, so marginal
EWSB (D = 0) means exactly  A2/B2 = 3/4 = eta(3)/zeta(3).  In Komori-Maru that is
impossible because the ANTIPERIODIC gauge coefficient of their eq. (68) is 7, an ODD
integer, so B2 is a half-odd-integer.  su7_gauge_from_group.py shows where the 7 comes
from: weight 2 per (P6=+1) pair and 1/2 per (P6=-1) pair, so

    W_anti = 2 n+  +  n- / 2 ,   odd  <=>  n- = 2 (mod 4).

n- is the number of |q|=1/2 pairs of the adjoint with P5P5' = -1 and P6 = -1.  This
script asks whether n- = 2 (mod 4) is forced or a coincidence, over every diagonal
+-1 parity assignment of SU(7) obeying their own precondition (P5P5' and P6 must take
equal values on the two Wilson-line indices 5 and 7, else the charge eigenbasis is not
a parity eigenbasis -- their unstated pi_5 = pi_7 condition).

Adjoint component (i,j): q = (w_i - w_j)/2 with w_k = delta_{k,5} - delta_{k,7};
P5P5' and P6 act as products of the diagonal entries.
"""
import itertools
from collections import Counter

P = lambda *a: print(*a, flush=True)
N = 7
I5, I7 = 4, 6                      # 0-based indices 5 and 7


def w(k):
    return (1 if k == I5 else 0) - (1 if k == I7 else 0)


def classify(pp, p6):
    """counts of |q|=1/2 PAIRS by (P5P5', P6). pp, p6 are tuples of +-1."""
    c = Counter()
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            q2 = w(i) - w(j)                     # = 2q
            if abs(q2) != 1:
                continue
            if q2 < 0:                            # count each +-pair once
                continue
            c[(pp[i] * pp[j], p6[i] * p6[j])] += 1
    return c


REF_PP = tuple(a * b for a, b in zip([1, 1, 1, 1, 1, -1, -1], [1, 1, 1, -1, -1, -1, 1]))
REF_P6 = (1, 1, 1, -1, -1, -1, -1)

P("=" * 84)
P("CONTROL -- their own assignment must give n- = 6 and n+ = 2, hence W_anti = 7")
P("=" * 84)
c = classify(REF_PP, REF_P6)
np_, nm = c[(-1, +1)], c[(-1, -1)]
P("  P5P5' = %s" % (REF_PP,))
P("  P6    = %s" % (REF_P6,))
P("  pairs with P5P5'=-1:  P6=+1 -> n+ = %d ,  P6=-1 -> n- = %d" % (np_, nm))
W = 2 * np_ + nm / 2
P("  W_anti = 2*%d + %d/2 = %s   %s" % (np_, nm, W, "ODD -> D != 0" if W % 2 == 1 else "EVEN"))
P("  matches their eq. (68) coefficient 7 : %s" % (W == 7))

P("")
P("=" * 84)
P("THE SWEEP -- every diagonal +-1 (P5P5', P6) obeying their pi_5 = pi_7 precondition")
P("=" * 84)
tally, tally_odd = Counter(), Counter()
tot = 0
for pp in itertools.product([1, -1], repeat=N):
    if pp[I5] != pp[I7]:                 # precondition, on P5P5'
        continue
    for p6 in itertools.product([1, -1], repeat=N):
        if p6[I5] != p6[I7]:             # and on P6
            continue
        c = classify(pp, p6)
        np_, nm = c[(-1, +1)], c[(-1, -1)]
        W = 2 * np_ + nm / 2
        tot += 1
        tally[nm % 4] += 1
        if float(W).is_integer() and int(W) % 2 == 1:
            tally_odd[nm % 4] += 1

P("  assignments obeying the precondition: %d" % tot)
P("")
P("  %-14s %10s %14s" % ("n- mod 4", "count", "of which D != 0"))
for r in sorted(tally):
    P("  %-14d %10d %14d" % (r, tally[r], tally_odd.get(r, 0)))
P("")
odd = sum(tally_odd.values())
P("  assignments with W_anti odd (so D can NEVER vanish): %d of %d  = %.1f %%" %
  (odd, tot, 100 * odd / tot))
P("  and every one of them has n- = 2 (mod 4): %s" %
  (all(r == 2 for r in tally_odd if tally_odd[r]),))


# ============================================================================
# SEGUNDA PASADA SOCRATICA: n- = 2m, y quien es m
# ============================================================================
def m_count(pp, p6):
    """indices outside the Wilson-line pair that are doubly flipped w.r.t. it."""
    eps, dlt = pp[I5], p6[I5]
    return sum(1 for k in range(N) if k not in (I5, I7)
               and pp[k] == -eps and p6[k] == -dlt)


P("")
P("=" * 84)
P("D -- the mod-4 rule collapses:  pairs (5,k) and (7,k) always share parities")
P("=" * 84)
bad = 0
for pp in itertools.product([1, -1], repeat=N):
    if pp[I5] != pp[I7]:
        continue
    for p6 in itertools.product([1, -1], repeat=N):
        if p6[I5] != p6[I7]:
            continue
        c = classify(pp, p6)
        nm = c[(-1, -1)]
        mm = m_count(pp, p6)
        W = 2 * c[(-1, +1)] + nm / 2
        if nm != 2 * mm or ((W % 2 == 1) != (mm % 2 == 1)):
            bad += 1
P("  over all 4096 assignments:  n- = 2m  AND  (D never vanishes  <=>  m odd)")
P("  counterexamples: %d   %s" % (bad, "PASS" if bad == 0 else "FAILED"))

P("")
P("  their own assignment: m = %d" % m_count(REF_PP, REF_P6))
eps, dlt = REF_PP[I5], REF_P6[I5]
who = [k + 1 for k in range(N) if k not in (I5, I7)
       and REF_PP[k] == -eps and REF_P6[k] == -dlt]
P("  and the indices that make it up are %s" % (who,))
P("  -- indices 1,2,3 are the SU(3)_C colour indices of their eq. (78),(79).")
P("")
P("  >> D CANNOT VANISH BECAUSE THE NUMBER OF COLOURS IS ODD.")
