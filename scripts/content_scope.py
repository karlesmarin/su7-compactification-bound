#!/usr/bin/env python3
"""content_scope.py -- "arbitrary bulk content" is arbitrary in the multiplicities, not the reps.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE WORD THAT WAS DOING WORK WITHOUT DECLARING IT.  This paper says "for arbitrary bulk content"
and means: arbitrary NON-NEGATIVE MULTIPLICITIES of the eight (representation, parity) pairs of
\\cite{KM25}'s Table 1.  It never says the second half.  And the omission is not cosmetic, because
the paper spent a whole section removing a cap on the SIZE of the content -- the 2.18 TeV that was
really an enumeration to six multiplets -- while leaving a cap on its TYPE unstated.

WHAT SETS THE CAP.  The Wilson line sees a state through its charge, and the charge of a weight is
q = (boxes)/2 in \\cite{KM25}'s normalisation, so the potential's integer charge c = 2|q| runs from
0 up to the NUMBER OF BOXES of the Young diagram.  Section 1 checks that against the term tables
themselves rather than against a memory of how they were derived:

    7 -> 1 box -> c <= 1     28 (sym^2) -> 2 -> c <= 2
    48 (adjoint) -> 2 -> c <= 2     84 -> 3 -> c <= 3

So c <= 3 in this paper because \\cite{KM25}'s largest representation has three boxes.  Nothing in
the model, the orbifold or the potential caps the number of boxes; SU(7) has representations with
four and more, and a bulk field in one of them would bring c = 4.

WHAT WOULD CHANGE IF ONE WERE ALLOWED, and it is more than the ceiling.  Section 2:

  - the ceiling.  A_4 collects m c^4, so a charge-four state is worth 256 where a charge-three
    state is worth 81.  The whole optimisation is over a different cone.
  - the DIMENSION.  With c in {1,2,3} the six functions span five, because the duplication formula
    relates g(+,1)+g(-,1) to g(+,2) and there is no c=4 to continue with.  Admit c=4 and the
    chain continues: g(+,2)+g(-,2) = g(+,4)/16 becomes a second relation, and eight functions span
    six.  Theorem (complete invariants) is a statement about THIS charge alphabet, and its five
    coordinates would become six.

Neither is a defect of the results; both are scope, and scope has to be printed.

Run:  python content_scope.py
"""
import json
import math
import pathlib
import sys
from fractions import Fraction as Fr

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
BOXES = {"7": 1, "28": 2, "48": 2, "84": 3}     # index-7 boxes; the adjoint is 7 x 7bar

# ================================================================= 1
P("=" * 100)
P("1 -- THE CHARGE CEILING IS THE BOX COUNT, CHECKED AGAINST THE TERM TABLES")
P("=" * 100)
P("  Part VI's charge assignment gives q = 7 -> -/+1/2, 28 -> -/+1, 84 -> -/+3/2, 48 -> 0, and")
P("  the potential's integer charge is c = 2|q|.  If that reading is right, the largest c in a")
P("  multiplet's term table must equal its number of boxes.  The term tables are the arbiter:")
P("")
P("  %-10s %8s %14s %10s %s" % ("multiplet", "boxes", "max c in table", "agree", "charges present"))
ok = True
for j, (r, e, p) in enumerate(REPS):
    tt = terms(r, e, p)
    cs = sorted({int(round(c)) for m, s, c in tt})
    good = max(cs) == BOXES[r]
    ok &= good
    P("  %-10s %8d %14d %10s %s" % (NAMES[j], BOXES[r], max(cs), good, cs))
P("")
P("  CONTROL -- max c equals the box count on all eight : %s" % ok)
assert ok, "the charge ceiling is not the box count -- the scope argument below does not follow"
P("")
P("  So c <= 3 holds because the largest representation in \\cite{KM25}'s Table 1 has three")
P("  boxes.  It is a property of THEIR LIST, not of the model: SU(7) has representations with")
P("  four boxes and more, and nothing in the orbifold, the potential or the ceiling forbids a")
P("  bulk field in one.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- WHAT A FOURTH BOX WOULD COST, AND IT IS NOT ONLY THE CEILING")
P("=" * 100)
P("  (a) THE CEILING.  A_4 collects m c^4 and 8D collects 8 m c^2, so the value of one periodic")
P("      state to the optimisation is c^4 against 8 c^2:")
P("")
P("  %8s %10s %10s %14s" % ("c", "c^4", "8 c^2", "A_4 per 8D"))
for c in (1, 2, 3, 4, 5):
    P("  %8d %10d %10d %14.3f" % (c, c ** 4, 8 * c ** 2, c ** 4 / (8.0 * c ** 2)))
P("")
P("      The ratio is c^2/8, so it grows without bound in the charge.  The ceiling is a maximum")
P("      of A_4 at fixed 8D, so admitting a higher charge admits a strictly better direction and")
P("      the certificate is over a different cone.  How much better is not computed here: it")
P("      needs the parity decomposition of a representation \\cite{KM25} do not use, and that is")
P("      their assignment to make, not ours.")
P("")
P("  (b) THE DIMENSION, which is the part that would be missed.  With charges {1,2,3} the six")
P("      functions g(s,c) span five, because Li_s(z)+Li_s(-z) = 2^{1-s} Li_s(z^2) relates")
P("      g(+,1)+g(-,1) to g(+,2) and the chain then STOPS -- continuing it would need c = 4.")
P("      Admit c = 4 and it does not stop.  Counting, with the relation available whenever 2c is")
P("      also a charge:")
P("")
KEYS_ALL = [(s, c) for s in (1, -1) for c in (1, 2, 3, 4, 5, 6)]
ys = np.linspace(0.0, 1.0, 2001)
P("  %-22s %10s %10s %12s %10s" % ("charge alphabet", "functions", "relations", "dimension",
                                   "this paper"))
for cmax in (3, 4, 5, 6):
    ch = list(range(1, cmax + 1))
    keys = [(s, c) for s in (1, -1) for c in ch]
    G = np.array([basis(ys, s, c) for s, c in keys])
    sv = np.linalg.svd(G, compute_uv=False)
    rk = int(np.sum(sv > sv[0] * 1e-10))
    P("  %-22s %10d %10d %12d %10s"
      % ("c in {1..%d}" % cmax, len(keys), len(keys) - rk, rk, "yes" if cmax == 3 else ""))
P("")
P("      Each doubling that lands back inside the alphabet costs one dimension, and the count is")
P("      the number of charges c with 2c also a charge.  For {1,2,3} that is one, c = 1, hence")
P("      six functions and five dimensions.  For {1,2,3,4} it is two, c = 1 and c = 2, hence")
P("      eight functions and six dimensions.")
P("")
P("      So the complete-invariants theorem is a statement about THIS charge alphabet.  Its five")
P("      coordinates are five because the alphabet stops at three, and a fourth box would make")
P("      them six.  That is scope, not error -- but unstated it would read as a general law.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- AND THE ASYMMETRY THAT MADE THIS WORTH CHECKING")
P("=" * 100)
P("  The paper corrected one cap on the content and left another standing:")
P("")
P("     the cap on SIZE  : 'at most six multiplets' -- CORRECTED, and by a factor 2.5 in the")
P("                        answer.  Nothing in the model caps how MANY multiplets there are.")
P("     the cap on TYPE  : 'the eight (rep, parity) pairs of their Table 1' -- NOT STATED.")
P("                        Nothing in the model caps WHICH representations there are either.")
P("")
P("  Both are properties of \\cite{KM25}'s table rather than of the theory, and the first was")
P("  worth a section.  The second is worth a sentence, and it now has one.")

out = dict(boxes=BOXES, charge_cap_is_box_count=bool(ok),
           dims={str(c): int(np.sum(np.linalg.svd(
               np.array([basis(ys, s, cc) for s in (1, -1) for cc in range(1, c + 1)]),
               compute_uv=False) > 1e-10 * np.linalg.svd(
               np.array([basis(ys, s, cc) for s in (1, -1) for cc in range(1, c + 1)]),
               compute_uv=False)[0])) for c in (3, 4, 5, 6)})
(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "content_scope.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
P("")
P("archived: outputs/content_scope.json")
