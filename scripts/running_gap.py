#!/usr/bin/env python3
"""running_gap.py -- the four-dimensional running below 1/R5, and the anchor extrapolations.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS FILE EXISTS, AND IT IS A REPAIR
-----------------------------------------
Section 10 states that below 1/R5 the theory is the Standard Model, so the gap
alpha_2^-1 - alpha_3^-1 "is 19.2 at 2 TeV and 18.2 at the certified ceiling -- it barely moves
across the whole allowed range, because four-dimensional running has under five e-folds to work
with and closes the gap at 0.61 per e-fold".  And section 8's table quotes 8.2 and 20.8 TeV as
extrapolations of 3.97 and 10.03 by the largest measured anchor ratio.

None of those five numbers was computed by anything in the repository.  `check_numbers.py` found
them unbacked in the archive AND in a fresh run of every script, so this is the repair: the Data
availability section promises every displayed number regenerates from the ancillary scripts, and
until now that promise was false for these.

THE COMPUTATION, which is one loop and elementary
-------------------------------------------------
    alpha_i^-1(mu) = alpha_i^-1(M_Z) - (b_i / 2 pi) ln(mu / M_Z)

with the Standard-Model one-loop coefficients b = (41/10, -19/6, -7) for the three factors, and
the measured inputs at M_Z.  Nothing here is fitted; the inputs are PDG numbers and they are
printed so a reader can put their own in.

Run:  python running_gap.py
"""
import math
import sys

P = lambda *a: print(*a, flush=True)

MZ = 91.1876                 # GeV
ALPHA_EM_INV = 127.951       # alpha_em^-1(M_Z), MS-bar
SIN2W = 0.23121              # sin^2 theta_W(M_Z), MS-bar
ALPHA_S = 0.1179             # alpha_s(M_Z)

B2, B3 = -19.0 / 6.0, -7.0   # Standard-Model one-loop, SU(2)_L and SU(3)_C

A2_INV_MZ = ALPHA_EM_INV * SIN2W
A3_INV_MZ = 1.0 / ALPHA_S

# The largest measured anchor ratio, READ from the archive rather than typed.  The first draft
# typed 2.08 -- the rounded value that appears in the prose -- and the extrapolations came out at
# 20.9 and 8.3 against the paper's 20.8 and 8.2.  The paper is right and the typed constant was
# wrong: the ratio is 2.076 on row (3), and rounding it before multiplying moves the last digit.
def _anchor_ratio():
    import pathlib
    import re
    p = pathlib.Path(__file__).resolve().parent / "outputs" / "amin_closed_form.txt"
    rows = re.findall(r"^\(\d\)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+([\d.]+)",
                      p.read_text(encoding="utf-8"), re.M)
    if not rows:
        raise SystemExit("FATAL: cannot read the anchor ratios from %s" % p)
    return max(float(r) for r in rows)


RATIO = _anchor_ratio()


def inv(a_inv_mz, b, mu):
    return a_inv_mz - b / (2 * math.pi) * math.log(mu / MZ)


def gap(mu):
    return inv(A2_INV_MZ, B2, mu) - inv(A3_INV_MZ, B3, mu)


def main():
    P("=" * 96)
    P("1 -- THE GAP  alpha_2^-1 - alpha_3^-1  UNDER FOUR-DIMENSIONAL RUNNING")
    P("=" * 96)
    P("inputs at M_Z = %.4f GeV:  alpha_em^-1 = %.3f, sin^2 theta_W = %.5f, alpha_s = %.4f"
      % (MZ, ALPHA_EM_INV, SIN2W, ALPHA_S))
    P("   =>  alpha_2^-1(M_Z) = %.3f,  alpha_3^-1(M_Z) = %.3f,  gap = %.3f"
      % (A2_INV_MZ, A3_INV_MZ, A2_INV_MZ - A3_INV_MZ))
    P("one-loop Standard-Model coefficients:  b_2 = %s, b_3 = %s" % ("-19/6", "-7"))
    P("")
    P("%14s %12s %12s %10s %12s" % ("scale", "a2^-1", "a3^-1", "gap", "e-folds"))
    for label, mu in (("M_Z", MZ), ("1 TeV", 1000.0), ("2 TeV", 2000.0),
                      ("3.97 TeV", 3967.0), ("10.03 TeV", 10034.0)):
        P("%14s %12.3f %12.3f %10.2f %12.3f"
          % (label, inv(A2_INV_MZ, B2, mu), inv(A3_INV_MZ, B3, mu), gap(mu),
             math.log(mu / MZ)))
    P("")
    g2, gc = gap(2000.0), gap(10034.0)
    P("gap at 2 TeV                : %.2f" % g2)
    P("gap at the certified ceiling: %.2f" % gc)
    P("it moves by %.2f over %.2f e-folds -- %.2f per e-fold."
      % (g2 - gc, math.log(10034.0 / 2000.0), (g2 - gc) / math.log(10034.0 / 2000.0)))
    P("")
    P("  So four-dimensional running cannot close it: the whole burden falls on the")
    P("  Kaluza-Klein tower, and how far the gap has to be closed is %.1f, not %.1f."
      % (gc, A2_INV_MZ - A3_INV_MZ))

    P("")
    P("=" * 96)
    P("2 -- THE ANCHOR EXTRAPOLATIONS, which are extrapolations and not bounds")
    P("=" * 96)
    P("The right-hand column of section 8's table multiplies by %.2f, the largest ratio measured"
      % RATIO)
    P("over the five published rows.  Printed here so the numbers in that table regenerate.")
    P("")
    P("%-46s %12s %12s" % ("", "ours", "x %.2f" % RATIO))
    for label, v in (("content of their own size (<= 6 multiplets)", 2.18),
                     ("ceiling over all contents", 10.03),
                     ("ceiling for contents that can afford the escape", 3.97)):
        P("%-46s %10.2f T %10.2f T" % (label, v, v * RATIO))
    P("")
    P("  EXTRAPOLATIONS, NOT BOUNDS.  A ratio measured on five rows is a measurement of a")
    P("  discrepancy, not a weight, and a two-loop term is under no obligation to be a uniform")
    P("  rescaling of anything.  The paper marks them with a dagger for that reason.")

    P("")
    P("=" * 96)
    P("CONTROLS")
    P("=" * 96)
    ok1 = abs(gap(2000.0) - 19.2) < 0.1
    ok2 = abs(gap(10034.0) - 18.2) < 0.1
    P("  gap(2 TeV)   reproduces the 19.2 in the text : %s  (%.3f)" % (ok1, gap(2000.0)))
    P("  gap(ceiling) reproduces the 18.2 in the text : %s  (%.3f)" % (ok2, gap(10034.0)))
    # a control that can fail: with the WRONG sign of b_3 the gap would run the other way
    bad = inv(A2_INV_MZ, B2, 10034.0) - inv(A3_INV_MZ, -B3, 10034.0)
    P("  control -- flip the sign of b_3 and the gap becomes %.2f, i.e. it moves the other way"
      % bad)
    P("             and by %.1f instead of %.1f: the sign of the QCD beta function is doing the"
      % (abs(bad - (A2_INV_MZ - A3_INV_MZ)), abs(gc - (A2_INV_MZ - A3_INV_MZ))))
    P("             work, and the test can tell.")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main())
