#!/usr/bin/env python3
"""rational_holonomies -- is the third minimum at y = 0.6809 a deformation of y = 2/3?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  semi_infinite.py found that the true-vacuum programme, run with no screen, generates a second
  cut on the rung 8D = 3, at y = 0.680868, forced by one content.  The paper records the number
  and leaves open why it sits there.  This file asks whether it sits there because of an exact
  rational holonomy nearby.

  THE BASIS.  Write P_q(y) = Re Li_5(e^{i q pi y}).  The antiperiodic functions are not new
  directions: Li_5(z) + Li_5(-z) = 2^-4 Li_5(z^2) gives

      Re Li_5(-e^{i c pi y})  =  (1/16) P_{2c}(y)  -  P_c(y),

  so with the charge alphabet c <= 3 every term of the potential lives in the span of
  P_1, P_2, P_3, P_4, P_6.  Five functions -- which is Theorem 3 seen from the phases rather
  than from the moments, and is a check on both.

  WHAT IS EXACT AT y = 2/3.  There q pi y = 2 pi q / 3, so the q = 3 and q = 6 modes sit at
  2 pi and 4 pi: their phases are trivial and their derivatives vanish.  For the rest,
  sum_n cos(2 pi n / 3) / n^5 = -(40/81) zeta(5) exactly, and S_4(4 pi/3) = -S_4(2 pi/3).  Hence

      F'(2/3)      = -pi S_4(2pi/3) [ C_1 - 2 C_2 + 4 C_4 ] ,
      F(2/3)-F(0)  = -(121/81) zeta(5) [ C_1 + C_2 + C_4 ] .

  Both are EXACT linear functionals of the content, in the same sense as F(1) - F(0) = (31/16)
  zeta(5) W is.  If the second one already violates the depth allowance on the offending
  content, the numerical cut at 0.680868 can be replaced by an exact one at 2/3.

Run:  python rational_holonomies.py > outputs/rational_holonomies.txt
"""

import cmath
import fractions
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

Z5 = 1.0369277551433699263
QS = (1, 2, 3, 4, 6)


def coeffs(tt):
    """(m, s, c) term table -> the C_q in F = sum_q C_q P_q."""
    C = {q: 0.0 for q in QS}
    for (m, s, c) in tt:
        c = int(c)
        if s > 0:
            C[c] = C.get(c, 0.0) + m
        else:
            C[c] = C.get(c, 0.0) - m
            C[2 * c] = C.get(2 * c, 0.0) + m / 16.0
    return C


def Fs(content, y):
    """F at a scalar alpha.  amin_closed_form's basis() is vectorised and hands back an array
    even for a scalar argument, which float() will not take."""
    import numpy as _np
    return float(_np.asarray(F(content, _np.asarray([float(y)]))).ravel()[0])


def P(q, y, nmax=200000):
    """Re Li_5(e^{i q pi y}) by its Fourier series."""
    th = q * math.pi * y
    return sum(math.cos(n * th) / n ** 5 for n in range(1, nmax + 1))


def S4(th, nmax=200000):
    return sum(math.sin(n * th) / n ** 4 for n in range(1, nmax + 1))


def line(ch="-", n=94):
    print(ch * n)


def main():
    fails = []
    line("=")
    print("THE THIRD MINIMUM AT y = 0.680868, AND THE EXACT HOLONOMY AT y = 2/3")
    line("=")

    # the content semi_infinite.py names, retyped from its archived output
    CONTENT = [("7", 1, 1, 6), ("7", 1, -1, 1), ("28", 1, -1, 5),
               ("48", 1, -1, 1), ("84", 1, 1, 2)]
    tt = table(CONTENT)
    C = coeffs(tt)

    print("\n[0] THE CONTENT, AND ITS COORDINATES IN THE PHASE BASIS")
    print("    6x7(+,+) + 7(+,-) + 5x28(+,-) + 48(+,-) + 2x84(+,+)   [semi_infinite.py]")
    mo = moments(CONTENT)
    print("    8D = %.0f   A_4 = %.0f" % (8 * mo["D"], mo["A4"]))
    print("\n    %-6s %14s" % ("q", "C_q"))
    for q in QS:
        print("    %-6d %14.6f" % (q, C[q]))

    # ---- control: the phase basis must reproduce the potential we already have -------------
    print("\n[1] CONTROL: DOES THE PHASE BASIS REPRODUCE F ?")
    worst = 0.0
    for y in (0.05, 0.2, 0.5, 0.680868, 0.9):
        direct = Fs(CONTENT, y)
        viaP = sum(C[q] * P(q, y) for q in QS)
        worst = max(worst, abs(direct - viaP))
        print("    y = %-10.6f  F direct = %14.6f   via P_q = %14.6f   diff %.2e"
              % (y, direct, viaP, abs(direct - viaP)))
    ok = worst < 1e-6
    print("\n   C1  the two representations agree ................. %s (worst %.1e)"
          % ("PASS" if ok else "FAIL", worst))
    if not ok:
        fails.append("C1")
    print("       so the Li_5(-z) identity is right and the basis is the same five functions")
    print("       Theorem 3 counts from the moments.")

    # ---- the exact values at 2/3 ------------------------------------------------------------
    print("\n[2] THE EXACT VALUES AT y = 2/3")
    for q in QS:
        want = Z5 if q in (3, 6) else -40.0 / 81.0 * Z5
        got = P(q, 2.0 / 3.0)
        tag = "zeta(5)" if q in (3, 6) else "-(40/81) zeta(5)"
        print("    P_%d(2/3) = %14.9f   expected %-18s = %14.9f" % (q, got, tag, want))
        if abs(got - want) > 1e-7:
            fails.append("C2")
    print("\n   C2  every P_q(2/3) is the predicted rational multiple of zeta(5)  %s"
          % ("PASS" if "C2" not in fails else "FAIL"))

    # ---- the two linear functionals ---------------------------------------------------------
    print("\n[3] THE TWO EXACT FUNCTIONALS AT y = 2/3")
    grad = C[1] - 2 * C[2] + 4 * C[4]
    depth = -(121.0 / 81.0) * Z5 * (C[1] + C[2] + C[4])
    print("\n    stationarity :  C_1 - 2 C_2 + 4 C_4  = %+.6f    (zero <=> 2/3 is an extremum)" % grad)
    print("    depth        :  F(2/3) - F(0)        = %+.6f" % depth)

    # every C_q here is a dyadic rational, so the depth is an EXACT rational multiple of zeta(5)
    # and deserves to be printed as one rather than as a float.
    Cf = {q: fractions.Fraction(C[q]).limit_denominator(10 ** 6) for q in QS}
    rat = -fractions.Fraction(121, 81) * (Cf[1] + Cf[2] + Cf[4])
    print("    and exactly   :  F(2/3) - F(0)        = %s zeta(5) = %+.9f"
          % (rat, float(rat) * Z5))
    ok = abs(float(rat) * Z5 - depth) < 1e-12
    print("\n   C3a the exact rational agrees with the float ...... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3a")

    fp = -math.pi * S4(2 * math.pi / 3.0) * grad
    fp_num = (Fs(CONTENT, 2.0 / 3.0 + 1e-5) - Fs(CONTENT, 2.0 / 3.0 - 1e-5)) / 2e-5
    print("\n    F'(2/3) from the formula   : %+.6f" % fp)
    print("    F'(2/3) by finite difference: %+.6f" % fp_num)
    ok = abs(fp - fp_num) < 1e-3 * max(1.0, abs(fp))
    print("\n   C3  the closed form for F'(2/3) is right ........... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C3")

    d_num = Fs(CONTENT, 2.0 / 3.0) - Fs(CONTENT, 0.0)
    ok = abs(depth - d_num) < 1e-6
    print("   C4  and so is the one for F(2/3) - F(0) ........... %s (%.6f vs %.6f)"
          % ("PASS" if ok else "FAIL", depth, d_num))
    if not ok:
        fails.append("C4")

    # ---- does the exact point already do the job? -------------------------------------------
    print("\n[4] WOULD A CUT AT 2/3 HAVE DONE THE WORK OF THE CUT AT 0.680868 ?")
    d_true = Fs(CONTENT, 0.680868) - Fs(CONTENT, 0.0)
    print("\n    F(0.680868) - F(0) = %+.6f   <- the numerical cut semi_infinite.py found" % d_true)
    print("    F(2/3)      - F(0) = %+.6f   <- the exact rational point" % depth)
    print("    the two differ by    %.4f%%" % (100 * abs(depth - d_true) / abs(d_true)))
    # Theta is PER VERTEX, not a constant of the paper.  This cut was forced at (A_4, 8D) =
    # (147, 3), not at the deciding vertex (104, 1) -- and the two differ by a factor of eight,
    # so the allowance has to be recomputed here rather than copied across rungs.  The closed
    # form is semi_infinite.py's depth_bound(); C5a checks the copy against its archived table
    # before C5 leans on it.
    MW, G4 = 80.4, 0.63
    KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
    MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
    Z3 = 1.2020569031595942854
    x_of = lambda t, k: math.sqrt(12 * Z3 * (k / 8.0) / (6 * MU(127.0) + t))
    theta_of = lambda t, k: 10.0 * (Z3 * (k / 8.0) * x_of(t, k) ** 2 / 8.0
                                    + MU(125.0) * x_of(t, k) ** 4 / 16.0)

    print("\n    Theta is per vertex.  Against semi_infinite.py's archived table on 8D = 1:")
    ok = True
    for t, want in ((215, 7.979e-4), (212, 8.026e-4), (104, 1.016e-3), (92, 1.046e-3)):
        got = theta_of(t, 1)
        good = abs(got - want) < 5e-7
        ok = ok and good
        print("       A_4 = %-5d Theta = %.4e   archived %.4e   %s"
              % (t, got, want, "ok" if good else "MISMATCH"))
    print("\n   C5a the copied closed form reproduces that table ... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C5a")

    THETA = theta_of(147, 3)
    print("\n    and at the vertex this cut was forced at, (A_4, 8D) = (147, 3):")
    print("       Theta = %.4e   -- eight times the 1.016e-03 of the deciding vertex" % THETA)
    ok = depth < -THETA
    print("       the depth at 2/3 is %.0f times it" % (abs(depth) / THETA))
    print("\n   C5  the EXACT point at 2/3 already violates -Theta . %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C5")
    print("       so the cut does not need the numerical minimum: the rational holonomy suffices.")

    # ---- and the same game at 1/2 -----------------------------------------------------------
    print("\n[5] AND THE SAME GAME AT y = 1/2")
    print("""
    At y = 1/2 the phases are q pi / 2, so the even charges q = 2, 4, 6 land on multiples of pi
    and drop out, while S_4(3pi/2) = -S_4(pi/2) = -beta(4), the Dirichlet beta at four.  Hence""")
    B4 = S4(math.pi / 2.0)
    grad2 = C[1] - 3 * C[3]
    fp2 = (Fs(CONTENT, 0.5 + 1e-5) - Fs(CONTENT, 0.5 - 1e-5)) / 2e-5
    pred2 = -math.pi * B4 * grad2
    print("        F'(1/2) = -pi beta(4) [ C_1 - 3 C_3 ],     beta(4) = %.9f" % B4)
    print("\n    C_1 - 3 C_3 = %+.6f" % grad2)
    print("    F'(1/2) from the formula    : %+.6f" % pred2)
    print("    F'(1/2) by finite difference: %+.6f" % fp2)
    ok = abs(pred2 - fp2) < 1e-3 * max(1.0, abs(pred2))
    print("\n   C6  the closed form for F'(1/2) is right .......... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C6")
    print("""
    So there are at least two exact rational holonomies with exact linear stationarity criteria:

        y = 2/3   is an extremum  <=>  C_1 - 2 C_2 + 4 C_4 = 0
        y = 1/2   is an extremum  <=>  C_1 - 3 C_3         = 0

    and y = 0, 1 are the two symmetric points the ceiling already uses.  With the charge alphabet
    c <= 3 the candidates are the rationals whose denominators divide a charge, which is the
    finite set {0, 1/3, 1/2, 2/3, 1}.  Whether the semi-infinite programme's cuts always sit near
    one of them is NOT shown here -- one content is one data point.""")

    line("=")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(sorted(set(fails)))))
        line("=")
        return 1
    print("VERDICT: the numerical cut at 0.680868 is a deformation of the exact rational")
    print("         holonomy y = 2/3, and the exact point already carries the cut.")
    line("=")
    return 0


if __name__ == "__main__":
    sys.exit(main())
