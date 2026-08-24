#!/usr/bin/env python3
"""holonomy_skeleton -- do the true-vacuum cuts land on the rational skeleton, or anywhere?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  rational_holonomies.py shows that ONE content -- the one whose cut semi_infinite.py placed at
  y = 0.680868 -- is a half-per-cent deformation of the exact holonomy y = 2/3.  One content is
  one data point.  This file asks the question that makes it a claim or kills it: over many
  contents, where does the global minimiser of F actually sit?

  With the charge alphabet c <= 3, the phases q pi y become trivial at the rationals whose
  denominator divides a charge: y in {0, 1/3, 1/2, 2/3, 1}.  That is the skeleton.

  THE CONTROL THAT MATTERS.  A distance to the nearest of three points is small by construction.
  So the measured mean is compared against a UNIFORM null on the same range.  Without that the
  test could not fail and would prove nothing.

  WHY IT IS FAST.  A first version called amin_closed_form.basis once per term per content, and
  basis builds a (points x 600) outer product and cosines it -- about 1.7 s a content, three
  hours for the lattice.  All of that is waste: by the identity Li_5(z) + Li_5(-z) = 2^-4
  Li_5(z^2), every term lies in the span of P_q(y) = Re Li_5(e^{i q pi y}) for q in 1,2,3,4,6,
  so the five functions are computed ONCE on the grid and each content is a linear combination
  of them.  Control C0 checks the fast path against the slow one before using it.

  It checkpoints as it goes, so a run that is interrupted leaves usable partial data.

Run:  python holonomy_skeleton.py > outputs/holonomy_skeleton.txt
"""

import contextlib
import io as _io
import itertools
import json
import math
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

QS = (1, 2, 3, 4, 6)
INTERIOR = [1.0 / 3.0, 0.5, 2.0 / 3.0]
NPTS = 8001
GRID = np.linspace(1e-4, 1.0, NPTS)
CKPT = HERE / "outputs" / "holonomy_skeleton.partial.json"


def coeffs(tt):
    """(m, s, c) term table -> the C_q of F = sum_q C_q P_q, via Li_5(-z) = 2^-4 Li_5(z^2) - Li_5(z)."""
    C = {q: 0.0 for q in QS}
    for (m, s, c) in tt:
        c = int(c)
        if s > 0:
            C[c] = C.get(c, 0.0) + m
        else:
            C[c] = C.get(c, 0.0) - m
            C[2 * c] = C.get(2 * c, 0.0) + m / 16.0
    return C


def phase_basis(ys, nmax=4000):
    """P_q(y) = Re Li_5(e^{i q pi y}) on a grid, computed once."""
    n = np.arange(1, nmax + 1)
    w = 1.0 / n ** 5
    out = {}
    for q in QS:
        out[q] = (np.cos(np.outer(ys, q * math.pi * n)) * w).sum(axis=1)
    return out


def lattice_tables(nmax=2):
    ns = {"__file__": str(HERE / "ceiling_ilp.py"), "__name__": "ceiling_ilp"}
    with contextlib.redirect_stdout(_io.StringIO()):
        exec(compile((HERE / "ceiling_ilp.py").read_text(encoding="utf-8"),
                     "ceiling_ilp.py", "exec"), ns)
    TT = [ns["terms"](r, e, p) for (r, e, p) in ns["REPS"]]
    GA = ns["GAUGE"]
    for mult in itertools.product(range(nmax + 1), repeat=8):
        if not any(mult):
            continue
        t = list(GA)
        for j, mj in enumerate(mult):
            if mj:
                t += [(m * mj, s, c) for (m, s, c) in TT[j]]
        yield mult, t


def dist_to(pts, y):
    return min(abs(y - p) for p in pts)


def line(ch="-", n=94):
    print(ch * n)


def main():
    fails = []
    line("=")
    print("WHERE THE GLOBAL MINIMISER SITS, OVER THE WHOLE LATTICE")
    line("=")

    PB = phase_basis(GRID)

    # ---- C0: the fast path must agree with amin_closed_form ---------------------------------
    print("\n[0] CONTROL: THE FAST PATH AGAINST THE SLOW ONE")
    worst = 0.0
    tested = 0
    for mult, tt in lattice_tables(nmax=1):
        C = coeffs(tt)
        fast = sum(C[q] * PB[q] for q in QS)
        slow = sum(m * basis(GRID, s, c) for m, s, c in tt)
        worst = max(worst, float(np.abs(fast - slow).max()))
        tested += 1
        if tested >= 8:
            break
    ok = worst < 1e-6
    print("    %d contents compared point by point on the whole grid" % tested)
    print("    worst absolute difference : %.2e" % worst)
    print("\n   C0  the linear-combination shortcut is exact ...... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C0")
        print("       refusing to use it.")
        line("=")
        print("VERDICT: ABORTED -- the shortcut does not reproduce the potential.")
        line("=")
        return 1

    # ---- the scan ---------------------------------------------------------------------------
    print("\n[1] THE SCAN")
    t0 = time.time()
    ew = at_one = 0
    interior = []
    total = 0
    for mult, tt in lattice_tables(nmax=2):
        C = coeffs(tt)
        v = sum(C[q] * PB[q] for q in QS)
        i = int(np.argmin(v))
        y = float(GRID[i])
        total += 1
        if y > 0.995:
            at_one += 1
        elif y < 0.20:
            ew += 1
        else:
            interior.append(y)
        if total % 500 == 0:
            CKPT.write_text(json.dumps({"done": total, "ew": ew, "at_one": at_one,
                                        "interior": interior}))
    dt = time.time() - t0
    CKPT.write_text(json.dumps({"done": total, "ew": ew, "at_one": at_one,
                                "interior": interior, "finished": True}))

    # THE ELAPSED TIME GOES TO stderr, NOT HERE.  An archived output that carries a wall-clock
    # reading can never reproduce byte for byte, so check_reproduces.py reported this file as
    # broken every run -- for a number that is not a result.  Whatever a gate compares must be
    # a function of the inputs alone.  [[stale-outputs-lie]]
    sys.stderr.write("    (scan took %.1f s)\n" % dt)
    print("\n    contents scanned                       : %d" % total)
    print("    minimiser at the electroweak end y<0.2 : %d" % ew)
    print("    minimiser at the far symmetric point   : %d" % at_one)
    print("    minimiser strictly INTERIOR            : %d" % len(interior))
    print("""
    The first two are already the skeleton -- y = 0 is where the closed form works and y = 1 is
    the point eq:truevac compares against.  The question is entirely about the third group.""")

    if len(interior) < 20:
        print("\n    too few interior minimisers to say anything.  Not a result.")
        line("=")
        print("VERDICT: INCONCLUSIVE -- the sample is too small.")
        line("=")
        return 2

    arr = np.array(interior)
    d = np.array([dist_to(INTERIOR, y) for y in arr])
    print("\n    interior minimisers: %d" % len(arr))
    print("    distance to the nearest of {1/3, 1/2, 2/3}:  mean %.5f  median %.5f  max %.5f"
          % (d.mean(), np.median(d), d.max()))

    rng = np.random.default_rng(20260822)
    null = rng.uniform(arr.min(), arr.max(), size=400000)
    dn = np.array([dist_to(INTERIOR, y) for y in null])
    print("    NULL, y uniform on [%.3f, %.3f]         :  mean %.5f  median %.5f"
          % (arr.min(), arr.max(), dn.mean(), np.median(dn)))

    ratio = dn.mean() / d.mean() if d.mean() > 0 else float("inf")
    print("\n    the minimisers are %.2f times closer to the skeleton than chance" % ratio)
    ok = ratio > 1.5
    print("\n   C1  interior minimisers crowd the skeleton ........ %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C1")
        print("       -- and that is the honest outcome if they do not.")

    print("\n    where they actually are:")
    hist, edges = np.histogram(arr, bins=24, range=(0.2, 1.0))
    # bars are scaled to the tallest bin.  A fixed cap would clip 365 and 74 to the same
    # width and draw a flat distribution that is not there -- the picture has to be audited
    # like any other claim.
    top = int(hist.max())
    for k in range(len(hist)):
        lo, hi = edges[k], edges[k + 1]
        mark = ""
        for p, nm in ((1.0 / 3, "1/3"), (0.5, "1/2"), (2.0 / 3, "2/3")):
            if lo <= p < hi:
                mark = "   <- %s" % nm
        bar = "#" * int(round(46.0 * hist[k] / top)) if top else ""
        print("     %.3f-%.3f  %-6d %-46s%s" % (lo, hi, hist[k], bar, mark))

    # ---- C2: the sharper falsifier ----------------------------------------------------------
    # The ratio in C1 is one number and it can be dulled by where the mass happens to sit.  The
    # local question does not average: if the skeleton attracted minimisers at all, each of the
    # three bins containing 1/3, 1/2, 2/3 would stand above BOTH of its neighbours.  This is the
    # test that would have caught a real effect that C1's mean washed out, and it is the test a
    # reader who distrusts the null can check off the printed histogram by eye.
    print("\n[2] IS ANY SKELETON BIN A LOCAL MAXIMUM ?")
    peaks = 0
    for p, nm in ((1.0 / 3, "1/3"), (0.5, "1/2"), (2.0 / 3, "2/3")):
        k = int(np.searchsorted(edges, p, side="right")) - 1
        k = min(max(k, 1), len(hist) - 2)
        lo_n, hi_n = int(hist[k - 1]), int(hist[k + 1])
        is_peak = hist[k] > lo_n and hist[k] > hi_n
        peaks += int(is_peak)
        print("     %-4s bin %-6d  neighbours %-6d %-6d   local maximum: %s"
              % (nm, hist[k], lo_n, hi_n, "YES" if is_peak else "no"))
    ok2 = peaks >= 2
    print("\n   C2  at least two of the three are peaks .......... %s (%d of 3)"
          % ("PASS" if ok2 else "FAIL", peaks))
    if not ok2:
        fails.append("C2")

    # and say plainly where the mass that is NOT at the skeleton actually went
    edge = int((arr < 0.30).sum())
    print("\n    for the record, %d of %d interior minimisers (%.1f%%) sit in [0.20,0.30), against"
          % (edge, len(arr), 100.0 * edge / len(arr)))
    print("    the %.1f%% a uniform draw would put there.  The excess is at the low boundary of the"
          % (100.0 * (0.30 - arr.min()) / (arr.max() - arr.min())))
    print("    interior window, next to the electroweak end -- not at any of the three rationals.")

    line("=")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line("=")
        return 1
    print("VERDICT: the interior minimisers are not spread over the interval; they crowd the")
    print("         rational holonomies the charge alphabet selects.")
    line("=")
    return 0


if __name__ == "__main__":
    sys.exit(main())
