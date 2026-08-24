#!/usr/bin/env python3
"""The closed form run on Haba-Yamashita's general SU(3) formula -- a PREDICTION bank.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Haba & Yamashita, JHEP 05 (2004) 059 (hep-ph/0401185), give the one-loop Wilson-line
potential of 5D SU(N) on S^1/Z_2 in closed form but never locate its minimum -- their
summary calls analysing the vacuum structure the hard part.  Their eq. (3.20), for the
SU(3) model with one Wilson line phase `a`, N_a^(+-) adjoint and N_f^(+-) fundamental
Dirac fermions and N_s^(+-) complex scalars in the bulk, reads

  V = C sum_n n^-5 [ (-3/2 + 2 Na+) cos(2 pi n a)
                   + 2 Na- cos(pi n (2a-1))
                   + (-3 + 4 Na+ - Ns+ + 2 Nf+) cos(pi n a)
                   + (4 Na- - Ns- + 2 Nf-) cos(pi n (a-1)) ]

Since cos(pi n (a-1)) = (-1)^n cos(pi n a) and cos(pi n (2a-1)) = (-1)^n cos(2 pi n a),
this is exactly our (m, s, c) table with c in {1,2} and s = +-1, and the closed form of
../part_vii/amin_closed_form.py applies verbatim.  Nothing here is transcribed from a
representation-theory reconstruction of ours: it is their published counting rule.

What is predicted is alpha_min and the curvature F''(alpha_min).  Converting F'' into a
Higgs mass needs the model's own normalisation (Komori-Maru supply theirs in their
eqs. (80),(82); Haba-Yamashita do not), so no m_h is quoted here.
"""
import itertools
import json
import math
import pathlib

P = lambda *a: print(*a, flush=True)
OUT = pathlib.Path(__file__).resolve().parent / "outputs"
Z3 = 1.2020569031595942854
H4 = 1.0 + 1 / 2 + 1 / 3 + 1 / 4
LN2 = math.log(2.0)


def hy_table(Nap, Nam, Nfp, Nfm, Nsp, Nsm):
    """their eq. (3.20), as (m, s, c)."""
    return [(-1.5 + 2 * Nap, +1, 2),
            (2.0 * Nam, -1, 2),
            (-3.0 + 4 * Nap - Nsp + 2 * Nfp, +1, 1),
            (4.0 * Nam - Nsm + 2 * Nfm, -1, 1)]


def moments(tab):
    A2 = B2 = A4 = B4 = A4L = 0.0
    for m, s, c in tab:
        c = float(c)
        if s > 0:
            A2 += m * c ** 2
            A4 += m * c ** 4
            A4L += m * c ** 4 * math.log(c)
        else:
            B2 += m * c ** 2
            B4 += m * c ** 4
    return dict(D=A2 - 0.75 * B2, A4=A4, G=A4 * H4 - A4L - LN2 * B4)


def alpha_closed(mo):
    D, A4, G = mo["D"], mo["A4"], mo["G"]
    if D <= 0 or A4 <= 0:
        return None
    x = 0.05
    for _ in range(600):
        den = 4 * G - A4 * (4 * math.log(x) + 1)
        if den <= 0:
            return None
        xn = math.sqrt(24 * Z3 * D / den)
        if abs(xn - x) < 1e-15:
            x = xn
            break
        x = xn
    return x / math.pi


def Fpp(mo, a):
    x = math.pi * a
    return math.pi ** 2 * (2 * Z3 * mo["D"] - mo["A4"] * x * x / 6.0)


# ---------------------------------------------------------------- controls first
P("=" * 92)
P("CONTROLS -- their own statements, before any prediction")
P("=" * 92)
tab = hy_table(0, 0, 0, 0, 0, 0)
mo = moments(tab)
P("  pure gauge (no bulk matter): their eq. (3.20) gives (m,s,c) = %s" %
  [(round(m, 3), s, c) for m, s, c in tab])
P("     -> D = %.3f  %s" % (mo["D"], "NEGATIVE: no electroweak symmetry breaking, as it must be"
                            if mo["D"] < 0 else "POSITIVE -- CONTROL FAILED"))
P("  their eq. (3.10), one adjoint with eta*eta' = + and no gauge sector, is")
P("     (C/2)[cos(2 pi n a) + 2 cos(pi n a)] -- our table for it: %s" %
  [(m, s, c) for m, s, c in [(1, 1, 2), (2, 1, 1)]])
P("     their (3.20) at Na+ = 1, gauge part removed, gives 2*Na+ = 2 and 4*Na+ = 4,")
P("     which is 4x the single-d.o.f. values 1 and 2 of (3.10) -- their own text says the")
P("     adjoint Dirac fermion enters with coefficient 4, so the two readings agree.")

P("")
P("  BLIND DIRECTION, exact and free: the potential sees Nf and Ns only through 2Nf - Ns,")
_base = hy_table(1, 1, 0, 0, 0, 0)
P("     (Nf,Ns) -> (Nf+1, Ns+2) at either parity leaves the table IDENTICAL: %s / %s" %
  (_base == hy_table(1, 1, 1, 0, 2, 0), _base == hy_table(1, 1, 0, 1, 0, 2)))
P("     -- Part V's blind class, in a second model, straight off their counting rule.")

# ---------------------------------------------------------------- the bank
P("")
P("=" * 92)
P("PREDICTIONS -- alpha_min for every content of their SU(3) model up to 4 of each species")
P("=" * 92)
rows = []
R = range(0, 5)
for Nap, Nam, Nfp, Nfm, Nsp, Nsm in itertools.product(R, repeat=6):
    mo = moments(hy_table(Nap, Nam, Nfp, Nfm, Nsp, Nsm))
    if mo["D"] <= 0:
        continue
    a = alpha_closed(mo)
    if a is None or not (0 < a < 0.5):
        continue
    f = Fpp(mo, a)
    if f <= 0:
        continue
    rows.append((a, mo["D"], mo["A4"], f, (Nap, Nam, Nfp, Nfm, Nsp, Nsm)))

rows.sort()
P("contents with D > 0 and a genuine interior minimum: %d" % len(rows))
P("")
P("%10s %8s %8s %11s   %s" % ("alpha_min", "D", "A4", "F''(min)", "(Na+,Na-,Nf+,Nf-,Ns+,Ns-)"))
P("  -- the ten smallest alpha_min: the large-hierarchy corner --")
for a, D, A4, f, n in rows[:10]:
    P("%10.5f %8.2f %8.1f %11.4f   %s" % (a, D, A4, f, n))
P("  -- the ten largest, where the closed form degrades (alpha not small) --")
for a, D, A4, f, n in rows[-10:]:
    P("%10.5f %8.2f %8.1f %11.4f   %s" % (a, D, A4, f, n))

inrange = [r for r in rows if r[0] <= 0.12]
P("")
P("inside the regime where the expansion is good (alpha <= 0.12): %d of %d contents" %
  (len(inrange), len(rows)))
if inrange:
    P("smallest alpha_min = %.5f  ->  if m_W = 80.4 GeV then 1/R5 = %.0f GeV" %
      (inrange[0][0], 2 * 80.4 / inrange[0][0]))

P("")
P("D is quantised here too: the realised values are")
P("   %s" % sorted({round(r[1], 4) for r in rows})[:14])

OUT.mkdir(exist_ok=True)
(OUT / "hy_predictions.json").write_text(json.dumps(
    [dict(alpha=a, D=D, A4=A4, Fpp=f,
          content=dict(zip(("Nap", "Nam", "Nfp", "Nfm", "Nsp", "Nsm"), n)))
     for a, D, A4, f, n in rows[:60]], indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "hy_predictions.json"))
