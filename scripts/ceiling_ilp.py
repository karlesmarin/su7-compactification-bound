#!/usr/bin/env python3
"""The compactification-scale ceiling, as an integer program over the multiplet lattice.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

STATE.md section 4 carried  1/R5 <= 2.18 TeV  with the Higgs window imposed.  That number was
an ENUMERATION over contents of at most 6 multiplets, not a ceiling: nothing in the model caps
the content.  This file replaces it with a ceiling that holds for ARBITRARY bulk content.

The whole thing rests on two exact relations.  At the stationary point, with x = pi*alpha,

    F''(a_min) = pi^2 [ 2 z3 D - A4 x^2 / 6 ]        (mh_closed_form.py: the log cancels)
    m_h        = K sqrt(F'') / alpha ,  K = (sqrt3/2pi^3) m_W g4

so, writing  mu := m_h^2 / (K pi^2)^2  -- the Higgs mass in the units the potential works in --
and eliminating the logarithm against the fixed point  x^2 = 24 z3 D / [4G - A4(4 ln x + 1)] :

    (I)   ln x  =  (G - 3 mu) / A4  -  3/4
    (II)  x^2   =  12 z3 D / (6 mu + A4)

Neither has a logarithm left in it; both are verified below to 1e-12 on the five published rows.
(II) is the useful one: ONCE THE HIGGS MASS IS PINNED, alpha_min is algebraic in the two moments,
and  1/R5 = 2 pi m_W / x  is maximised by making D small and A4 large.  D cannot go below 1/8
(the odd-eighths theorem), so the ceiling is a statement about how large A4 can be while (I) still
admits a content -- and (I) is a LINEAR-FRACTIONAL constraint on the multiplicity vector.  That is
the integer program.

Sections:
  0  the lattice data, checked against moments() so it cannot silently desynchronise
  1  the two exact relations, on the five published rows
  2  the mod-6 law  8D = 2 A4 + 3 (mod 6), which is the odd-eighths theorem plus a mod-3 half
  3  the certified ceiling: an exact rational LP dual over the moment cone
  4  the ceiling as a function of the content size N -- the curve 2.18 TeV -> the asymptote
  5  controls, including the numerical minimisation of F for the champion content
  6  the ceiling once Part VI's escape must be affordable -- the bill moves the rung it lives on
"""
import itertools
import json
import math
import pathlib
from fractions import Fraction as Fr

import numpy as np

exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4          # K = (sqrt3/2pi^3) m_W g4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_LO, MH_HI = 125.0, 127.0                                  # their own window
MU_LO, MU_HI = MU(MH_LO), MU(MH_HI)
LN2, LN3, H4Q = math.log(2), math.log(3), Fr(25, 12)

# ---------------------------------------------------------------- 0  the lattice data
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]

_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]          # A4  (integer)
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]      # 8D  (integer)
GV = [moments([(r, e, p, 1)])["G"] - _g["G"] for r, e, p in REPS]                   # G
A0, K0, G0 = round(_g["A4"]), round(8 * _g["D"]), _g["G"]

# Part VI's escape needs a host, and the host has a price -- both read off the lattice above, so
# section 4 can count the contents that can pay it and section 6 can certify what it costs the
# ceiling.  Nothing here is typed: COST8 is the host's own 8D.
HOST = NAMES.index("84(+,+)")
COST8 = KV[HOST]
K_ESC = COST8 + 1 + (COST8 % 2)       # the smallest ODD rung strictly above the cost
assert K_ESC % 2 == 1 and K_ESC > COST8 and K_ESC - 2 <= COST8

# G is rational in H4=25/12, ln2, ln3: g = q*25/12 - r*ln3 - s*ln2.  Typed out so section 3 can
# work in exact arithmetic; CHECKED against the floats just built.
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
GQ_SYM0 = (-18, 0, Fr(-39, 2))
_sym = lambda t: float(t[0] * H4Q) - float(t[1]) * LN3 - float(t[2]) * LN2
assert all(abs(_sym(GQ_SYM[j]) - GV[j]) < 1e-12 for j in range(8)), "G symbolic/float desync"
assert abs(_sym(GQ_SYM0) - G0) < 1e-12, "G0 symbolic/float desync"

P("")
P("=" * 100)
P("0 -- THE LATTICE.  every content is  gauge + sum_j n_j * multiplet_j,  n in Z^8_{>=0}")
P("=" * 100)
P("%-10s %6s %6s %12s %10s" % ("multiplet", "A4", "8D", "G", "G/A4"))
for j in range(8):
    P("%-10s %6d %6d %12.5f %10.4f" % (NAMES[j], AV[j], KV[j], GV[j],
                                       GV[j] / AV[j] if AV[j] else float("nan")))
P("%-10s %6d %6d %12.5f" % ("gauge", A0, K0, G0))
P("all three moments are LINEAR in the multiplicities; A4 and 8D are integers.")

# ---------------------------------------------------------------- 1  the two relations
P("")
P("=" * 100)
P("1 -- THE TWO EXACT RELATIONS, on the five published rows (m_h is OUR m_h, not their column)")
P("=" * 100)
P("%-5s %10s %9s %10s | %12s %10s | %12s %10s" %
  ("row", "x", "m_h", "mu", "x from (I)", "err", "x from (II)", "err"))
rel = []
for label, content, a_them, mh_them, invR in T1:
    a, mo = closed_form(content)
    x, D, A4, G = math.pi * a, mo["D"], mo["A4"], mo["G"]
    fpp = math.pi ** 2 * (2 * Z3 * D - A4 * x * x / 6.0)
    mh = KK * math.sqrt(fpp) / a
    mu = MU(mh)
    xI = math.exp((G - 3 * mu) / A4 - 0.75)
    xII = math.sqrt(12 * Z3 * D / (6 * mu + A4))
    P("%-5s %10.6f %9.4f %10.4f | %12.6f %10.1e | %12.6f %10.1e" %
      (label, x, mh, mu, xI, abs(xI / x - 1), xII, abs(xII / x - 1)))
    rel.append(max(abs(xI / x - 1), abs(xII / x - 1)))
P("")
P("max relative error over the five rows: %.1e  -- these are identities, not fits." % max(rel))

# ---------------------------------------------------------------- 2  the mod-6 law
P("")
P("=" * 100)
P("2 -- THE MOD-6 LAW:   8D = 2 A4 + 3   (mod 6)   for every bulk content")
P("=" * 100)
P("Proof.  Per term (m, s, c) of the potential, with c and m integers for matter:")
P("   A4 gets  m c^4 (s=+1 only);   8D = 8A2 - 6B2  gets  8 m c^2 (s=+1) or -6 m c^2 (s=-1).")
P("   mod 3:  s=-1 gives -6mc^2 = 0;  s=+1 gives  m c^2 (8 + c^2) = m c^2 (c^2 + 2) = 0,")
P("           because c^4 = c^2 (mod 3) for EVERY integer c (Fermat).  So 8D + A4 = 0 (mod 3).")
P("   mod 2:  matter contributes 8A2 - 6B2, even; the gauge sector contributes -27, odd.")
P("   Together:  8D + A4 = 0 (mod 3)  and  8D odd  =>  8D = 2A4 + 3 (mod 6).")
P("")
P("Check on the five published rows and on the generators:")
P("   %-12s %6s %6s %10s %10s" % ("content", "8D", "A4", "8D+A4 %3", "8D %2"))
for label, content, a_them, mh_them, invR in T1:
    mo = moments(content)
    k, t = round(8 * mo["D"]), round(mo["A4"])
    P("   %-12s %6d %6d %10d %10d" % (label, k, t, (k + t) % 3, k % 2))
_bad = [(NAMES[j], AV[j], KV[j]) for j in range(8) if (AV[j] + KV[j]) % 3]
P("   generators violating  a + 8d = 0 (mod 3): %s" % (_bad or "NONE (all eight)"))
P("   gauge: A4 = %d, 8D = %d, sum mod 3 = %d, 8D mod 2 = %d" % (A0, K0, (A0 + K0) % 3, K0 % 2))
rng = np.random.default_rng(20260808)
_viol = 0
for _ in range(20000):
    n = rng.integers(0, 9, size=8)
    t, k = A0 + int(n @ AV), K0 + int(n @ KV)
    if (k + t) % 3 or k % 2 == 0:
        _viol += 1
P("   random control: 20000 contents with multiplicities in [0,8], violations = %d" % _viol)
P("")
P("The odd-eighths theorem of section 3 is the mod-2 half of this.  The mod-3 half is new: it")
P("removes two thirds of the (A4, 8D) plane, and the ceiling below feels it -- the real LP optimum")
P("falls at (216, 1), where 216 + 1 = 217 is not divisible by 3, so no content sits there.  The")
P("effect on the number is small (control (c) below measures it: 7 GeV); the effect on the")
P("LATTICE is not, and it is what makes the ceiling a statement about contents rather than about")
P("a relaxation.")

# ---------------------------------------------------------------- 3  the certified ceiling
P("")
P("=" * 100)
P("3 -- THE CEILING.  exact rational LP dual over the moment cone")
P("=" * 100)
P("For a content at (A4, 8D) = (t, k) whose Higgs mass sits at the top of the window, (II) fixes")
P("x, hence 1/R5, with no freedom left.  (I) then DEMANDS a specific value of G:")
P("")
P("     G*(t,k) = t (ln x + 3/4) + 3 mu ,      x = sqrt( 12 z3 (k/8) / (6 mu + t) )")
P("")
P("so the pair (t,k) is realisable only if some content at (t,k) has G <= G*(t,k).  Relaxing the")
P("multiplicities to the reals turns 'the smallest G at (t,k)' into a two-constraint LP whose DUAL")
P("has only two variables -- so its vertices can be enumerated EXACTLY in rationals.")

# ln2 and ln3 rounded in the direction that keeps every GQ[j] a LOWER bound on the true g_j, so
# the LP bound built from them stays a bound.  30-digit brackets from mpmath, verified below.
import mpmath as mp
mp.mp.dps = 60
_DEN = 10 ** 30
_lo = lambda v: Fr(int(mp.floor(v * _DEN)), _DEN)
_hi = lambda v: Fr(int(mp.ceil(v * _DEN)), _DEN)
LN2_LO, LN2_HI = _lo(mp.log(2)), _hi(mp.log(2))
LN3_LO, LN3_HI = _lo(mp.log(3)), _hi(mp.log(3))
_bnd = lambda q, r, s: (Fr(q) * H4Q - Fr(r) * (LN3_HI if r >= 0 else LN3_LO)
                        - Fr(s) * (LN2_HI if s >= 0 else LN2_LO))
GQ = [_bnd(*GQ_SYM[j]) for j in range(8)]
G0Q = _bnd(*GQ_SYM0)
_mpq = lambda v: mp.mpf(Fr(v).numerator) / Fr(v).denominator
_exact = lambda q, r, s: _mpq(q) * mp.mpf(25) / 12 - _mpq(r) * mp.log(3) - _mpq(s) * mp.log(2)
assert all(Fr(GQ_SYM[j][2]) >= 0 and Fr(GQ_SYM[j][1]) >= 0 for j in range(8))
for j in range(8):
    assert mp.mpf(GQ[j].numerator) / GQ[j].denominator <= _exact(*GQ_SYM[j]), NAMES[j]
assert mp.mpf(G0Q.numerator) / G0Q.denominator <= _exact(*GQ_SYM0)

VERTS = []
for i, j in itertools.combinations(range(8), 2):
    det = AV[i] * KV[j] - AV[j] * KV[i]
    if det == 0:
        continue
    lam = Fr(GQ[i] * KV[j] - GQ[j] * KV[i], det)
    nu = Fr(AV[i] * GQ[j] - AV[j] * GQ[i], det)
    if all(lam * AV[m] + nu * KV[m] <= GQ[m] for m in range(8)) and (lam, nu) not in VERTS:
        VERTS.append((lam, nu))
P("")
P("dual feasible region  { (lam,nu) : lam*a_j + nu*(8d)_j <= g_j  for all eight j } -- vertices:")
for lam, nu in VERTS:
    tight = [NAMES[m] for m in range(8) if lam * AV[m] + nu * KV[m] == GQ[m]]
    P("   lam = %-22s nu = %-22s   tight on %s" % (str(lam)[:22], str(nu)[:22], ", ".join(tight)))
P("   (verified in exact rational arithmetic; ln2 and ln3 rounded so the bound stays a bound)")


def gmin_cone(t, k):
    """exact lower bound on G over the real relaxation at (A4,8D)=(t,k); None if unreachable."""
    T, Q = t - A0, k - K0
    if T < 0 or Q > 8 * T:                       # the moment cone is generated by (0,-6),(1,8)
        return None
    return G0Q + max(l * T + n * Q for l, n in VERTS)


def x_of(t, k, mu):
    return math.sqrt(12 * Z3 * (k / 8.0) / (6 * mu + t))


def gstar(t, k, mu):
    return t * (math.log(x_of(t, k, mu)) + 0.75) + 3 * mu


def tmax(k, mu, slack=400):
    """largest A4 at which the cone still admits a content, respecting the mod-6 law.

    The whole range is scanned -- no early break -- because nothing proves the feasible set of
    A4 is an interval; and the scan is only trusted if it ends in `slack` consecutive failures.
    """
    hi = max(4000, 80 * k)
    out = None
    for t in range(hi):
        if (t + k) % 3:
            continue
        gm = gmin_cone(t, k)
        if gm is not None and float(gm) <= gstar(t, k, mu):
            out = t
    assert out is None or out < hi - slack, "scan bound too low at 8D=%d" % k
    return out


P("")
P("%8s %10s %12s %14s %14s %10s" % ("8D", "max A4", "alpha_min", "1/R5 (GeV)", "m_h forced", "TeV"))
ceil_rows, ceiling = [], None
for k in [1, 3, 5, 7, 9, 15, 21, 33, 45, 65, 99, 129, 201, 301, 501]:
    t = tmax(k, MU_HI)
    if t is None:
        P("%8d   no content can sit in the window at any A4" % k)
        continue
    x = x_of(t, k, MU_HI)
    inv = 2 * math.pi * MW / x
    ceil_rows.append(dict(k8D=k, A4=t, alpha=x / math.pi, invR=inv))
    # the TeV column is not decoration: the paper quotes the per-D ceiling in TeV (10.03, 6.27,
    # 4.21, 3.14, 2.80) and printing only GeV left those renderings unbacked in the archive.
    P("%8d %10d %12.6f %14.1f %14.1f %10.2f" % (k, t, x / math.pi, inv, MH_HI, inv / 1000))
    if ceiling is None or inv > ceiling[0]:
        ceiling = (inv, t, k)
P("")
P("CEILING, valid for arbitrary bulk content:  1/R5 <= %.0f GeV = %.2f TeV" % (ceiling[0], ceiling[0] / 1000))
P("   attained at  A4 = %d,  D = %s = %.4f,  alpha_min = %.5f" %
  (ceiling[1], "%d/8" % ceiling[2], ceiling[2] / 8, x_of(ceiling[1], ceiling[2], MU_HI) / math.pi))
_mono = all(ceil_rows[i]["invR"] > ceil_rows[i + 1]["invR"] for i in range(len(ceil_rows) - 1))
P("   the per-D ceiling is monotone decreasing in D: %s.  So the maximum sits on the QUANTUM" % _mono)
P("   D = 1/8 -- the hierarchy is generated by sitting on the smallest curvature the lattice has,")
P("   which is the odd-eighths theorem doing physical work.  At large D it flattens out towards")
P("   ~%.1f TeV; D large is not what kills the ceiling, it is what makes it ordinary." %
  (ceil_rows[-1]["invR"] / 1000))

# ---------------------------------------------------------------- 4  the size curve
P("")
P("=" * 100)
P("4 -- WHAT THE CAP WAS WORTH.  the ceiling as a function of the content size N")
P("=" * 100)
AVn, KVn, GVn = np.array(AV, float), np.array(KV, float), np.array(GV)


def alpha_closed_vec(D, A4, G):
    """the fixed point x^2 = 24 z3 D/[4G - A4(4 ln x + 1)], vectorised; nan where it fails."""
    x = np.full(D.shape, 0.05)
    ok = D > 0
    for _ in range(400):
        den = 4 * G - A4 * (4 * np.log(x) + 1)
        good = ok & (den > 0)
        xn = np.where(good, np.sqrt(np.abs(24 * Z3 * D / np.where(good, den, 1))), np.nan)
        if np.nanmax(np.abs(xn - x)) < 1e-14:
            x = xn
            break
        x = xn
    return x


def sweep(N):
    """every content of at most N multiplets, by the closed form."""
    vecs = np.zeros((1, 8), dtype=np.int32)
    seen = {(0,) * 8}
    cur = [np.zeros(8, dtype=np.int32)]
    allv = []
    for _ in range(N):
        nxt, nl = [], []
        for v in cur:
            for j in range(8):
                w = v.copy()
                w[j] += 1
                tk = tuple(w)
                if tk not in seen:
                    seen.add(tk)
                    nxt.append(w)
                    nl.append(w)
        allv.extend(nl)
        cur = nxt
    V = np.array(allv, dtype=np.int32)
    D = (K0 + V @ KVn) / 8.0
    A4 = A0 + V @ AVn
    G = G0 + V @ GVn
    x = alpha_closed_vec(D, A4, G)
    with np.errstate(invalid="ignore"):
        fpp = math.pi ** 2 * (2 * Z3 * D - A4 * x * x / 6.0)
        mh = KK * np.sqrt(np.where(fpp > 0, fpp, np.nan)) * math.pi / x
        inv = 2 * math.pi * MW / x
    win = np.isfinite(mh) & (mh >= MH_LO) & (mh <= MH_HI)
    return V, D, A4, inv, mh, win


P("The last three columns are section 6's: how many of the window contents hold Part VI's host, how")
P("many can still afford to donate it (8D >= %d), and the best hierarchy among those that can." % K_ESC)
P("%4s %12s %10s %6s %7s %11s %10s %s" %
  ("N", "contents", "in window", "host", "can pay", "best 1/R5", "D", "content"))
curve = []
NMAX = 14
for N in range(4, NMAX + 1):
    V, D, A4, inv, mh, win = sweep(N)
    if not win.any():
        P("%4d %12d %10d %6s %7s %11s" % (N, len(V), 0, "--", "--", "--"))
        continue
    K8 = np.rint(8 * D).astype(int)
    host = win & (V[:, HOST] >= 1)
    pays = host & (K8 >= K_ESC)
    i = int(np.nanargmax(np.where(win, inv, -np.inf)))
    txt = " + ".join("%dx%s" % (m, NAMES[j]) for j, m in enumerate(V[i]) if m)
    ip = int(np.nanargmax(np.where(pays, inv, -np.inf))) if pays.any() else None
    P("%4d %12d %10d %6d %7d %11.0f %10s %s" %
      (N, len(V), int(win.sum()), int(host.sum()), int(pays.sum()), inv[i],
       "%d/8" % K8[i], txt))
    curve.append(dict(N=N, contents=int(len(V)), in_window=int(win.sum()),
                      invR=float(inv[i]), D=float(D[i]), A4=float(A4[i]),
                      content=[int(z) for z in V[i]],
                      with_host=int(host.sum()), can_pay=int(pays.sum()),
                      invR_paying=(float(inv[ip]) if ip is not None else None),
                      D_paying=(float(D[ip]) if ip is not None else None),
                      content_paying=([int(z) for z in V[ip]] if ip is not None else None)))
best = curve[-1]
P("")
P("The 2.18 TeV in STATE.md section 4 is the N=6 row.  It is a ceiling on the SIZE OF THE CONTENT,")
P("not on the model: at N=%d the same window already allows %.2f TeV, a factor %.2f more, and the" %
  (best["N"], best["invR"] / 1000, best["invR"] / 2182.0))
P("certified sup over all N is %.2f TeV." % (ceiling[0] / 1000))

# ---------------------------------------------------------------- 5  controls
P("")
P("=" * 100)
P("5 -- CONTROLS")
P("=" * 100)
P("(a) EVERY champion of the curve, re-minimised numerically -- the closed form is not trusted to")
P("    police itself.  m_h at the numeric alpha must still land inside the window.")
P("    %-52s %10s %10s %8s %10s" % ("content", "closed", "numeric", "err %", "m_h(num)"))
ctrl_ok = True
for row in curve:
    champ = [(REPS[j][0], REPS[j][1], REPS[j][2], m) for j, m in enumerate(row["content"]) if m]
    a_num = numeric_min(champ)
    a_cf, mo = closed_form(champ)
    fpp_n = math.pi ** 2 * (2 * Z3 * mo["D"] - mo["A4"] * (math.pi * a_num) ** 2 / 6.0)
    mh_n = KK * math.sqrt(fpp_n) / a_num if fpp_n > 0 else float("nan")
    txt = " + ".join("%dx%s" % (m, NAMES[j]) for j, m in enumerate(row["content"]) if m)
    P("    %-52s %10.6f %10.6f %+8.3f %10.4f%s" %
      (txt[:52], a_cf, a_num, 100 * (a_cf - a_num) / a_num, mh_n,
       "" if MH_LO - 0.5 <= mh_n <= MH_HI + 0.5 else "   *** OUT OF WINDOW ***"))
    row["alpha_numeric"], row["mh_numeric"] = a_num, mh_n
    ctrl_ok &= MH_LO - 0.5 <= mh_n <= MH_HI + 0.5
P("    every champion survives its own numerical minimisation: %s" % ctrl_ok)

P("")
P("(b) the ceiling must dominate every content actually found in the window:")
V, D, A4, inv, mh, win = sweep(NMAX)
worst = float(np.nanmax(np.where(win, inv, -np.inf)))
P("    max 1/R5 over the %d window contents at N<=%d : %.0f GeV" % (int(win.sum()), NMAX, worst))
P("    certified ceiling                              : %.0f GeV   -> %s" %
  (ceiling[0], "HOLDS" if worst <= ceiling[0] else "*** VIOLATED ***"))

P("")
P("(c) a control that must FAIL -- drop the mod-3 law and the ceiling lands on a phantom point:")
t_ph = ceiling[1] + 1
P("    (A4,8D) = (%d,%d):  cone bound G >= %.3f, G*(127) = %.3f  -> LP says feasible, but" %
  (t_ph, ceiling[2], float(gmin_cone(t_ph, ceiling[2]) or 0), gstar(t_ph, ceiling[2], MU_HI)))
P("    (A4 + 8D) mod 3 = %d != 0, so NO content sits there.  Ignoring it inflates the ceiling to" %
  ((t_ph + ceiling[2]) % 3))
P("    %.0f GeV." % (2 * math.pi * MW / x_of(t_ph, ceiling[2], MU_HI)))

P("")
P("(d) the ceiling scales with the Higgs window, as it must (1/R5 ~ m_h at fixed D):")
for mh_try in (125.0, 127.0, 200.0):
    mu = MU(mh_try)
    t = tmax(1, mu)
    P("    m_h <= %5.0f GeV  ->  max A4 = %5d,  1/R5 <= %8.0f GeV" %
      (mh_try, t, 2 * math.pi * MW / x_of(t, 1, mu)))

# ------------------------------------------------- 6  the ceiling a content that pays Part VI gets
P("")
P("=" * 100)
P("6 -- WHAT PART VI's ESCAPE COSTS THE CEILING")
P("=" * 100)
P("Part VI's surviving escape from proton decay is a family-dependent charge, and it must be hosted")
P("by an 84(+,+) donated to the brane -- which removes it from the Higgs potential.  Section 0's")
P("lattice already prices that: the host carries 8D = %+d, so donating it takes %d/8 off D." %
  (COST8, COST8))
P("")
P("Electroweak symmetry breaks only for D > 0, so a content can pay and survive iff  8D > %d." % COST8)
P("8D is ODD (section 2), so the condition is an integer one with no slack in it:")
P("")
P("     a content can afford the escape  <=>  8D >= %d" % K_ESC)
P("")
P("and the per-rung ceiling of section 3 is monotone decreasing in D, so the escape does not merely")
P("forbid contents -- it moves the ceiling, by forbidding the rungs the ceiling lives on.")
P("")
P("%8s %10s %12s %14s" % ("8D", "max A4", "alpha_min", "1/R5 (GeV)"))
esc_rows, esc_ceiling = [], None
K_SCAN = list(range(K_ESC, 62, 2)) + [65, 71, 79, 91, 99, 111, 129, 151, 175, 201, 251, 301, 401, 501]
for k in K_SCAN:
    t = tmax(k, MU_HI)
    if t is None:
        P("%8d   no content can sit in the window at any A4" % k)
        continue
    x = x_of(t, k, MU_HI)
    inv = 2 * math.pi * MW / x
    esc_rows.append(dict(k8D=k, A4=t, alpha=x / math.pi, invR=inv))
    if k <= 23 or k in (33, 51, 99, 201, 501):
        P("%8d %10d %12.6f %14.1f" % (k, t, x / math.pi, inv))
    if esc_ceiling is None or inv > esc_ceiling[0]:
        esc_ceiling = (inv, t, k)
_emono = all(esc_rows[i]["invR"] > esc_rows[i + 1]["invR"] for i in range(len(esc_rows) - 1))
P("   ... %d rungs scanned, from 8D = %d to %d; monotone decreasing throughout: %s" %
  (len(esc_rows), esc_rows[0]["k8D"], esc_rows[-1]["k8D"], _emono))
P("")
P("CEILING for any content that can AFFORD the escape:  1/R5 <= %.0f GeV = %.2f TeV" %
  (esc_ceiling[0], esc_ceiling[0] / 1000))
P("   attained at  A4 = %d,  D = %d/8,  alpha_min = %.5f" %
  (esc_ceiling[1], esc_ceiling[2], x_of(esc_ceiling[1], esc_ceiling[2], MU_HI) / math.pi))
P("   against %.0f GeV with the escape not required -- the escape divides the ceiling by %.2f." %
  (ceiling[0], ceiling[0] / esc_ceiling[0]))
P("")
P("The two halves of the argument come from different papers and neither states this: Part VI prices")
P("the escape and says in its own words that the modified potential is not recomputed there; Part VII")
P("certifies the ceiling and imposes no anomaly condition.  The quantum does the work twice -- the")
P("unconstrained ceiling sits ON 8D = 1, and the escape's bill is exactly what makes that rung")
P("unaffordable.")

P("")
P("(e) CONTROL -- the enumeration must respect it.  Every content of at most %d multiplets that" % NMAX)
P("    lands in the window AND holds a host AND still has D > 0 after donating it:")
V, D, A4, inv, mh, win = sweep(NMAX)
K8 = np.rint(8 * D).astype(int)
pays = win & (V[:, HOST] >= 1) & (K8 >= K_ESC)
P("    %d of the %d window contents can afford it (%d hold a host at all)" %
  (int(pays.sum()), int(win.sum()), int((win & (V[:, HOST] >= 1)).sum())))
if pays.any():
    i = int(np.nanargmax(np.where(pays, inv, -np.inf)))
    txt = " + ".join("%dx%s" % (m, NAMES[j]) for j, m in enumerate(V[i]) if m)
    P("    best 1/R5 among them : %.0f GeV   8D = %d -> %d after donating   m_h = %.1f" %
      (inv[i], K8[i], K8[i] - COST8, mh[i]))
    P("    %s" % txt)
    P("    certified escape ceiling : %.0f GeV   -> %s" %
      (esc_ceiling[0], "HOLDS" if inv[i] <= esc_ceiling[0] else "*** VIOLATED ***"))
    esc_best = dict(invR=float(inv[i]), k8D=int(K8[i]), mh=float(mh[i]),
                    content=[int(z) for z in V[i]])
else:
    esc_best = None
_ch = curve[-1]
_chk, _chhost = round(8 * _ch["D"]), _ch["content"][HOST]
P("    the champion of the unconstrained ceiling sits at 8D = %d and holds %d host(s): %s" %
  (_chk, _chhost,
   "it is %d/8 short of the bill and cannot donate one" % (COST8 - _chk + 1) if _chk < K_ESC
   else "it can afford the escape after all -- the two ceilings do not separate at N=%d" % NMAX))
if _chk < K_ESC:
    P("    -- the content that generates the largest hierarchy is exactly the content that cannot pay.")

P("")
P("SCOPE, and it has to be said whole.  This bounds the content BEFORE the donation -- the row that")
P("pays.  It does NOT bound the world after: the post-donation bulk needs only 8D >= 1, which is the")
P("unconstrained ceiling again.  The statement is 'a bulk content in their Table 1's sense that can")
P("implement the escape by donating one of its own multiplets cannot generate more than %.2f TeV'," %
  (esc_ceiling[0] / 1000))
P("not 'a universe with the escape in it cannot'.  It also inherits Part VI's anchor caveat: it is")
P("our F, whose alpha_min column does not reproduce theirs.  What does NOT inherit it is the")
P("arithmetic: 8D >= %d is a statement about integers." % K_ESC)

OUT.mkdir(exist_ok=True)
(OUT / "ceiling_ilp.json").write_text(json.dumps(dict(
    ceiling_GeV=ceiling[0], ceiling_A4=ceiling[1], ceiling_8D=ceiling[2],
    mh_window=[MH_LO, MH_HI], g4=G4,
    dual_vertices=[[str(l), str(n)] for l, n in VERTS],
    per_k=ceil_rows, size_curve=curve,
    escape=dict(host="84(+,+)", cost8=COST8, min_8D=K_ESC,
                ceiling_GeV=esc_ceiling[0], ceiling_A4=esc_ceiling[1], ceiling_8D=esc_ceiling[2],
                ratio_to_unconstrained=ceiling[0] / esc_ceiling[0],
                per_k=esc_rows, best_enumerated=esc_best)), indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "ceiling_ilp.json"))
