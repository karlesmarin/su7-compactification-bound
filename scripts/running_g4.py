#!/usr/bin/env python3
"""running_g4.py -- g_4 stops being a choice and becomes a fixed point.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

CARLES'S POINT 7.  Every absolute scale in this paper passes through K = (sqrt3/2pi^3) m_W g_4, so
mu ~ g_4^-2 and the ceiling depends on it.  The paper uses g_4 = 0.63 and says, correctly, that it
is the Standard-Model SU(2)_L coupling run up to the compactification scale, while the sibling
anchor uses 0.653 at the weak scale.  Both are defensible, the difference is 3.7 %, and the choice
is left to the reader.

But g_4 need not be chosen at all.  The scale it should be evaluated at is 1/R_5 = 2 pi m_W / x,
and x depends on mu which depends on g_4.  Imposing

    g_4 = g_2( 1/R_5 ),        1/R_5 = 2 pi m_W / x,        x^2 = 12 zeta(3) D / (A_4 + 6 mu),
    mu  = ( m_h / (K pi^2) )^2,   K = (sqrt3 / 2 pi^3) m_W g_4

closes the loop and turns a convention into a self-consistent condition.  The running is weak --
under five e-folds -- so this cannot destroy the arithmetic; what it can do is say which of 0.63
and 0.653 the model itself picks, and whether the answer moves the ceiling's vertex.

The four-dimensional running is NOT re-typed here: inv() and the PDG inputs are imported from
running_gap.py, which already archives the same one-loop table, so the two cannot drift apart.

Sections:
  1  where 0.63 and 0.653 sit on the running curve -- the paper's open item, as a number
  2  the fixed point, and that it is a contraction
  3  the ceiling recomputed self-consistently: does the vertex move?
  4  controls, including one that must come out different

Run:  python running_g4.py
"""
import json
import math
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import running_gap as RG                       # noqa: E402  -- the running, already archived

exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW = 80.4
MH_LO, MH_HI = 125.0, 127.0
G4_PAPER, G4_AHMN = 0.63, 0.653


def g2_at(mu):
    """the Standard-Model SU(2)_L coupling at scale mu, one loop, from the PDG inputs at M_Z."""
    return math.sqrt(4 * math.pi / RG.inv(RG.A2_INV_MZ, RG.B2, mu))


def KK_of(g4):
    return math.sqrt(3.0) / (2 * math.pi ** 3) * MW * g4


def mu_of(mh, g4):
    return (mh / (KK_of(g4) * math.pi ** 2)) ** 2


# ================================================================= 1
P("=" * 100)
P("1 -- WHERE THE TWO VALUES SIT ON THE RUNNING CURVE")
P("=" * 100)
P("  running imported from running_gap.py: alpha_2^-1(M_Z) = %.3f, b_2 = %s, M_Z = %s GeV"
  % (RG.A2_INV_MZ, RG.B2, RG.MZ))
P("")
P("  %14s %14s %12s %10s" % ("scale", "alpha_2^-1", "g_2", "e-folds"))
for lab, s in (("M_Z", RG.MZ), ("1 TeV", 1e3), ("2 TeV", 2e3), ("3.97 TeV", 3967.0),
               ("9.22 TeV", 9218.0), ("10.03 TeV", 10034.0), ("100 TeV", 1e5)):
    P("  %14s %14.3f %12.5f %10.3f"
      % (lab, RG.inv(RG.A2_INV_MZ, RG.B2, s), g2_at(s), math.log(s / RG.MZ)))
P("")
P("  So the two numbers the paper leaves side by side are the SAME coupling at two scales:")
P("     %.3f  (\\cite{AHMN}, the weak scale)   against  g_2(M_Z) = %.5f" % (G4_AHMN, g2_at(RG.MZ)))
P("     %.3f  (this paper, 'run up to the compactification scale')" % G4_PAPER)
P("")


def scale_where(target):
    """the scale at which g_2 equals target, by bisection."""
    lo, hi = RG.MZ, 1e12
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if g2_at(mid) > target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


s63, s653 = scale_where(G4_PAPER), scale_where(G4_AHMN)
P("  g_2 = %.3f at %.4g GeV      g_2 = %.3f at %.4g GeV" % (G4_PAPER, s63, G4_AHMN, s653))
P("  and the ceiling itself sits at 9218-10034 GeV, where g_2 = %.5f to %.5f."
  % (g2_at(10034.0), g2_at(9218.0)))
P("")
P("  ==> 0.653 is the coupling at M_Z to four digits, and 0.63 is the coupling at the ceiling")
P("      rounded to two.  The paper's sentence is right; what follows makes it a fixed point")
P("      instead of a rounding.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- THE FIXED POINT")
P("=" * 100)
P("  Iterate  g_4 -> g_2( 2 pi m_W / x(g_4) )  at the ceiling's vertex, A_4 = 104, 8D = 1,")
P("  with m_h at the top of the window.  If it is a contraction the loop closes on itself.")
P("")


def x_of(t, k, mh, g4):
    return math.sqrt(12 * Z3 * (k / 8.0) / (6 * mu_of(mh, g4) + t))


def fixed_point_scaled(t, k, mh, fac, g0=G4_PAPER, n=200):
    """the same closure with the matching scale taken at fac / R_5 instead of 1 / R_5."""
    g = g0
    for _ in range(n):
        gn = g2_at(fac * 2 * math.pi * MW / x_of(t, k, mh, g))
        if abs(gn - g) < 1e-14:
            return gn, True
        g = gn
    return g, False


def fixed_point(t, k, mh, g0=G4_PAPER, n=60, verbose=False):
    g, hist = g0, []
    for i in range(n):
        x = x_of(t, k, mh, g)
        inv_r = 2 * math.pi * MW / x
        gn = g2_at(inv_r)
        hist.append((i, g, inv_r, gn))
        if abs(gn - g) < 1e-14:
            g = gn
            break
        g = gn
    return g, hist


gstar, hist = fixed_point(104, 1, MH_HI)
P("  %6s %14s %16s %14s" % ("iter", "g_4 in", "1/R5 (GeV)", "g_4 out"))
for i, g, inv_r, gn in hist[:6]:
    P("  %6d %14.9f %16.2f %14.9f" % (i, g, inv_r, gn))
if len(hist) > 6:
    P("  %6s %14s %16s %14s" % ("...", "", "", ""))
    i, g, inv_r, gn = hist[-1]
    P("  %6d %14.9f %16.2f %14.9f" % (i, g, inv_r, gn))
P("")
P("  fixed point : g_4* = %.9f   reached in %d iterations" % (gstar, len(hist)))
d0, d1 = fixed_point(104, 1, MH_HI, g0=0.50)[0], fixed_point(104, 1, MH_HI, g0=0.80)[0]
P("  starting from 0.50 and from 0.80 : %.9f and %.9f -- the same point" % (d0, d1))
assert abs(d0 - gstar) < 1e-9 and abs(d1 - gstar) < 1e-9, "the iteration is not a contraction"
h = 1e-4
slope = (fixed_point(104, 1, MH_HI, g0=gstar + h, n=1)[1][0][3]
         - fixed_point(104, 1, MH_HI, g0=gstar - h, n=1)[1][0][3]) / (2 * h)
P("  |d(map)/dg| at the fixed point : %.2e   (a contraction needs < 1)" % abs(slope))
assert abs(slope) < 1, "the map is not contracting -- the fixed point is not the limit"
P("")
P("  and the paper's 0.63 against it: %.5f vs %.5f, a difference of %.2f %%"
  % (G4_PAPER, gstar, 100 * abs(G4_PAPER - gstar) / gstar))

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- WHAT IT DOES TO THE CEILING")
P("=" * 100)
P("  mu ~ g_4^-2, so a smaller g_4 raises mu, and (II) then lowers x and RAISES 1/R5.  The")
P("  question is whether it raises it past a lattice step, because A_4 is quantised.")
P("")
P("  %-28s %10s %14s %14s %12s" % ("g_4", "value", "mu", "1/R5 at 104", "shift"))
base = None
for lab, g in (("this paper", G4_PAPER), ("self-consistent", gstar),
               ("\\cite{AHMN}, weak scale", G4_AHMN)):
    inv_r = 2 * math.pi * MW / x_of(104, 1, MH_HI, g)
    if base is None:
        base = inv_r
    P("  %-28s %10.5f %14.4f %14.1f %11.2f %%"
      % (lab, g, mu_of(MH_HI, g), inv_r, 100 * (inv_r - base) / base))
P("")
P("  So the self-consistent value moves the number at the SAME vertex by well under a percent.")
P("  Whether the VERTEX itself moves is a separate question and it is the one that matters,")
P("  because A_4 steps by three:")
P("")
_g = moments([])
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])
LN2, LN3 = math.log(2), math.log(3)
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
WV = [1, -1, 3, -3, 6, -6, 9, -9]
U_OFF, W_OFF = -19.5, -1.5


def best_at(Atgt, Ktgt, g4):
    """max logarithmic budget at (A_4, 8D) with W > 0 -- the section-ceiling search, exhaustive."""
    best = [None]
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, u, v, w, n):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            ww = w + m * WV[0] + W_OFF
            if ww <= 0:
                return
            cand = (u + m * GQ_SYM[0][2] + U_OFF) * LN2 + v * LN3
            if best[0] is None or cand > best[0]:
                best[0] = cand
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            n.append(c)
            rec(i + 1, a - c * AV[j], k - c * KV[j], u + c * GQ_SYM[j][2],
                v + c * GQ_SYM[j][1], w + c * WV[j], n)
            n.pop()

    rec(0, Atgt, Ktgt, 0, 0, 0, [])
    return best[0]


def vertex(g4):
    """the largest A_4 on 8D = 1 whose budget covers what (I) demands, with W > 0."""
    mu = mu_of(MH_HI, g4)
    for t in range(215, 59, -1):
        if (t + 1) % 3:
            continue
        x = math.sqrt(12 * Z3 * (1 / 8.0) / (6 * mu + t))
        need = float(25 * t) / 12 - (t * (math.log(x) + 0.75) + 3 * mu)
        have = best_at(t - A0, 1 - K0, g4)
        if have is not None and have >= need:
            return t, 2 * math.pi * MW / x
    return None, None


P("  %-28s %10s %10s %14s" % ("g_4", "value", "A_4 vertex", "1/R5 (GeV)"))
rows = []
for lab, g in (("this paper", G4_PAPER), ("self-consistent", gstar),
               ("\\cite{AHMN}, weak scale", G4_AHMN)):
    t, inv_r = vertex(g)
    rows.append((lab, g, t, inv_r))
    P("  %-28s %10.5f %10s %14.1f" % (lab, g, t, inv_r))
P("")
same = rows[0][2] == rows[1][2]
P("  the vertex moves between the paper's value and the self-consistent one : %s" % (not same))
P("  the vertex moves between the self-consistent value and the weak scale  : %s"
  % (rows[1][2] != rows[2][2]))

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- CONTROLS")
P("=" * 100)
P("  (i) THE RUNNING MUST BE THE ONE ALREADY ARCHIVED.  running_gap.py prints alpha_2^-1 at")
P("      five scales; this file imports its inv() rather than retyping it, so the check is that")
P("      the archived table is reproduced digit for digit:")
arch = (HERE / "outputs" / "running_gap.txt").read_text(encoding="utf-8", errors="replace")
ok = True
for lab, s in (("M_Z", RG.MZ), ("1 TeV", 1e3), ("2 TeV", 2e3), ("3.97 TeV", 3967.0),
               ("10.03 TeV", 10034.0)):
    v = "%.3f" % RG.inv(RG.A2_INV_MZ, RG.B2, s)
    hit = v in arch
    ok &= hit
    P("      %-12s alpha_2^-1 = %-10s in the archive : %s" % (lab, v, hit))
P("      CONTROL -- every value found in outputs/running_gap.txt : %s" % ok)
assert ok, "the running has drifted from the archived run"

P("")
P("  (ii) AND THE FIXED POINT MUST BE ABLE TO COME OUT SOMEWHERE ELSE.  Run the same iteration")
P("       with the running switched off -- g_2 frozen at M_Z -- and it must land on the weak")
P("       scale value instead, or the loop is not doing anything:")
_saved = RG.B2
try:
    RG.B2 = 0.0
    gfrozen = fixed_point(104, 1, MH_HI)[0]
finally:
    RG.B2 = _saved
P("       running on  : g_4* = %.5f" % gstar)
P("       running off : g_4* = %.5f   (must be g_2(M_Z) = %.5f)" % (gfrozen, g2_at(RG.MZ)))
assert abs(gfrozen - g2_at(RG.MZ)) < 1e-9 and abs(gfrozen - gstar) > 1e-3, \
    "switching the running off changes nothing -- the fixed point is not measuring it"

P("")
P("  (iii) THE NUMBERS AS THE PAPER WOULD PRINT THEM:")
P("       g_2(M_Z)                        : %.4f" % g2_at(RG.MZ))
P("       g_4 self-consistent at the vertex: %.4f" % gstar)
P("       the paper's value               : %.2f" % G4_PAPER)
P("       1/R5 at A_4 = 104, self-consistent : %.0f GeV" % rows[1][3])
P("       against the paper's value          : %.0f GeV" % rows[0][3])
P("       shift                              : %.2f %%"
  % (100 * (rows[1][3] - rows[0][3]) / rows[0][3]))
P("       at the weak-scale g_4, the vertex   : A_4 = %d" % rows[2][2])
P("       and the ceiling there              : %.0f GeV" % rows[2][3])
P("       contraction factor of the map      : %.1e" % abs(slope))
P("       scale at which g_2 = 0.63          : %.0f GeV" % s63)

# ================================================================= 5
P("")
P("=" * 100)
P("5 -- WHAT THE FIXED POINT DOES AND DOES NOT OPEN")
P("=" * 100)
P("  Carles asked whether this points at more couplings at other scales.  Two readings, and")
P("  both are measurable rather than matters of opinion.")
P("")
P("  (a) IS THE FIXED POINT UNIQUE IN THE PHYSICAL RANGE?  A contraction is local; the basin")
P("      is not.  Sweeping the starting value across everything a reader might try:")
P("")
P("      %10s %16s %16s" % ("g_4 start", "g_4* reached", "1/R5 (GeV)"))
lands = set()
for g0 in (0.30, 0.40, 0.50, 0.60, 0.6275, 0.63, 0.653, 0.75, 0.90, 1.10, 1.30):
    gf = fixed_point(104, 1, MH_HI, g0=g0)[0]
    lands.add(round(gf, 9))
    P("      %10.4f %16.9f %16.1f" % (g0, gf, 2 * math.pi * MW / x_of(104, 1, MH_HI, gf)))
P("")
P("      distinct fixed points reached : %d" % len(lands))
assert len(lands) == 1, "the map has more than one attractor in the physical range"
P("      One.  The map is g -> g_2(M_KK(g)) with g_2 monotone decreasing in the scale and M_KK")
P("      decreasing in g, so the composition is monotone and its slope is 6.5e-3: there is no")
P("      room for a second crossing.  (The formal g -> 0 solution needs M_KK -> infinity and is")
P("      not a physical branch.)")
P("")
P("  (b) DOES ANYTHING ELSE IN THE PAPER SIT AT A SCALE THAT IS ITSELF g_4-DEPENDENT?  Yes, one")
P("      thing: the gauge-coupling gap of section pred is quoted AT the ceiling, and the ceiling")
P("      moves with g_4.  So the same closure applies to it, and the size is worth having:")
P("")
P("      %-34s %12s %14s %14s" % ("evaluated at", "scale (GeV)", "a2^-1 - a3^-1", "shift"))
_g0 = None
for lab, s in (("the paper's ceiling, 10034", 10034.0),
               ("the constrained vertex, 9218", 9218.0),
               ("the same, self-consistent", rows[1][3])):
    gp = RG.inv(RG.A2_INV_MZ, RG.B2, s) - RG.inv(RG.A3_INV_MZ, RG.B3, s)
    if _g0 is None:
        _g0 = gp
    P("      %-34s %12.1f %14.3f %13.3f" % (lab, s, gp, gp - _g0))
P("")
P("      The gap moves by hundredths.  So the answer to the question is narrow and worth saying")
P("      narrowly: closing the loop removes ONE convention, it does not reveal a new coupling,")
P("      and the only other quantity that inherits the self-consistency barely feels it.")
P("")
P("  (c) AND THE QUESTION THAT MATTERS, WHICH IS WHETHER THE CONVENTION WAS REMOVED OR MOVED.")
P("      The closure fixes g_4 GIVEN that it is evaluated at 1/R_5.  But why 1/R_5?  The KK")
P("      masses are n/R_5, the fundamental domain has a 2 pi in it, and a matching scale is a")
P("      convention like any other.  If varying it moves g_4 by more than the 0.39 % the closure")
P("      buys, then nothing was removed -- the arbitrariness was relabelled, and saying")
P("      otherwise would be an overclaim.  So:")
P("")
P("      %-22s %13s %12s %10s %14s %10s %s" % ("matching scale", "scale (GeV)", "g_4*",
                                               "A_4 vertex", "1/R5 (GeV)", "vs 1/R5", "SM?"))
# EL BETA DE CUATRO DIMENSIONES SOLO VALE POR DEBAJO DEL PRIMER NIVEL KK.  Es lo que dice este
# mismo fichero al final, y lo que dice la S13 del articulo: por debajo de 1/R_5 la teoria es el
# Modelo Estandar.  2 pi/R_5 son 6.28 M_KK, con seis niveles abiertos por debajo, y ahi el
# running de la torre es de potencia y no logaritmico.  Los puntos de arriba se calculan y se
# imprimen ---son informativos--- pero NO entran en el sistematico.
for lab, fac in (("1/(2 R_5)", 0.5), ("1/R_5", 1.0), ("2/R_5", 2.0),
                 ("pi/R_5", math.pi), ("2 pi/R_5", 2 * math.pi)):
    g, _ = fixed_point_scaled(104, 1, MH_HI, fac)
    tv, rv = vertex(g)
    P("      %-22s %13.1f %12.6f %10s %14.1f %9.2f %% %s"
      % (lab, fac * 2 * math.pi * MW / x_of(104, 1, MH_HI, g), g, tv, rv,
         100 * (g - gstar) / gstar, "si" if fac <= 1.0 else "NO: por encima de M_KK"))
P("")
_lo, _ = fixed_point_scaled(104, 1, MH_HI, 0.5)
_hi, _ = fixed_point_scaled(104, 1, MH_HI, 2 * math.pi)
_spread_valido = 100 * abs(_lo - gstar) / gstar
_spread_todo = 100 * abs(_hi - _lo) / gstar
P("      EL TRAMO EN EL QUE ESTE RUNNING ES EL RUNNING, de 1/(2R_5) a 1/R_5:")
P("        spread en g_4                                    : %.2f %%" % _spread_valido)
P("      lo que la vuelta de tuerca gana sobre el 0.63 del paper : %.2f %%"
  % (100 * abs(G4_PAPER - gstar) / gstar))
P("")
_vt_valido = sorted({vertex(g)[0] for g in (_lo, gstar, G4_PAPER)})
_vt_todo = sorted({vertex(g)[0] for g in (_lo, gstar, _hi, G4_PAPER)})
P("      vertices dentro del tramo valido, g_4 en [%.4f, %.4f] : A_4 = %s"
  % (min(_lo, gstar), max(_lo, gstar), ", ".join(str(v) for v in _vt_valido)))
P("      vertices si se estira hasta 2 pi/R_5                  : A_4 = %s"
  % ", ".join(str(v) for v in _vt_todo))
P("")
P("      LO QUE ESTO MIDE, Y LO QUE NO.  Dentro del tramo en el que el articulo afirma que la")
P("      teoria es el Modelo Estandar, mover la escala de matching un factor dos cambia g_4 un")
P("      %.2f %% y NO mueve el vertice.  Eso es todo lo que un running cuatridimensional puede" % _spread_valido)
P("      decir aqui.  Estirarlo hasta 2 pi/R_5 daria %.4f y movería el vertice a A_4 = %d, pero"
  % (_hi, vertex(_hi)[0]))
P("      esa escala esta por encima de seis niveles KK: el numero queda impreso como marca de lo")
P("      que NO esta zanjado, no como sistematico.  De modo que g_4 sigue siendo un DATO ---la")
P("      escala de matching es un convenio y esto no lo quita--- pero el tamanyo de esa libertad")
P("      no queda establecido aqui.  Establecerlo pide el matching de la torre, que es la segunda")
P("      condicion del final de este fichero.")
P("")
P("      Lo que la vuelta de tuerca SI establece, y merece quedarse: dada la escala de matching,")
P("      el valor tampoco es libre ---es un punto fijo que contrae fuerte, no un ajuste--- y los")
P("      dos numeros que el articulo deja uno al lado del otro son un mismo acoplo a dos escalas.")
P("")
# EL CONTROL, y tiene dos mitades porque la conclusion tiene dos mitades.  Si el tramo valido
# empezara a mover el vertice, o si la extension dejara de moverlo, la frase de arriba ya no
# seria la que describe estos numeros.  [[a-control-that-cannot-fail]]
assert len(_vt_valido) == 1, (
    "el tramo valido ya mueve el vertice: la conclusion de arriba hay que reescribirla")
assert len(_vt_todo) > 1, (
    "estirar hasta 2 pi/R_5 ya no mueve el vertice: entonces la razon para no citarlo es solo "
    "que el beta no vale ahi, y la frase de arriba dice dos cosas donde ahora hay una")
P("      AND A SECOND CONDITION THAT WOULD BE REAL, not attempted here: g_4 is the FOUR-")
P("      dimensional SU(2)_L coupling and this closure runs it inside the Standard Model. Above")
P("      1/R_5 the theory is six-dimensional and the running is power law, so matching")
P("      g_4(1/R_5) to the SU(7) coupling of the bulk is an independent condition on the same")
P("      number -- and unlike this one it does not have five gentle e-folds to work in.")

(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "running_g4.json").write_text(json.dumps(dict(
    g2_MZ=g2_at(RG.MZ), g4_fixed=gstar, g4_paper=G4_PAPER, g4_ahmn=G4_AHMN,
    scale_of_063=s63, scale_of_0653=s653,
    rows=[dict(label=l, g4=g, A4=t, invR=r) for l, g, t, r in rows]), indent=1), encoding="utf-8")
P("")
P("archived: outputs/running_g4.json")
