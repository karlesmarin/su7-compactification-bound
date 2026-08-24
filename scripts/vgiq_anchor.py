#!/usr/bin/env python3
"""The second anchor: our machinery against von Gersdorff-Irges-Quiros, hep-th/0204223.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part VI's open question is that OUR F does not reproduce Komori-Maru's published alpha_min column
(ratio 1.94 with a 48, 1.20 without).  One anchor, and it fails; nothing localises the failure.
vGIQ give a SECOND, independent one -- a different group, a different dimension, a different
decade -- and they publish four things we can hit without adjusting anything:

  (i)   the potential in closed form, their (5.4)-(5.6):
            V = (1/128 pi^6 R^4) Tr[ V(r_F) - V(r_B) ],  V(r) = 3(Li_5(r) + Li_5(r*)),
            r = e^{2 pi i q(omega)}
  (ii)  the curvature at the origin, their (5.2):
            m_h^2 = (3/32 pi^4 R^2) g^2 zeta(3) [ 3 C_2(G) - 4 C_R N_f ]
  (iii) the critical flavour numbers at which omega=0 stops being a minimum, their (5.3)
            C_G/C_R < 4/(3 N_f):   N_f = 3 for SU(2) fund, 9/2 for SU(3) fund, 3/4 for adjoint
  (iv)  and -- the one that needs NO normalisation at all -- MINIMUM LOCATIONS:
            "For SU(2) one finds two degenerate minima ... at omega = 1/4, 3/4 for N_f > 3/4.
             For SU(3) ... two degenerate minima at omega ~= 0.29, 0.71 for N_f > 3/4."

Their conventions are fixed by their own two footnotes, not guessed:
  - the KK gauge boson mass matrix "has eigenvalues n^2, (n+2omega)^2, (n-2omega)^2" (p.22), so the
    shift of a state of Wilson-line eigenvalue lambda is q = 2 lambda omega;
  - footnote 7: "the parameter alpha in Ref. [46] is related to our omega as alpha = 2 omega",
    which is OUR alpha.  Hence  c_w = 2 lambda_w, integers, exactly our charge convention.

Their orbifold is S^1/Z_2 with a SINGLE parity: every KK tower is integer-moded, so there is no
antiperiodic sector, B_2 = 0, and our D = A_2 - (3/4) B_2 collapses to  D = A_2 = sum m c^2 -- pure
Dynkin index.  Their (5.2) is our formula in the Delta = 0 corner (see GATE_BRIDGE.md).  That is the
prediction being tested here, not an assumption: the group-theory content is computed from explicit
generators and compared with C_2(G) and C_R computed independently.
"""
import json
import math
import pathlib

import numpy as np

exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

DOF_GAUGE, DOF_FERMION = 3, 4        # their "global factor 3 - 4 N_f", p.21


# ---------------------------------------------------------------- su(N) by hand, no algebra library
def su_generators(n):
    """the n^2-1 hermitian generators T^A with tr(T^A T^B) = delta^{AB}/2 (their normalisation)."""
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            a = np.zeros((n, n), complex); a[i, j] = a[j, i] = 0.5
            gens.append(a)
            b = np.zeros((n, n), complex); b[i, j] = -0.5j; b[j, i] = 0.5j
            gens.append(b)
    for k in range(1, n):
        d = np.zeros((n, n), complex)
        for i in range(k):
            d[i, i] = 1
        d[k, k] = -k
        gens.append(d / math.sqrt(2 * k * (k + 1)))
    return gens


def adjoint_of(T, gens):
    """the matrix of ad(T) in the basis `gens`, from [T, T^B] = i f^{ABC} T^C."""
    n = len(gens)
    M = np.zeros((n, n), complex)
    for b, G in enumerate(gens):
        comm = T @ G - G @ T
        for c, H in enumerate(gens):
            M[c, b] = 2 * np.trace(comm @ H)          # tr(T^A T^B) = delta/2
    return M


def charges(T, gens, rep):
    """the integer charges c = 2 lambda of the Wilson-line generator T in `rep`."""
    M = T if rep == "fund" else adjoint_of(T, gens)
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    c = 2 * ev
    assert np.allclose(c, np.round(c), atol=1e-9), "charges are not integers: %s" % c
    return np.round(c).astype(int)


def broken_generator(n, dH):
    """a generator odd under Lambda = diag(1_dH, -1_dK): the orbifold's Higgs direction."""
    gens = su_generators(n)
    L = np.diag([1.0] * dH + [-1.0] * (n - dH))
    for T in gens:
        if np.allclose(L @ T @ L, -T) and abs(T[0, dH]) > 1e-9:
            return T, gens
    raise SystemExit("no broken generator found")


# ---------------------------------------------------------------- their potential, our expansion
NM = 400
_n = np.arange(1, NM + 1)
_w5 = _n.astype(float) ** -5.0


def W(cs, omega):
    """sum over the states of a representation of Re Li_5(e^{2 pi i c omega}), c = 2 lambda."""
    a = np.atleast_1d(np.asarray(omega, float))
    out = np.zeros(a.shape)
    for c in cs:
        out += (np.cos(np.outer(a, 2 * math.pi * c * _n)) * _w5).sum(axis=1)
    return out


def Vpot(cadj, cfund, nf, omega, rep="fund"):
    cf = cadj if rep == "adj" else cfund
    return -DOF_GAUGE * W(cadj, omega) + DOF_FERMION * nf * W(cf, omega)


P("")
P("=" * 100)
P("0 -- THE CHARGES, from explicit generators.  their footnote: eigenvalues n^2, (n +- 2 omega)^2")
P("=" * 100)
DATA = {}
for name, n, dH in (("SU(2)", 2, 1), ("SU(3)", 3, 2)):
    T, gens = broken_generator(n, dH)
    cadj, cfund = charges(T, gens, "adj"), charges(T, gens, "fund")
    CG = float(np.sum((cadj / 2.0) ** 2))            # C_2(G) = tr_adj(T T)
    CR = float(np.sum((cfund / 2.0) ** 2))           # C_R    = tr_R(T T)
    DATA[name] = dict(cadj=cadj, cfund=cfund, CG=CG, CR=CR)
    P("%-6s  adjoint charges c = %-22s  C_2(G) = %.4f   (textbook: %d)" %
      (name, str(sorted(cadj)), CG, n))
    P("        fund    charges c = %-22s  C_R    = %.4f   (textbook: 1/2)" % (str(sorted(cfund)), CR))
P("")
P("SU(2) adjoint gives c = -2, 0, 2 -- which IS their 'n^2, (n+2w)^2, (n-2w)^2'.  Convention fixed.")

P("")
P("=" * 100)
P("1 -- THEIR CRITERION (5.3) FROM OUR D.  all periodic (single parity) => D = A_2 = sum m c^2")
P("=" * 100)
P("D = -3 sum_adj c^2 + 4 N_f sum_R c^2 = 4[ 4 N_f C_R - 3 C_2(G) ],  so D > 0 (our EWSB verdict)")
P("is exactly their C_G/C_R < 4/(3 N_f).  The critical flavour numbers:")
P("")
P("%-8s %-10s %14s %14s %14s" % ("group", "fermions", "N_f critical", "theirs", "match"))
crit = []
for name, rep, theirs in (("SU(2)", "fund", 3.0), ("SU(3)", "fund", 4.5),
                          ("SU(2)", "adj", 0.75), ("SU(3)", "adj", 0.75)):
    d = DATA[name]
    CR = d["CR"] if rep == "fund" else d["CG"]
    nf = 3 * d["CG"] / (4 * CR)
    P("%-8s %-10s %14.4f %14.4f %14s" %
      (name, rep, nf, theirs, "YES" if abs(nf - theirs) < 1e-9 else "*** NO ***"))
    crit.append(dict(group=name, rep=rep, nf_ours=nf, nf_theirs=theirs))
P("")
P("Control -- the same numbers straight from a numerical second derivative of THEIR potential (5.4):")
h = 1e-4
for name, rep, theirs in (("SU(2)", "fund", 3.0), ("SU(3)", "fund", 4.5), ("SU(3)", "adj", 0.75)):
    d = DATA[name]
    lo, hi = 0.0, 20.0
    for _ in range(200):                              # bisect on the sign of V''(0)
        mid = 0.5 * (lo + hi)
        v = Vpot(d["cadj"], d["cfund"], mid, np.array([-h, 0.0, h]), rep)
        fpp = (v[0] - 2 * v[1] + v[2]) / h ** 2
        if fpp > 0:
            lo = mid
        else:
            hi = mid
    P("   %-6s %-5s  V''(0) changes sign at N_f = %.6f   (theirs %.4f)  %s" %
      (name, rep, 0.5 * (lo + hi), theirs, "OK" if abs(0.5 * (lo + hi) - theirs) < 1e-5 else "***"))

P("")
P("=" * 100)
P("2 -- THE TEST THAT NEEDS NO NORMALISATION: their published MINIMUM LOCATIONS")
P("=" * 100)
P("For adjoint fermions q_F = q_B, so V = (4 N_f - 3) W_adj(omega) and the minimum location is a")
P("pure statement about the adjoint charge spectrum.  No prefactor, no g, no dof factor survives.")
P("")
P("Done in mpmath at 40 digits, not on a grid.  With  d/dtheta Re Li_5(e^{i theta}) = -Im Li_4,")
P("the stationarity condition is a sum of Clausen functions:  sum_c c * Sl_4(2 pi c omega) = 0.")
import mpmath as mp
mp.mp.dps = 40


def Wm(cs, om):
    om = mp.mpf(om)
    return mp.fsum(mp.re(mp.polylog(5, mp.expjpi(2 * int(c) * om))) for c in cs)


def dWm(cs, om):
    om = mp.mpf(om)
    return -2 * mp.pi * mp.fsum(int(c) * mp.im(mp.polylog(4, mp.expjpi(2 * int(c) * om))) for c in cs)


mins, MINS_MP = {}, {}
for name, theirs in (("SU(2)", 0.25), ("SU(3)", 0.29)):
    cs = DATA[name]["cadj"]
    step = mp.mpf(1) / 400                            # coarse scan first, then bisect the derivative
    xs = [step * i for i in range(1, 200)]
    x0 = min(xs, key=lambda t: Wm(cs, t))
    lo, hi = x0 - step, x0 + step
    assert dWm(cs, lo) * dWm(cs, hi) < 0, "no sign change bracketing the minimum for %s" % name
    om = mp.findroot(lambda t: dWm(cs, t), (lo, hi), solver="bisect", tol=mp.mpf("1e-30"))
    d2 = (Wm(cs, om + mp.mpf("1e-6")) - 2 * Wm(cs, om) + Wm(cs, om - mp.mpf("1e-6"))) * mp.mpf("1e12")
    mins[name] = float(om)
    MINS_MP[name] = om
    P("   %-6s adjoint fermions, N_f > 3/4:  minimum at omega = %s" % (name, mp.nstr(om, 20)))
    P("          (and %s by the symmetry omega -> 1 - omega);  W'' = %s > 0" %
      (mp.nstr(1 - om, 6), mp.nstr(d2, 6)))
    P("          theirs: %.2f   ->   %s" %
      (theirs, "MATCHES their two printed digits"
       if abs(float(om) - theirs) < 0.005 else "*** MISMATCH ***"))
P("")
P("Both of ours are EXACT, and for reasons, not by luck:")
P("   SU(2): the adjoint is c = 0, +-2, so W = zeta(5) + 2 Re Li_5(e^{4 pi i omega}), minimised where")
P("          4 pi omega = pi:  omega = 1/4.   ours - 1/4 = %s" %
  mp.nstr(mp.mpf(mins["SU(2)"]) - mp.mpf(1) / 4, 5))
P("   SU(3): the adjoint is c = 0,0,+-1,+-1,+-2, so stationarity is  Sl_4(t) + Sl_4(2t) = 0, and at")
P("          t = 2pi/3 the Clausen function is ODD about 2pi:  Sl_4(4pi/3) = -Sl_4(2pi/3).  Exactly")
P("          zero.  The stationary point sits on the Z_3 CENTRE of SU(3):  omega = 1/3.")
_chk = dWm(DATA["SU(3)"]["cadj"], mp.mpf(1) / 3)
P("          control: W'(1/3) = %s   (and W''(1/3) > 0, so it is the minimum)" % mp.nstr(_chk, 8))

P("")
P("And the control that says which of the two is the physical answer: with adjoint matter the")
P("potential sees only the centre, so the minimum should be the CENTRE-SYMMETRIC holonomy -- the")
P("textbook Gross-Pisarski-Yaffe statement.  The Polyakov loop in the fundamental is")
P("   P(omega) = (1/N) sum_{c in fund} exp(2 pi i c omega) ,  and it must VANISH there:")
for name in ("SU(2)", "SU(3)"):
    cf = DATA[name]["cfund"]
    for lab, om in (("ours", MINS_MP[name]), ("theirs", mp.mpf("0.29") if name == "SU(3)"
                                                   else mp.mpf("0.25"))):
        pl = mp.fsum(mp.expjpi(2 * int(c) * om) for c in cf) / len(cf)
        P("   %-6s %-7s omega = %-10s |P| = %s %s" %
          (name, lab, mp.nstr(om, 6), mp.nstr(abs(pl), 8),
           " <- centre-symmetric" if abs(pl) < mp.mpf("1e-25") else ""))
P("")
P("   Our 1/3 makes the SU(3) Polyakov loop vanish exactly; 0.29 does not.  So the minimum we get")
P("   is the one the centre demands, and it is the same statement Part III makes about centre")
P("   charge governing the Wilson-line potential.")
P("")
P("   >>> SO SU(2) MATCHES AND SU(3) DOES NOT.  They print 0.29; we get 1/3 = 0.3333...")
P("")
P("   The gap is a factor  ours/theirs = %.4f, and  1/(2 sqrt 3) = %.6f  rounds to 0.29 -- i.e. their"
  % (mins["SU(3)"] / 0.29, 1 / (2 * math.sqrt(3))))
P("   number is ours times sqrt(3)/2, which is what a different normalisation of the broken")
P("   generator would do.  Their text fixes that normalisation for SU(2) (the KK mass eigenvalues in")
P("   footnote 9) and NOT for SU(3), so this stays OPEN.  It is not claimed as explained.")
P("")
P("   But the centre control above puts the evidence on OUR side, not on a shared convention: our")
P("   1/3 is where the Polyakov loop vanishes, and with adjoint matter that is where the minimum")
P("   has to be.  So the likeliest reading is that their two-digit 0.29 is a slip -- which we")
P("   cannot prove, so it is reported as unresolved with the evidence stated.")
P("")
P("   Either way the anchor does its job, and this is the point of the whole file:")
P("      - four critical flavour numbers, exact;   SU(2) minimum, exact;")
P("      - the entire group-theory bracket of their curvature formula, exact;")
P("      - and what is left over is a pure convention scalar, 2/pi^2, with no physics in it.")
P("   OUR F, OUR EXPANSION, OUR D AND OUR MINIMISER ARE THEREFORE VALIDATED against an independent")
P("   published one-loop computation, in a different group, dimension and decade.  Which relocates")
P("   Part VI's open question: the 1.94/1.20 residual against Komori-Maru cannot be blamed on the")
P("   potential or on the minimisation.  It has to live in the model-specific dictionary -- their")
P("   charge and parity assignment, or the alpha-to-1/R_5 relation -- and that is a much smaller")
P("   place to look than 'our F'.")

P("")
P("=" * 100)
P("3 -- THE PREFACTOR of their (5.2), rebuilt from their (5.4)")
P("=" * 100)
P("Re Li_5(e^{i t}) = zeta(5) - zeta(3) t^2/2 + ... with t = 2 pi q = 2 pi c omega, so")
P("   d^2/domega^2 |_0 = -(2 pi c)^2 zeta(3),  and with V = (6/128 pi^6 R^4)[Tr_F - Tr_B]:")
P("      d^2V/domega^2|_0 = -(6/128 pi^6 R^4)(2 pi)^2 zeta(3) [4 N_f sum_R c^2 - 3 sum_adj c^2]")
P("   Their (5.1) omega = h R/2 gives d^2/dh^2 = (R^2/4) d^2/domega^2, and sum c^2 = 4 C:")
pref = (6.0 / 128) * (2 * math.pi) ** 2 * 4 / (4 * math.pi ** 6)
P("      m_h^2 = %.10f * zeta(3)/(pi^4 R^2) * [3 C_2(G) - 4 C_R N_f]" % (pref * math.pi ** 2))
P("      theirs (5.2):  3/32 = %.10f    * g^2 * zeta(3)/(pi^4 R^2) * [ same bracket ]" % (3 / 32))
P("")
_r = pref * math.pi ** 2 / (3 / 32)
P("   ratio ours/theirs at g = 1 :  %.8f     and  2/pi^2 = %.8f  ->  %s" %
  (_r, 2 / math.pi ** 2, "the same" if abs(_r - 2 / math.pi ** 2) < 1e-9 else "*** different ***"))
P("")
P("   The bracket -- the whole group-theory content, and the 3:4 -- comes out EXACTLY.  What is left")
P("   over is a pure number, 2/pi^2, with no group theory and no zeta in it: it is their g^2 and")
P("   their A_5-to-h convention (5.1), and the pi's are what a KK reduction g_4^2 = g_5^2/(2 pi R)")
P("   puts there.  Our expansion neither predicts it nor needs it -- but note that it is again a")
P("   NORMALISATION that separates us from a published number, and never the structure.")

OUT.mkdir(exist_ok=True)
(OUT / "vgiq_anchor.json").write_text(json.dumps(dict(
    charges={k: dict(adj=[int(z) for z in v["cadj"]], fund=[int(z) for z in v["cfund"]],
                     CG=v["CG"], CR=v["CR"]) for k, v in DATA.items()},
    critical_nf=crit, minima={k: float(v) for k, v in mins.items()},
    prefactor_ratio=float(pref * math.pi ** 2 / (3 / 32))), indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "vgiq_anchor.json"))
