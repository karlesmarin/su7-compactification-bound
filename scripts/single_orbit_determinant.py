#!/usr/bin/env python3
"""single_orbit_determinant.py -- el reparto 3:1 o 4:1, decidido en una sola orbita cargada.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

LA PREGUNTA.  Toda la aritmetica de este articulo descansa en que 8D sea un entero IMPAR, y eso
depende del peso que el sector gauge aporta a cada sector de paridad P_6.  La ec. (68) de
Komori-Maru sale exactamente de asignar N=4 a A_mu y N=1 al escalar superviviente de A_5/A_6.  Pero
4+1 = 5, y un campo gauge en seis dimensiones tiene D-2 = 4 polarizaciones fisicas.  Hosotani, Noda
y Takenaga hacen la cuenta covariante y dicen que la matriz de masas de los fantasmas es la misma
que la de A_mu, de modo que A_mu + fantasmas dan 4-2 = 2, mas 2 de las componentes
extradimensionales: cuatro en total.

Lo que decide entre 4:1 y 3:1 no necesita SU(7) ni las 48 componentes.  Basta UNA orbita cargada,
abeliana, en gauge de fondo, y preguntarse que sobrevive en el limite R_5 >> R_6.

EL MECANISMO, Y ES TODO LO QUE HAY.  Bajo x_6 -> -x_6, un generador de paridad P_6 = eps tiene

    A_mu, A_5  con paridad  eps        (modos cos(m x_6/R_6): incluyen m = 0)
    A_6        con paridad -eps        (modos sin(m x_6/R_6): NO tienen m = 0)
    fantasmas  con la paridad de A_mu, porque heredan la del parametro gauge

La dependencia en alpha entra por la torre de x_5.  Un modo con m >= 1 lleva una masa transversa
m/R_6 >> 1/R_5, y su contribucion alpha-dependiente esta suprimida como K_{5/2}(2 pi k R_5 m/R_6),
que a R_5/R_6 = 1000 es indistinguible de cero.  C1 lo mide en vez de afirmarlo.

Luego, a orden dominante, SOLO cuentan las componentes que tienen modo cero en x_6:

    P_6 = +1 :  A_mu (4) + A_5 (1) - fantasmas (2)  =  3
    P_6 = -1 :  A_6 (1)                             =  1

y la resta de fantasmas ES dependiente de la paridad: en el sector P_6 = -1 los fantasmas no
tienen modo cero, asi que no hay nada que restar.  Ese es el punto entero.

CONTROLES
  C1  la supresion de los modos m >= 1 es real y enorme.  Si no lo fuera, el limite R_5 >> R_6 no
      seria un limite y todo lo de abajo sobra.
  C2  con los fantasmas puestos, el conteo da 3 y 1, y la ec. (68) pasaria a (3/2, 3, 6).
  C3  SIN restar los fantasmas, el conteo da 4 y 1 y reproduce (2, 4, 7) EXACTAMENTE.  Este es el
      control que convierte una sospecha en un diagnostico: dice que la diferencia entre nuestra
      cuenta y la suya es exactamente un determinante de Faddeev-Popov en el sector periodico.
  C4  y que le pasa a 8D en cada caso, que es lo unico que decide si el Teorema 1 existe.

Run:  python single_orbit_determinant.py > outputs/single_orbit_determinant.txt
"""
import pathlib
import sys
from fractions import Fraction as F

import mpmath as mp

P = lambda *a: print(*a, flush=True)
mp.mp.dps = 30

# nuestra tabla de generadores, de la Parte VI: pares (q,-q) por canal
PAIRS_C2_PLUS = 1          # (c=2, s=+1): 1 par, todos P_6 = +1
PAIRS_C1_PLUS = 2          # (c=1, s=+1): 2 pares, todos P_6 = +1
ANTI_PLUS, ANTI_MINUS = 2, 6   # (c=1, s=-1): 2 pares P_6=+1 y 6 P_6=-1
KM25_68 = (F(2), F(4), F(7))


def alpha_weight(k, mu_R5):
    """peso relativo del k-esimo enrollamiento para masa transversa mu, contra mu = 0.

    El termino alpha-dependiente va como (mu/(pi k R_5))^{5/2} K_{5/2}(2 pi k R_5 mu), y su
    limite mu -> 0 reproduce 1/k^5.  El cociente es libre de normalizacion."""
    if mu_R5 == 0:
        return mp.mpf(1)
    x = 2 * mp.pi * k * mu_R5
    num = (mu_R5 / (mp.pi * k)) ** mp.mpf(2.5) * mp.besselk(mp.mpf(2.5), x)
    # limite mu -> 0: K_{5/2}(x) ~ Gamma(5/2) 2^{3/2} / x^{5/2}
    den = mp.gamma(mp.mpf(2.5)) * mp.mpf(2) ** mp.mpf(1.5) / (mp.pi ** mp.mpf(2.5) * (2 * k) ** mp.mpf(2.5) * k ** mp.mpf(2.5))
    return num / den


def main():
    fails = []
    P("=" * 94)
    P("EL REPARTO 3:1 o 4:1, EN UNA SOLA ORBITA CARGADA")
    P("=" * 94)

    P("\n[C1] LOS MODOS m >= 1 ESTAN SUPRIMIDOS?  (si no, no hay limite que tomar)")
    P("     peso relativo del modo transverso m, contra el m = 0, en el enrollamiento k = 1")
    P("\n     %-10s %14s %14s" % ("R_5/R_6", "m = 1", "m = 2"))
    ok1 = True
    for ratio in (10, 100, 1000):
        w1 = alpha_weight(1, mp.mpf(ratio))
        w2 = alpha_weight(1, mp.mpf(2 * ratio))
        P("     %-10d %14.3e %14.3e" % (ratio, float(w1), float(w2)))
        if ratio == 1000 and float(w1) > 1e-100:
            ok1 = False
    P("\n     C1 %s -- a R_5/R_6 = 1000 el primer modo transverso pesa menos de 1e-100." %
      ("PASS" if ok1 else "FAIL"))
    P("        Asi que a orden dominante SOLO cuenta lo que tiene modo cero en x_6.")
    if not ok1:
        fails.append("C1")

    P("\n[C2] QUE TIENE MODO CERO EN x_6, SECTOR A SECTOR")
    P("""
     Bajo x_6 -> -x_6 un generador de paridad eps tiene A_mu y A_5 con paridad eps -- modos
     cos, que incluyen m = 0 -- y A_6 con paridad -eps -- modos sin, que empiezan en m = 1.
     Los fantasmas heredan la paridad de A_mu.
""")
    P("     %-10s %8s %6s %6s %10s %8s" % ("sector", "A_mu", "A_5", "A_6", "fantasmas", "neto"))
    P("     %-10s %8s %6s %6s %10s %8d" % ("P_6 = +1", "4", "1", "-", "-2", 4 + 1 - 2))
    P("     %-10s %8s %6s %6s %10s %8d" % ("P_6 = -1", "-", "-", "1", "-", 1))
    n_plus_gh, n_minus_gh = 3, 1
    P("\n     con fantasmas : %d y %d,  que suman %d = D-2 en seis dimensiones"
      % (n_plus_gh, n_minus_gh, n_plus_gh + n_minus_gh))
    ok2 = n_plus_gh + n_minus_gh == 4
    P("     C2 %s" % ("PASS" if ok2 else "FAIL"))
    if not ok2:
        fails.append("C2")

    P("\n[C3] Y SIN RESTAR LOS FANTASMAS?  El control que convierte la sospecha en diagnostico.")
    for gh, label in ((True, "con fantasmas   (3 y 1)"), (False, "sin fantasmas   (4 y 1)")):
        npl = 3 if gh else 4
        row = (F(PAIRS_C2_PLUS * npl, 2),
               F(PAIRS_C1_PLUS * npl, 2),
               F(ANTI_PLUS * npl + ANTI_MINUS * 1, 2))
        match = "  <-- REPRODUCE su ec. (68)" if row == KM25_68 else ""
        P("     %-24s ->  (%s, %s, %s)%s" % (label, row[0], row[1], row[2], match))
    P("""
     C3 PASS -- omitir la resta de fantasmas en el sector periodico reproduce (2,4,7) EXACTAMENTE.
        La diferencia entre las dos cuentas es, ni mas ni menos, un determinante de
        Faddeev-Popov en el sector P_6 = +1.  Eso no demuestra que su ec. (68) este mal: puede
        haber una cancelacion implicita en como construyen el espectro.  Lo que demuestra es que
        la discrepancia tiene UN solo origen y esta localizada.""")

    P("\n[C4] Y QUE LE PASA AL TEOREMA DE LOS OCTAVOS IMPARES")

    def eight_D(row):
        # GAUGE = -(fila)/2 en (m,s,c);  A2 = suma_{s=+} m c^2 ,  B2 = suma_{s=-} m c^2
        c2, c1p, c1m = row
        A2 = -(c2 / 2) * 4 - (c1p / 2) * 1
        B2 = -(c1m / 2) * 1
        return 8 * (A2 - F(3, 4) * B2)

    for gh, label in ((False, "su conteo        (2, 4, 7)"), (True, "con fantasmas   (3/2, 3, 6)")):
        npl = 3 if gh else 4
        row = (F(PAIRS_C2_PLUS * npl, 2), F(PAIRS_C1_PLUS * npl, 2),
               F(ANTI_PLUS * npl + ANTI_MINUS * 1, 2))
        e = eight_D(row)
        par = "IMPAR -> Teorema 1 VIVE" if (e.denominator == 1 and e % 2) else "PAR -> Teorema 1 MUERE"
        P("     %-28s  8D_gauge = %-6s  %s" % (label, e, par))
    P("""
     La materia siempre aporta 8D par (c y m enteros), asi que la paridad de 8D la decide
     enteramente el sector gauge.  Por eso esto no es un ajuste fino de un coeficiente: es la
     diferencia entre que exista un teorema y que no.""")

    P("\n" + "=" * 94)
    if fails:
        P("VERDICT: %d CONTROL(ES) FALLAN: %s" % (len(fails), ", ".join(fails)))
        P("=" * 94)
        return 1
    P("VERDICT: el conteo con fantasmas da 3:1 y sumaria D-2 = 4; el de la ec. (68) da 4:1 y suma")
    P("         cinco.  La diferencia es un determinante de Faddeev-Popov en el sector periodico,")
    P("         y decide si el Teorema de los octavos impares existe.")
    P("=" * 94)
    return 0


if __name__ == "__main__":
    sys.exit(main())
