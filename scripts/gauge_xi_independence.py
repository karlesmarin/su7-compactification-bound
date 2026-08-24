#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gauge_xi_independence.py -- .depende de xi el reparto gauge por paridad?  No.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

La CONSULTA_2026-08-23_fantasmas_gauge.md hace tres preguntas.  Esta contesta la tercera:
si el coeficiente del sector gauge puede depender del parametro de gauge xi.  Las otras dos
---las condiciones de contorno de los fantasmas, y si [KM25] tienen una cancelacion
implicita--- siguen abiertas, y el veredicto del final lo dice.

METODO: se enuncia la afirmacion y se intenta ROMPER, no confirmar.  Cinco objeciones, cada
una contestada o aceptada.  La quinta ---la suma sobre modos del termino en ln xi--- es la
que por poco se escapa, y se contesta por dos caminos independientes.

Complementa a single_orbit_determinant.py, que CUENTA modos cero; esto calcula el operador.
"""
import numpy as np
from fractions import Fraction as F

np.random.seed(20260824)
D = 6
I_A6 = 5                       # indice de A_6
rng = np.random.default_rng(20260824)

def delta(p, xi):
    return np.eye(D) * (p @ p) - (1 - 1 / xi) * np.outer(p, p)

def linea(t):
    print(t)

linea("=" * 78)
linea("AFIRMACION.  El reparto por paridad del sector gauge es 3 (sector eps) y 1 (sector")
linea("             -eps), y no depende de xi en lo unico que le importa al potencial: su")
linea("             dependencia en alpha.")
linea("=" * 78)
linea("")

# ---------------------------------------------------------------- objecion 1
linea("OBJECION 1.  'El determinante no es tan simple; hay que calcularlo de verdad.'")
linea("-" * 78)
malos = 0
for _ in range(200):
    p = rng.normal(size=D)
    xi = float(np.exp(rng.normal() * 2))          # xi en varios ordenes de magnitud
    ev = np.sort(np.linalg.eigvalsh(delta(p, xi)))
    p2 = p @ p
    esperado = np.sort(np.array([p2 / xi] + [p2] * (D - 1)))
    if not np.allclose(ev, esperado, rtol=1e-10):
        malos += 1
linea("  200 momentos y xi aleatorios, autovalores contra {p^2/xi, p^2 x5}: %d fallos" % malos)
linea("  RESPUESta: Delta es p^2 I menos una actualizacion de RANGO UNO.  Un autovalor")
linea("  longitudinal p^2/xi y D-1 = 5 transversos p^2.  det = (p^2)^D / xi.")
linea("  La objecion no se sostiene: el determinante es exacto y cerrado.")
linea("")

# ---------------------------------------------------------------- objecion 2
linea("OBJECION 2.  'El termino de xi MEZCLA los dos sectores de paridad, asi que no hay")
linea("             reparto que valga.'  --- Y esta objecion es CORRECTA en general.")
linea("-" * 78)
p = rng.normal(size=D)
xi = 3.7
Dl = delta(p, xi)
linea("  La entrada que conectaria A_6 con A_0 vale  Delta_{60} = -(1-1/xi) p_6 p_0")
linea("     numericamente: %.6f     (p_6 = %.4f)" % (Dl[I_A6, 0], p[I_A6]))
linea("  No es cero.  Con xi distinto de 1 el operador NO se parte en bloques de paridad")
linea("  definida.  La objecion se sostiene, y hay que contestarla, no esquivarla.")
linea("")

# ---------------------------------------------------------------- objecion 3
linea("OBJECION 3.  'Vale, pero .que modos llevan la dependencia en alpha?'")
linea("-" * 78)
linea("  alpha entra por la torre de x_5.  La direccion x_6 da p_6 = m/R_6.  C1 de")
linea("  single_orbit_determinant.py MIDE que con R_5/R_6 = 1000 la parte alpha-dependiente")
linea("  de los modos m >= 1 pesa menos de 1e-100.  De modo que la dependencia en alpha la")
linea("  llevan los modos con p_6 = 0, y solo ellos.")
linea("")
linea("  *** Y AQUI HAY QUE SER HONESTO: eso NO dice que los modos m >= 1 no contribuyan.")
linea("  Contribuyen una constante enorme, independiente de alpha.  La afirmacion es sobre")
linea("  la DEPENDENCIA EN ALPHA, que es lo que es el potencial.  Dicho asi y no mas.")
linea("")
p0 = p.copy(); p0[I_A6] = 0.0
D0 = delta(p0, xi)
linea("  Con p_6 = 0:  Delta_{6a} para a != 6  ->  %s" % np.array2string(
    np.array([D0[I_A6, a] for a in range(D) if a != I_A6]), precision=1))
linea("                Delta_{66}            ->  %.6f   contra p^2 = %.6f" % (D0[I_A6, I_A6], p0 @ p0))
linea("  El termino de xi lleva un factor p_6 en toda entrada que toque el indice 6.")
linea("  Con p_6 = 0 el operador SE PARTE, y A_6 queda con p^2 limpio, sin xi.")
linea("")

# ---------------------------------------------------------------- objecion 4
linea("OBJECION 4.  'Aun asi A_6 y A_mu estan los dos ahi y podrian mezclarse.'")
linea("-" * 78)
linea("  No lo estan, y esta es la parte estructural.  Un generador de paridad P_6 = eps:")
linea("     eps = +1 :  A_mu, A_5 ~ cos(m x_6/R_6)  ->  TIENEN modo m = 0")
linea("                 A_6       ~ sin(m x_6/R_6)  ->  NO tiene modo m = 0")
linea("     eps = -1 :  al reves.")
linea("  De modo que en p_6 = 0 solo esta presente UNA de las dos clases de paridad.")
linea("  No hay nada que mezclar, ni siquiera antes de mirar el factor p_6.")
linea("  Dos razones independientes para lo mismo, que es como conviene que sea.")
linea("")

# ---------------------------------------------------------------- el calculo
linea("=" * 78)
linea("EL CALCULO, SECTOR A SECTOR, CON xi GENERAL Y p_6 = 0")
linea("=" * 78)

def sector(nombre, presentes, fantasmas, xi, p0):
    """Gamma del sector, en unidades de (1/2) ln p^2, y el coeficiente de ln xi."""
    sub = np.array(presentes)
    M = delta(p0, xi)[np.ix_(sub, sub)]
    s, logdet = np.linalg.slogdet(M)
    p2 = p0 @ p0
    # ln det = k ln p^2 - c ln xi  ->  se despeja c comparando con xi = 1
    M1 = delta(p0, 1.0)[np.ix_(sub, sub)]
    _, logdet1 = np.linalg.slogdet(M1)
    c = -(logdet - logdet1) / np.log(xi)
    k = logdet1 / np.log(p2)
    gamma_lnp2 = 0.5 * k - fantasmas * 1.0        # cada par de fantasmas resta ln p^2
    linea("  %s" % nombre)
    linea("     componentes presentes : %d      pares de fantasmas : %d" % (len(sub), fantasmas))
    linea("     ln det = %.4f ln p^2  -  %.4f ln xi" % (k, c))
    linea("     Gamma  = %.4f ln p^2  -  %.4f ln xi" % (gamma_lnp2, 0.5 * c))
    linea("     grados de libertad (un escalar real = 1/2 ln p^2) : %.4f" % (2 * gamma_lnp2))
    return 2 * gamma_lnp2, 0.5 * c

xi = 7.3
g_mas, c_mas = sector("eps = +1  ->  A_mu (4) + A_5 (1), y los dos fantasmas TIENEN modo cero",
                      [0, 1, 2, 3, 4], 1, xi, p0)
linea("")
g_menos, c_menos = sector("eps = -1  ->  solo A_6; los fantasmas heredan eps y NO tienen modo cero",
                          [I_A6], 0, xi, p0)
linea("")
linea("  TOTAL : %.4f + %.4f = %.4f  grados  =  D-2 en seis dimensiones" % (g_mas, g_menos, g_mas + g_menos))
linea("")

# ---------------------------------------------------------------- objecion 5
linea("=" * 78)
linea("OBJECION 5, LA SERIA.  'Queda un -1/2 ln xi en el sector eps=+1.  .Y si la SUMA")
linea("             SOBRE MODOS lo vuelve alpha-dependiente?'")
linea("=" * 78)
linea("  Esta objecion es la buena, y por poco se me escapa.  El termino es")
linea("     -(1/2) ln xi  x  (numero de modos)")
linea("  y en regularizacion zeta el 'numero de modos' de una torre DESPLAZADA es")
linea("     sum_{n in Z} 1  =  zeta_H(0, a) + zeta_H(0, 1-a)")
linea("  con a el desplazamiento que mete la linea de Wilson.  Y zeta_H(0,a) = 1/2 - a")
linea("  DEPENDE de a.  Asi que hay que hacer la cuenta, no suponerla:")
for a in (F(0), F(1, 4), F(1, 3), F(1, 2), F(7, 10)):
    tot = (F(1, 2) - a) + (F(1, 2) - (1 - a))
    linea("     a = %-5s :  zeta_H(0,a) + zeta_H(0,1-a) = (%s) + (%s) = %s"
          % (a, F(1, 2) - a, F(1, 2) - (1 - a), tot))
linea("")
linea("  Se cancela IDENTICAMENTE, para todo a.  Las dos mitades de la torre aportan")
linea("  desplazamientos opuestos y su suma no ve el desplazamiento.")
linea("")
linea("  Y hay una segunda via, mas robusta porque no depende del regulador elegido: el")
linea("  coeficiente de ln xi no lleva p, de modo que su integral sobre el momento")
linea("  cuadridimensional es una integral SIN ESCALA, y en regularizacion dimensional")
linea("  una integral sin escala es cero.")
linea("")
linea("  La objecion se contesta por dos caminos independientes.  Pero merecia hacerse:")
linea("  'no lleva p, luego es constante' habria sido un paso en falso.")
linea("")

# ---------------------------------------------------------------- falsacion
linea("=" * 78)
linea("FALSACION.  Si la conclusion no dependiera de sus hipotesis, no probaria nada.")
linea("=" * 78)

pn = p.copy(); pn[I_A6] = 0.9
Dn = delta(pn, 4.0)
linea("  (a) con p_6 != 0 el bloque de paridad definida deja de existir:")
linea("      Delta_{60} = %.4f  ->  %s" % (Dn[I_A6, 0], "NO se anula, correcto" if abs(Dn[I_A6, 0]) > 1e-9 else "*** se anula: el control no puede fallar ***"))

linea("  (b) DECOY: si los fantasmas heredaran la paridad de A_6 en vez de la de A_mu,")
linea("      el reparto seria otro.  Se calcula, y tiene que salir DISTINTO:")
g2, _ = sector("      decoy  eps=+1 sin fantasmas", [0, 1, 2, 3, 4], 0, xi, p0)
linea("      decoy: %.1f y 1 -> total %.1f, que NO es D-2 = 4." % (g2, g2 + 1))
linea("      La hipotesis sobre los fantasmas ESTA haciendo trabajo.  Bien.")
linea("")

linea("  (c) en CINCO dimensiones la misma cuenta tiene que dar D-2 = 3:")
D5 = 5
def delta5(p, xi):
    return np.eye(D5) * (p @ p) - (1 - 1 / xi) * np.outer(p, p)
p5 = rng.normal(size=D5)
_, ld5 = np.linalg.slogdet(delta5(p5, 1.0))
k5 = ld5 / np.log(p5 @ p5)
g5 = 2 * (0.5 * k5 - 1)
linea("      D = 5:  (1/2)(%d) - 1 = %.1f  ->  %.1f grados = D-2" % (k5, g5 / 2, g5))
linea("")

# ---------------------------------------------------------------- veredicto
linea("=" * 78)
linea("VEREDICTO, Y LO QUE NO CIERRA")
linea("=" * 78)
linea("  CIERRA: el reparto es 3 y 1, suma D-2 = 4, y no depende de xi ---ni por mezcla de")
linea("  sectores, porque en p_6 = 0 solo hay una clase de paridad presente y ademas el")
linea("  termino de xi lleva p_6; ni por la suma sobre modos, porque el resto se cancela")
linea("  identicamente en zeta y es una integral sin escala en dim reg.")
linea("")
linea("  Eso es (a, b) = (3/2, 1/2): la semilla CANDIDATA.")
linea("")
linea("  NO CIERRA, y son las otras dos preguntas de la consulta:")
linea("   (a) que los fantasmas hereden la paridad de A_mu se ASUME aqui.  Es lo estandar")
linea("       ---el parametro gauge no tiene indice vectorial, luego transforma como A_mu---")
linea("       y es lo que hacen [HNT], pero asumirlo no es demostrarlo, y el decoy (b) de")
linea("       arriba muestra que la conclusion depende de ello.")
linea("   (b) si [KM25] tienen una cancelacion implicita en como construyen su espectro que")
linea("       reponga el grado que falta.  Eso solo lo puede decir quien escribio ese")
linea("       espectro.")
linea("")
linea("  Y una tercera, que es MIA y no estaba en la consulta:")
linea("   (c) todo esto es a orden dominante en R_5/R_6, que es lo que C1 mide.  La")
linea("       correccion de los modos m >= 1 a la DEPENDENCIA EN ALPHA es < 1e-100 a")
linea("       R_5/R_6 = 1000, pero el articulo no fija ese cociente: [KM25] solo piden")
linea("       R_5 >> R_6.  A cociente moderado la cuenta habria que rehacerla.")
linea("")
linea("  De modo que esto NO convierte el Teorema 1 en incondicional.  Reduce el hueco de")
linea("  tres preguntas a dos, y quita de en medio la unica que era puramente tecnica.")
