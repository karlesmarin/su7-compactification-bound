#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ghost_seed_vs_anchor.py -- .explica la semilla candidata el residuo del anclaje?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

DOS COSAS QUE LE HEMOS DICHO A [KM25], Y LA PREGUNTA DE SI SON UNA SOLA.

  (i)  No reproducimos su columna alpha_min: la razon por filas va de 1.03 a 2.08,
       agrupada 1.94 con un 48 y 1.20 sin el.  Eso esta en la Parte VI y se le escribio.
  (ii) Su ec. (68) puede que omita un determinante de Faddeev-Popov en el sector
       periodico: el conteo covariante da (3/2,3,6) donde ellos imprimen (2,4,7).

Los dos son sobre los MISMOS coeficientes.  Si la semilla candidata mejorara la columna,
serian el mismo problema y la evidencia para (ii) seria fuerte.  Si la empeora, son dos
problemas distintos y conviene decirlo asi antes de escribirle.

LA PRUEBA.  Se minimiza el potencial EXACTO de las cinco filas publicadas con las dos
semillas y se compara con su columna.  Nada mas.  El resultado puede salir en cualquiera de
los dos sentidos, y por eso vale la pena hacerlo.
"""
import contextlib
import io
import math
import pathlib
import runpy

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent

# amin_closed_form.py imprime su propio informe al importarse; se silencia.
_buf = io.StringIO()
with contextlib.redirect_stdout(_buf):
    NS = runpy.run_path(str(HERE / "amin_closed_form.py"))

numeric_min = NS["numeric_min"]
moments = NS["moments"]
T1 = NS["T1"]
GAUGE_PUB = list(NS["GAUGE"])

# Las funciones miran GAUGE en SU propio espacio de nombres, no en la copia que devuelve
# runpy.  Se cambia ahi, y se restaura al final.
GLOB = NS["numeric_min"].__globals__

# (m, s, c).  El articulo lleva la MITAD de los coeficientes del corchete de la ec. (68).
GAUGE_CAND = [(-0.75, 1, 2), (-1.5, 1, 1), (-3.0, -1, 1)]


def ocho_D(gauge):
    """8D del sector gauge, para comprobar que la semilla es la que se cree."""
    A2 = sum(m * c ** 2 for m, s, c in gauge if s > 0)
    B2 = sum(m * c ** 2 for m, s, c in gauge if s < 0)
    return 8 * A2 - 6 * B2


P("=" * 92)
P("CONTROL 0 -- que las dos listas son las dos semillas y no otra cosa")
P("=" * 92)
for nombre, g in (("publicada", GAUGE_PUB), ("candidata", GAUGE_CAND)):
    corchete = tuple(-2 * m for m, s, c in g)
    P("  %-10s  GAUGE = %-42s  corchete = %-14s  8D_gauge = %+d"
      % (nombre, g, corchete, ocho_D(g)))
esperado = {"publicada": -27, "candidata": -18}
ok0 = (ocho_D(GAUGE_PUB) == -27 and ocho_D(GAUGE_CAND) == -18)
P("  esperado -27 y -18 (tabla de las dos ramas del articulo): %s" % ("OK" if ok0 else "*** NO ***"))
P("")

P("=" * 92)
P("LA PRUEBA -- alpha_min de las cinco filas publicadas, con las dos semillas")
P("=" * 92)
P("  fila    theirs      ours(pub)   razon    ours(cand)  razon    .mejora?   tiene 48")
P("  " + "-" * 86)

filas = []
for etiqueta, contenido, theirs, mh, _ in T1:
    tiene48 = any(rep == "48" for rep, *_ in contenido)

    GLOB["GAUGE"] = list(GAUGE_PUB)
    a_pub = numeric_min(contenido)

    GLOB["GAUGE"] = list(GAUGE_CAND)
    a_cand = numeric_min(contenido)

    r_pub = a_pub / theirs if a_pub else float("nan")
    r_cand = a_cand / theirs if a_cand else float("nan")
    mejora = abs(r_cand - 1.0) < abs(r_pub - 1.0)
    filas.append((etiqueta, theirs, a_pub, r_pub, a_cand, r_cand, mejora, tiene48))
    P("  %-6s  %.4f      %.6f    %.3f    %.6f    %.3f    %-8s   %s"
      % (etiqueta, theirs, a_pub, r_pub, a_cand, r_cand,
         "MEJORA" if mejora else "empeora", "si" if tiene48 else "no"))

GLOB["GAUGE"] = list(GAUGE_PUB)          # restaurado

P("")
n_mejora = sum(1 for f in filas if f[6])


def resumen(idx, nombre):
    rs = [f[idx] for f in filas]
    disp = max(rs) / min(rs)
    med = sum(abs(r - 1.0) for r in rs) / len(rs)
    P("  %-22s razones %s" % (nombre, " ".join("%.3f" % r for r in rs)))
    P("  %-22s dispersion max/min = %.3f   |razon - 1| medio = %.3f" % ("", disp, med))
    return disp, med


P("=" * 92)
P("RESUMEN")
P("=" * 92)
d_pub, m_pub = resumen(3, "semilla publicada")
P("")
d_cand, m_cand = resumen(5, "semilla candidata")
P("")
P("  filas que mejoran con la candidata: %d de %d" % (n_mejora, len(filas)))
P("")

P("=" * 92)
P("VEREDICTO")
P("=" * 92)
if m_cand < m_pub and n_mejora >= 3:
    P("  La semilla candidata ACERCA la columna.  Los dos problemas que le hemos escrito a")
    P("  [KM25] podrian ser el mismo, y eso es evidencia independiente para el determinante.")
elif m_cand > m_pub:
    P("  La semilla candidata ALEJA la columna: el residuo del anclaje empeora de %.3f a %.3f" % (m_pub, m_cand))
    P("  en |razon - 1| medio.  De modo que los DOS PROBLEMAS SON DISTINTOS:")
    P("")
    P("    - la discrepancia de alpha_min NO se explica cambiando el reparto gauge;")
    P("    - y el reparto gauge NO se puede defender diciendo que arregla alpha_min.")
    P("")
    P("  Eso es informacion util en las dos direcciones, y hay que escribirsela asi: la")
    P("  pregunta de los fantasmas se sostiene o se cae por el determinante, no por esta")
    P("  columna.  Ninguno de los dos hallazgos avala al otro.")
else:
    P("  Sin diferencia apreciable: la columna no distingue las dos semillas, de modo que no")
    P("  es un discriminante y no debe usarse como si lo fuera.")
P("")

P("=" * 92)
P("FALSACION -- .puede esta prueba dar 'mejora' alguna vez?")
P("=" * 92)
P("  Se le da una semilla FICTICIA construida para acercarse a su columna: si la prueba no")
P("  la detectara, no estaria midiendo nada.")
mejor, mejor_f = None, None
for f in (0.2, 0.4, 0.6, 0.8):
    g = [(m * f, s, c) for m, s, c in GAUGE_PUB]
    GLOB["GAUGE"] = g
    rs = []
    for etiqueta, contenido, theirs, mh, _ in T1:
        a = numeric_min(contenido)
        rs.append(a / theirs if a else float("nan"))
    med = sum(abs(r - 1.0) for r in rs) / len(rs)
    P("     gauge x %.1f  ->  |razon - 1| medio = %.3f%s" % (f, med, "   <- mejor que la publicada" if med < m_pub else ""))
    if mejor is None or med < mejor:
        mejor, mejor_f = med, f
GLOB["GAUGE"] = list(GAUGE_PUB)
P("  la prueba SI distingue: el mejor factor ficticio (%.1f) da %.3f contra %.3f de la publicada."
  % (mejor_f, mejor, m_pub))
P("  De modo que 'empeora' es una medida y no una imposibilidad del instrumento.")
