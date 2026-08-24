#!/usr/bin/env python3
"""seed_shift_character.py -- las dos semillas gauge se diferencian en UN SOLO VECTOR, y cada
fila de la tabla de dos ramas de la S13 es una contraccion contra el.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

POR QUE EXISTE.  El articulo demuestra rama a rama que unos resultados sobreviven al cambio de
semilla gauge y otros no: la ley mod 6 sobrevive, `2W` impar sobrevive, `8D` impar no, la cota
`|F(0)| >= zeta(5)/32` no.  Cada una se comprueba por separado y con su propio argumento.

Pero las dos semillas son dos PUNTOS BASE del mismo reticulo afin, de modo que su diferencia es
un vector, y una afirmacion lineal sobre el contenido sobrevive al cambio exactamente cuando su
funcional se anula sobre ese vector.  Eso convierte cuatro comprobaciones separadas en una sola
contraccion, y dice ademas QUE OTRAS afirmaciones sobrevivirian sin tener que probarlas.

QUE SE CALCULA, Y DE DONDE SALE.  Nada se copia de la tabla del articulo: las cinco coordenadas
del sector gauge se derivan de los TRES numeros del corchete de la ec. (68) ---(2,4,7) en la
semilla publicada, (3/2,3,6) en la candidata--- por las mismas definiciones de la S2, y el
resultado se contrasta DESPUES contra los valores impresos.  Si el guion se limitara a repetir
la tabla no probaria nada.  [[circular-artifact-measurement-returns-definition]]

Uso:  python seed_shift_character.py
"""
from fractions import Fraction as F
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

# El corchete gauge de la ec. (68): tres numeros, una fila por (carga, paridad).
#   g1 -> (c=2, s=+1)     g2 -> (c=1, s=+1)     g3 -> (c=1, s=-1)
# y la multiplicidad de la tabla (m,s,c) de la S2 es m = -g/2 en cada fila.
SEMILLAS = {
    "publicada  (a,b)=(2,1/2)":   (F(2), F(4), F(7)),
    "candidata  (a,b)=(3/2,1/2)": (F(3, 2), F(3), F(6)),
}

# Los valores IMPRESOS en la tabla de dos ramas de la S13, para contrastar al final.
IMPRESOS = {
    "publicada  (a,b)=(2,1/2)":   dict(A4=F(-18), D8=F(-27), U2=F(-39), V=F(0), W2=F(-3),
                                       F0=F(9)),
    "candidata  (a,b)=(3/2,1/2)": dict(A4=F(-27, 2), D8=F(-18), U2=F(-30), V=F(0), W2=F(-3),
                                       F0=F(18)),
}


def coordenadas(g1, g2, g3):
    """Las cinco coordenadas del sector gauge, desde el corchete y nada mas.

    Definiciones de la S2:
      A_2 = sum_{s=+1} m c^2      B_2 = sum_{s=-1} m c^2      D = A_2 - (3/4) B_2
      A_4 = sum_{s=+1} m c^4      B_4 = sum_{s=-1} m c^4      A_4^L = sum_{s=+1} m c^4 ln c
      G   = (25/12) A_4 - U ln2 - V ln3,  con  U ln2 + V ln3 = A_4^L + ln2 B_4
      W   = sum_{c impar} m (-s)
      32 F(0)/zeta(5) = S_0 + 31 Delta_0,  con S_0 = sum m, Delta_0 = sum m s
    """
    m1, m2, m3 = -g1 / 2, -g2 / 2, -g3 / 2          # (c=2,+), (c=1,+), (c=1,-)

    A2 = m1 * 4 + m2 * 1
    B2 = m3 * 1
    D = A2 - F(3, 4) * B2

    A4 = m1 * 16 + m2 * 1
    B4 = m3 * 1
    # A_4^L = m1*16*ln2 + m2*1*ln1 = -8 g1 ln2 ; luego U ln2 + V ln3 = (-8 g1 + B4) ln2
    # de modo que U = -8 g1 + B4 y V = 0: el sector gauge no llega a carga tres.
    U = m1 * 16 + B4
    V = F(0)

    # W: solo cargas impares, es decir las dos filas de c=1.  El signo es -s.
    W = m2 * (-1) + m3 * (+1)

    S0 = m1 + m2 + m3
    Delta0 = m1 * (+1) + m2 * (+1) + m3 * (-1)
    F0 = S0 + 31 * Delta0                            # = 32 F(0)/zeta(5)

    return dict(A4=A4, D8=8 * D, U2=2 * U, V=V, W2=2 * W, F0=F0)


P("=" * 88)
P("LAS CINCO COORDENADAS DEL SECTOR GAUGE, DERIVADAS DEL CORCHETE")
P("=" * 88)
calc = {}
malo = 0
for nombre, (g1, g2, g3) in SEMILLAS.items():
    c = coordenadas(g1, g2, g3)
    calc[nombre] = c
    esp = IMPRESOS[nombre]
    ok = all(c[k] == esp[k] for k in esp)
    malo += 0 if ok else 1
    P("  %s   corchete (%s, %s, %s)" % (nombre, g1, g2, g3))
    for k in ("A4", "D8", "U2", "V", "W2", "F0"):
        marca = "ok" if c[k] == esp[k] else "*** DIFIERE del impreso %s ***" % esp[k]
        P("      %-4s = %-8s   %s" % (k, c[k], marca))

if malo:
    P("\nel derivado no reproduce la tabla impresa: no hay nada que decir.")
    raise SystemExit(1)

# --------------------------------------------------------------------------- el vector
pub = calc["publicada  (a,b)=(2,1/2)"]
can = calc["candidata  (a,b)=(3/2,1/2)"]
COORD = ["A4", "D8", "U2", "V", "W2"]
# el reticulo entero es (2A_4, 8D, 2U, V, 2W): A_4 se dobla porque en la candidata es semientero
delta = {"2A4": 2 * (can["A4"] - pub["A4"])}
for k in ("D8", "U2", "V", "W2"):
    delta[k] = can[k] - pub[k]

P("")
P("=" * 88)
P("EL DESPLAZAMIENTO ENTRE LAS DOS SEMILLAS, EN LAS COORDENADAS ENTERAS")
P("=" * 88)
P("  Delta = (2A_4, 8D, 2U, V, 2W) = (%s, %s, %s, %s, %s)"
  % (delta["2A4"], delta["D8"], delta["U2"], delta["V"], delta["W2"]))
P("")
P("  Un enunciado lineal ell sobre el contenido sobrevive al cambio de semilla exactamente")
P("  cuando ell(Delta) se anula en el modulo en que el enunciado esta escrito.")

# --------------------------------------------------------------------------- las filas
# (nombre, funcional sobre (2A_4, 8D, 2U, V, 2W), modulo, sobrevive segun el articulo)
FILAS = [
    ("8D impar",                       (0, 1, 0, 0, 0), 2, False),
    ("A_4 entero  (2A_4 par)",         (1, 0, 0, 0, 0), 2, False),
    ("2U impar",                       (0, 0, 1, 0, 0), 2, False),
    ("2W impar",                       (0, 0, 0, 0, 1), 2, True),
    ("32F(0)/zeta(5) impar",           None,            2, False),   # aparte, ver abajo
    ("ley mod 6:  8D - 2A_4",          (-1, 1, 0, 0, 0), 6, True),
    ("congruencia 8D - 2U (mod 8)",    (0, 1, -1, 0, 0), 8, True),
]

P("")
P("=" * 88)
P("CADA FILA DE LA TABLA DE DOS RAMAS, COMO UNA CONTRACCION CONTRA Delta")
P("=" * 88)
P("  %-30s %-14s %-10s %s" % ("enunciado", "ell(Delta)", "mod", "veredicto"))
P("  " + "-" * 76)
vec = [delta["2A4"], delta["D8"], delta["U2"], delta["V"], delta["W2"]]
fallos = 0
for nombre, ell, mod, sobrevive_art in FILAS:
    if ell is None:
        val = can["F0"] - pub["F0"]
    else:
        val = sum(F(a) * b for a, b in zip(ell, vec))
    invariante = (val % mod == 0)
    acuerdo = (invariante == sobrevive_art)
    fallos += 0 if acuerdo else 1
    P("  %-30s %-14s %-10s %s   %s"
      % (nombre, val, mod,
         "invariante" if invariante else "lo detecta",
         "de acuerdo con el articulo" if acuerdo else "*** EN DESACUERDO ***"))

P("")
if fallos:
    P("%d fila(s) en desacuerdo con lo que dice el articulo." % fallos)
    raise SystemExit(1)

P("=" * 88)
P("Las %d filas salen de una sola contraccion, y ninguna hace falta probarla aparte." % len(FILAS))
P("La ley mod 6 sobrevive por la razon mas fuerte de todas: ell(Delta) es CERO, no cero mod 6,")
P("porque el desplazamiento mueve 2A_4 y 8D EXACTAMENTE LO MISMO, nueve.")
P("Y 2W sobrevive porque el desplazamiento no lo toca en absoluto: lee la diferencia de los dos")
P("pesos de carga uno, y el candidato los mueve a los dos igual.")
P("=" * 88)

# --------------------------------------------------------------------------- control
# Un control que PUEDE fallar: si el desplazamiento fuera otro ---por ejemplo si el candidato
# moviera tambien el peso antiperiodico--- 2W dejaria de ser invariante y el guion tiene que
# decirlo.  Sin esto, el resultado de arriba no distingue una estructura de una casualidad.
P("")
P("CONTROL, y puede fallar: una semilla ficticia que mueva TAMBIEN el peso antiperiodico")
falso = coordenadas(F(3, 2), F(3), F(5))          # g3: 6 -> 5
d2W = 2 * (falso["W2"] / 2 - pub["W2"] / 2) if False else (falso["W2"] - pub["W2"])
P("  corchete (3/2, 3, 5):  Delta(2W) = %s  ->  %s"
  % (d2W, "2W SIGUE invariante -- el control no sirve" if d2W % 2 == 0
     else "2W deja de ser invariante, como debe ser"))
if d2W % 2 == 0:
    raise SystemExit("el control no puede fallar: no distingue nada")
