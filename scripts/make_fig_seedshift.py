#!/usr/bin/env python3
"""make_fig_seedshift.py -- las dos semillas son dos puntos base, y lo que las separa es un vector.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

QUE DIBUJA, Y POR QUE EN TRES DIMENSIONES.  La S13 dice, en prosa, que la tabla de dos ramas
entera es una sola contraccion contra

    Delta_seed = (2A_4, 8D, 2U, V, 2W) = (9, 9, 9, 0, 0)

y que un enunciado sobrevive al cambio de semilla exactamente cuando su funcional se anula sobre
ese vector.  Eso es geometria, y en prosa cuesta media pagina.

La clave es que TODO GENERADOR DE MATERIA ES PAR en las tres coordenadas que deciden paridades
---(2A_4, 8D, 2W)---, de modo que ningun contenido de bulk puede cambiar ninguna de las tres: la
clase de paridad la fija el PUNTO BASE y nada mas.  Las clases de paridad son entonces los ocho
vertices de un cubo, cada semilla ocupa uno, y el desplazamiento candidato es una arista.

  panel izquierdo, 3D : el cubo de las ocho clases de paridad.  La semilla publicada ocupa el
      vertice (par, IMPAR, IMPAR) y la candidata el (IMPAR, par, IMPAR).  La arista entre las
      dos se mueve en 2A_4 y en 8D y NO en 2W: por eso el teorema del 2W sobrevive y el de los
      octavos impares no.  Los ocho generadores de materia se dibujan como lazos sobre el
      vertice, porque eso es literalmente lo que hacen: no lo mueven.
  panel derecho       : las cinco componentes de Delta y, debajo de cada una, que enunciado del
      articulo vive en ella.

NADA SE TRANSCRIBE.  Los generadores y el punto base salen de `congruences.py`, que ya los
construye desde las tablas de terminos, y las coordenadas de la semilla candidata salen del
corchete gauge por las mismas definiciones, igual que en `seed_shift_character.py`.

Uso:  python make_fig_seedshift.py [--es]
"""
import os
import sys
from fractions import Fraction as F

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ES = "--es" in sys.argv
T = lambda en, es: es if ES else en
SUF = "_es" if ES else ""
P = lambda *a: print(*a, flush=True)

BLUE, AMBER, RED = "#3A86C8", "#E0A030", "#C0392B"
STROKE, INK, MUTED = "#1F4E79", "#1F2933", "#6B7280"
GRID, SURF = "#E5E7EB", "#FCFCFB"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.facecolor": SURF, "figure.facecolor": SURF,
})


# --------------------------------------------------------------------------- los datos
def generadores_y_base():
    """(nombres, generadores, base publicada) en las CINCO coordenadas, desde congruences.py.

    Se importa en vez de retecleare: la tabla de terminos vive alli y ya esta comprobada contra
    las cinco filas publicadas.  Su salida se silencia porque este guion tiene la suya.
    """
    sys.path.insert(0, HERE)
    real, cwd = sys.stdout, os.getcwd()
    os.chdir(HERE)
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    try:
        import congruences as C
    finally:
        sys.stdout.close()
        sys.stdout, _ = real, os.chdir(cwd)
    return list(C.NAMES), [list(m) for m in C.MATTER], list(C.GAUGE5)


def gauge_desde_corchete(g1, g2, g3):
    """Las cinco coordenadas del sector gauge desde los tres numeros de la ec. (68).

    Mismas definiciones que seed_shift_character.py: m = -g/2 por fila, y las filas son
    (c=2, s=+1), (c=1, s=+1), (c=1, s=-1).
    """
    m1, m2, m3 = -F(g1) / 2, -F(g2) / 2, -F(g3) / 2
    A2, B2 = m1 * 4 + m2, m3
    A4, B4 = m1 * 16 + m2, m3
    U = m1 * 16 + B4
    W = m2 * (-1) + m3 * (+1)
    return [A4, 8 * (A2 - F(3, 4) * B2), 2 * U, F(0), 2 * W]


NOMBRES, MATTER, GAUGE_PUB = generadores_y_base()
GAUGE_CAN = gauge_desde_corchete(F(3, 2), F(3), F(6))

# el reticulo entero usa 2A_4, porque en la semilla candidata A_4 es semientero
def a5(v):
    return [2 * F(v[0]), F(v[1]), F(v[2]), F(v[3]), F(v[4])]


PUB, CAN = a5(GAUGE_PUB), a5(GAUGE_CAN)
MAT5 = [a5(m) for m in MATTER]
DELTA = [c - p for c, p in zip(CAN, PUB)]
COORD = ["2A_4", "8D", "2U", "V", "2W"]

P("=" * 88)
P("LAS DOS SEMILLAS, EN LAS COORDENADAS ENTERAS (2A_4, 8D, 2U, V, 2W)")
P("=" * 88)
P("  publicada : %s" % [str(x) for x in PUB])
P("  candidata : %s" % [str(x) for x in CAN])
P("  Delta     : %s" % [str(x) for x in DELTA])

# --------------------------------------------------------------------------- controles
# (1) TODO generador de materia par en las tres coordenadas de paridad. Si uno fuera impar, la
#     clase de paridad dependeria del contenido y no del punto base, y esta figura no diria nada.
TRES = [0, 1, 4]                                   # 2A_4, 8D, 2W
P("")
P("  CONTROL 1: .es par todo generador de materia en (2A_4, 8D, 2W)?")
impares = [(n, [str(m[i]) for i in TRES]) for n, m in zip(NOMBRES, MAT5)
           if any(m[i] % 2 for i in TRES)]
for n, v in impares:
    P("     %-10s %s   IMPAR" % (n, v))
P("     %s" % ("los ocho son pares: la clase de paridad la fija el punto base"
               if not impares else "*** hay generadores impares: la figura no vale ***"))
if impares:
    raise SystemExit(1)

# (2) UN CONTROL QUE PUEDE FALLAR.  Una semilla ficticia que mueva tambien el peso antiperiodico
#     tiene que caer en OTRO vertice del cubo, o el dibujo no distingue nada.
FALSA = a5(gauge_desde_corchete(F(3, 2), F(3), F(5)))
cel = lambda v: tuple(int(v[i]) % 2 for i in TRES)
P("")
P("  CONTROL 2, y puede fallar: una semilla ficticia de corchete (3/2, 3, 5)")
P("     celda de la publicada : %s" % (cel(PUB),))
P("     celda de la candidata : %s" % (cel(CAN),))
P("     celda de la ficticia  : %s   %s"
  % (cel(FALSA),
     "distinta de las dos, como debe ser" if cel(FALSA) not in (cel(PUB), cel(CAN))
     else "*** cae en la misma: el cubo no separa nada ***"))
if cel(FALSA) in (cel(PUB), cel(CAN)):
    raise SystemExit("el control no puede fallar")

CEL_PUB, CEL_CAN = cel(PUB), cel(CAN)
P("")
P("  la arista publicada -> candidata se mueve en: %s"
  % ", ".join(COORD[i] for i in TRES if CEL_PUB[TRES.index(i)] != CEL_CAN[TRES.index(i)]))
P("  y NO se mueve en   : %s"
  % ", ".join(COORD[i] for i in TRES if CEL_PUB[TRES.index(i)] == CEL_CAN[TRES.index(i)]))


# --------------------------------------------------------------------------- la figura
def draw():
    fig = plt.figure(figsize=(11.6, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.10)
    ax = fig.add_subplot(gs[0], projection="3d")
    ax.set_facecolor(SURF)

    # las doce aristas del cubo de clases de paridad
    verts = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    edges = [(u, v) for u in verts for v in verts if sum(abs(x - y) for x, y in zip(u, v)) == 1]
    ax.add_collection3d(Line3DCollection([(u, v) for u, v in edges],
                                         colors="#C9D3DD", linewidths=1.0, zorder=1))
    for v in verts:
        ocupada = v in (CEL_PUB, CEL_CAN)
        if not ocupada:
            ax.scatter(*v, s=26, c="white", edgecolors=MUTED, linewidths=0.8, zorder=2,
                       depthshade=False)

    # el segmento entre las dos semillas.  NO es una arista: cambia DOS coordenadas (2A_4 y 8D)
    # y deja 2W quieta, que es justamente lo que la figura tiene que ensenyar.
    ax.plot(*zip(CEL_PUB, CEL_CAN), color=AMBER, lw=3.6, zorder=5, solid_capstyle="round")

    for cel_, col, lab in ((CEL_PUB, BLUE, T("published seed", "semilla publicada")),
                           (CEL_CAN, RED, T("candidate seed", "semilla candidata"))):
        ax.scatter(*cel_, s=240, c=col, edgecolors="white", linewidths=1.8, zorder=6,
                   depthshade=False, label=lab)

    # LOS ROTULOS VAN EN COORDENADAS DE EJE, no de datos.  Colocados en 3D con ax.text() se
    # proyectan donde la vista los manda: el $\Delta_{\rm seed}$ caia encima del punto azul y
    # el «y no en 2W» encima de la etiqueta del eje.  Un rotulo que se pisa con otro no esta
    # dicho.  [[verify-visuals-headless]]
    ax.text2D(0.325, 0.745, r"$\Delta_{\rm seed}$", transform=ax.transAxes,
              color="#8a6512", fontsize=12.5, fontweight="bold", ha="left", va="center")
    ax.text2D(0.585, 0.90, T("moves in $2A_4$ and $8D$", "se mueve en $2A_4$ y $8D$"),
              transform=ax.transAxes, color="#8a6512", fontsize=8.6, style="italic", ha="left")
    ax.text2D(0.585, 0.845, T("and not in $2W$", "y no en $2W$"),
              transform=ax.transAxes, color=STROKE, fontsize=8.6, style="italic", ha="left")

    ax.set_xticks([0, 1]); ax.set_yticks([0, 1]); ax.set_zticks([0, 1])
    par, imp = T("even", "par"), T("odd", "impar")
    ax.set_xticklabels([par, imp], fontsize=8.2)
    ax.set_yticklabels([par, imp], fontsize=8.2)
    ax.set_zticklabels([par, imp], fontsize=8.2)
    ax.set_xlabel("$2A_4$", labelpad=-2, fontsize=10)
    ax.set_ylabel("$8D$", labelpad=-2, fontsize=10)
    ax.set_zlabel("$2W$", labelpad=-2, fontsize=10)
    ax.set_title(T("The eight parity classes, and the seed picks one",
                   "Las ocho clases de paridad, y la semilla elige una"),
                 fontsize=9.8, color=INK, fontweight="bold", pad=0)
    ax.view_init(elev=18, azim=-58)
    try:
        ax.set_box_aspect((1, 1, 0.92))
    except Exception:                                   # noqa: BLE001
        pass
    ax.legend(loc="upper left", fontsize=8.2, frameon=False, bbox_to_anchor=(-0.02, 0.97))
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor("white")
        a.pane.set_alpha(0.30)
        a._axinfo["grid"].update(color=GRID, linewidth=0.5)

    # la frase que explica por que la celda basta va DEBAJO del panel, no dentro: dentro se
    # cruzaba con las aristas y con el segmento.  Un rotulo que no se lee no esta dicho.
    ax.text2D(0.5, -0.045,
              T("every matter generator is even in all three coordinates, so no bulk content\n"
                "leaves its cell: the parity class is the base point's alone",
                "todo generador de materia es par en las tres coordenadas, así que ningún\n"
                "contenido sale de su celda: la clase de paridad es del punto base y de nadie más"),
              transform=ax.transAxes, ha="center", va="top", fontsize=7.8,
              color=MUTED, style="italic")

    # ---- derecha: las cinco componentes de Delta y que vive en cada una
    ax2 = fig.add_subplot(gs[1])
    vals = [int(x) for x in DELTA]
    etiquetas = [r"$2A_4$", r"$8D$", r"$2U$", r"$V$", r"$2W$"]
    vive = [T("$A_4$ integral", "$A_4$ entero"),
            T("odd eighths, Thm 1", "octavos impares, Teo. 1"),
            T("$2U$ odd", "$2U$ impar"),
            T("(nothing here)", "(nada aquí)"),
            T("$2W$ odd, and $|F(1)-F(0)|$", "$2W$ impar, y $|F(1)-F(0)|$")]
    cols = [RED if v % 2 else BLUE for v in vals]
    y = np.arange(len(vals))[::-1]
    ax2.barh(y, vals, height=0.52, color=cols, edgecolor=SURF, linewidth=1.4, zorder=3)
    for yy, v, w, c in zip(y, vals, vive, cols):
        ax2.text((v if v else 0) + 0.35, yy, "%d   %s" % (v, w), va="center", fontsize=8.4,
                 color=c if v % 2 else STROKE)
    ax2.set_yticks(y); ax2.set_yticklabels(etiquetas, fontsize=10.5)
    ax2.set_xlim(-0.6, 16.0)
    ax2.set_xlabel(T(r"component of $\Delta_{\rm seed}$",
                     r"componente de $\Delta_{\rm seed}$"))
    ax2.set_title(T("Odd component: the claim flips.  Zero: it survives.",
                    "Componente impar: el enunciado voltea.  Cero: sobrevive."),
                  fontsize=9.8, color=INK, fontweight="bold", pad=8)
    ax2.grid(axis="x", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.text(0.5, -0.20,
             T(r"and $8D-2A_4$ picks up $9-9=0$: the mod-six law survives for the strongest"
               "\n" r"reason of all --- not zero modulo six, but zero.",
               r"y $8D-2A_4$ recoge $9-9=0$: la ley módulo seis sobrevive por el motivo más"
               "\n" r"fuerte de todos: no cero módulo seis, sino cero."),
             transform=ax2.transAxes, ha="center", va="top", fontsize=7.9,
             color=MUTED, style="italic")

    stem = "fig_seedshift" + SUF
    for d in (HERE, os.path.join(HERE, "paper")):
        fig.savefig(os.path.join(d, stem + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, stem + ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    P("")
    P("  escrito: %s.pdf (aqui y en paper/) y %s.png" % (stem, stem))


draw()
