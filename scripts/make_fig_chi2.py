#!/usr/bin/env python3
"""make_fig_chi2.py -- el perfil de chi^2 del recast, con sus DOS referencias dibujadas.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Un arbitro pidio exactamente esto, y en una caja: dibujar chi^2(M) - chi^2_SM junto a
chi^2(M) - chi^2_min, para las dos escalas por separado y para el perfilado.  La razon es que el
limite se leia con un umbral de 3.84 medido desde el Modelo Estandar, y eso solo es el punto del
95 % de Wilks si el Modelo Estandar ES el mejor ajuste.  No lo es.

La figura ensena las tres cosas que una tabla no ensena de un vistazo:

  * el minimo esta en un M_KK FINITO, 16 TeV, y no en el infinito.  La preferencia es de
    Delta chi^2 = 4.12 sobre no tener torre, o 2.0 sigma nominales.
  * las dos elecciones de escala de CMS dan la misma forma, de manera que ese minimo no es un
    artefacto de interpolar entre columnas.  Era la objecion, y se contesta mirando.
  * el panel derecho es la misma curva sobre la seleccion M_jj > 3.6 TeV, que quita el bin del
    exceso: el minimo se aplana y la preferencia cae a 0.29.  ESO es lo que dice de donde venia.

  Y las dos lecturas del limite quedan dibujadas donde se cruzan, en vez de citadas.

Run:  python make_fig_chi2.py
"""
import json
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"

ES = "--es" in sys.argv


def T(en, es):
    return es if ES else en


# la paleta de la serie, leida de make_figures_vii.py para no inventar una segunda
HDR = "#1B6F8C"
RED = "#B5530F"
GREEN = "#008C36"
AMBER = "#B07814"
SURF = "#FBFAF7"
STROKE = "#2C2A26"

CEIL_MH = 9.09          # el techo a la masa medida del Higgs
CEIL_WIN = 10.03        # el techo sobre la ventana de KM


def load():
    d = json.loads((OUT / "scan_mkk.json").read_text(encoding="utf-8"))
    for k in ("nlo_all", "nlo_ci"):
        if k not in d:
            raise SystemExit("falta la clave '%s' en scan_mkk.json; corre scan_mkk.py" % k)
    c = json.loads((OUT / "ceiling_ilp.json").read_text(encoding="utf-8"))
    return d, c


def comb_teeth(ci, blk):
    """(1/R_5, Delta chi^2) en cada diente del peine aritmetico.

    El peine son los techos POR PELDANO de ceiling_ilp.py: el mayor 1/R_5 que el programa entero
    consiente en cada 8D.  No es una banda ni una estimacion -- es el conjunto de valores que un
    contenido del bulk puede alcanzar, y esta cuantizado.  Cruzarlo con la curva experimental
    contesta la pregunta propia de este articulo, que no es "cual es la cota inferior" sino
    "donde caen los modelos permitidos sobre la verosimilitud".  [[chain-the-new-brick-to-the-old-one]]

    SOBRE EL CHI^2 DEBIL, que es el de 77 puntos y el que produce el limite citado.  La primera
    version leyo el peine sobre `chi2_abs`, la curva de 70 puntos -- la version FUERTE, la que el
    articulo se niega a usar para el limite por ser la generosa -- y saco 13.4 donde la honesta da
    12.0.  El defecto no era grande pero iba entero en la direccion que nos convenia, que es la
    unica direccion en la que un error propio no se nota.  [[my-errors-all-favour-the-hypothesis]]"""
    if "scan_full" not in blk:
        raise SystemExit("scan_mkk.json no trae 'scan_full': vuelve a correr scan_mkk.py.  Leer "
                         "el peine sobre la curva de 70 puntos daria un numero mas favorable "
                         "que el que el articulo puede defender.")
    sf = np.array(blk["scan_full"])
    out = []
    for r in ci["per_k"]:
        t = r["invR"] / 1000.0
        if t < sf[:, 0].min() or t > sf[:, 0].max():
            continue
        out.append((r["k8D"], r["A4"], t, float(np.interp(t, sf[:, 0], sf[:, 1]))))
    return out


def panel(ax, blk, title, show_ceiling):
    M = np.array([r[0] for r in blk["chi2_abs"]])
    c = np.array([r[1] for r in blk["chi2_abs"]])
    ps = np.array(blk["chi2_per_scale"])
    sm = blk["chi2_SM"]
    cmin = float(c.min())
    Mmin = float(M[int(np.argmin(c))])
    sm0, sm1 = blk["chi2_SM_per_scale"]

    # las dos referencias, que es el asunto entero
    ax.plot(M, c - sm, color=HDR, lw=2.4, zorder=6,
            label=T(r"profiled, referenced to the SM", r"perfilado, referido al SM"))
    ax.plot(M, c - cmin, color=RED, lw=2.4, ls="--", zorder=6,
            label=T(r"profiled, referenced to the minimum",
                    r"perfilado, referido al mínimo"))
    # y las dos escalas de CMS por separado, para que se vea que la forma no es del interpolante
    ax.plot(ps[:, 0], ps[:, 1] - sm0, color=AMBER, lw=1.0, alpha=0.85, zorder=4,
            label=T(r"$\mu=m_{jj}$ alone", r"sólo $\mu=m_{jj}$"))
    ax.plot(ps[:, 0], ps[:, 2] - sm1, color=GREEN, lw=1.0, alpha=0.85, zorder=4,
            label=T(r"$\mu=\langle p_T\rangle$ alone", r"sólo $\mu=\langle p_T\rangle$"))

    # las etiquetas van dentro del recuadro.  Estaban en M.max()*0.995, que es 39.8 TeV mientras
    # el eje llega a 30: se dibujaban flotando fuera del panel.  Un rotulo colocado en
    # coordenadas de datos tiene que respetar el limite del eje, no el del array.
    # [[a-figure-is-an-unaudited-document]]
    for thr, nm, st in ((3.84, r"$3.84$", "-"), (2.71, r"$2.71$", ":")):
        ax.axhline(thr, color=STROKE, lw=0.8, ls=st, alpha=0.55, zorder=2)
        ax.text(0.985, thr, nm, transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=7.4, color=STROKE, alpha=0.8)
    ax.axhline(0.0, color=STROKE, lw=0.7, alpha=0.35, zorder=2)

    # el minimo finito, marcado
    ax.plot([Mmin], [cmin - sm], marker="o", ms=6, color=RED, zorder=8,
            markeredgecolor="white", markeredgewidth=1.1)
    # ARRIBA y no abajo: colgada a -34 puntos caia sobre las marcas del eje x y sobre su rotulo,
    # que dejaba de leerse.  La curva a la derecha del minimo es plana, asi que arriba hay sitio.
    ax.annotate(T("best fit  %.0f TeV,  $\\Delta\\chi^2=%.2f$ vs no tower"
                  % (Mmin, sm - cmin),
                  "mejor ajuste  %.0f TeV,  $\\Delta\\chi^2=%.2f$ frente a sin torre"
                  % (Mmin, sm - cmin)),
                xy=(Mmin, cmin - sm), xytext=(Mmin + 1.2, 16.5), textcoords="data",
                ha="left", va="center", fontsize=7.6, color=RED, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.8, alpha=0.6,
                                shrinkA=2, shrinkB=4,
                                connectionstyle="angle,angleA=0,angleB=90,rad=3"),
                bbox=dict(fc=SURF, ec=RED, lw=0.5, alpha=0.94, pad=2.6,
                          boxstyle="round,pad=0.32"))

    if show_ceiling:
        # la banda es la region que la ARITMETICA excluye, de 0 al techo de la ventana; la linea
        # es el techo a la masa medida del Higgs, que es el numero que un lector debe tomar.
        ax.axvspan(0, CEIL_WIN, color=HDR, alpha=0.09, zorder=1)
        ax.axvline(CEIL_MH, color=HDR, lw=1.2, alpha=0.75, zorder=3)
        ax.text(CEIL_MH - 0.45, 38.0, T("ceiling, measured $m_h$", "techo, $m_h$ medido"),
                rotation=90, va="top", ha="right", fontsize=7.0, color=HDR, alpha=0.95)
        ax.text(CEIL_WIN + 0.4, 38.0,
                T("no content reaches here", "ningún contenido llega aquí"),
                rotation=90, va="top", ha="left", fontsize=6.8, color=HDR, alpha=0.65)

    ax.set_title(title, fontsize=9.6, color=STROKE, fontweight="bold", pad=7)
    ax.set_xlabel(T(r"$1/R_5$  [TeV]", r"$1/R_5$  [TeV]"), fontsize=8.6)
    ax.set_xlim(3.0, 30.0)
    ax.set_ylim(-6.0, 40.0)
    ax.grid(alpha=0.15, lw=0.6)
    for s in ax.spines.values():
        s.set_color(STROKE)
        s.set_alpha(0.35)
    ax.tick_params(labelsize=7.8, colors=STROKE)
    return Mmin, sm - cmin


def main():
    d, ci = load()
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.6), facecolor=SURF)
    for ax in axes:
        ax.set_facecolor(SURF)

    a = panel(axes[0], d["nlo_all"],
              T("All seven mass bins", "Los siete bins de masa"), True)
    b = panel(axes[1], d["nlo_ci"],
              T(r"CMS's own range, $M_{jj}>3.6$ TeV",
                r"El rango propio de CMS, $M_{jj}>3.6$ TeV"), True)
    axes[1].get_legend()  # el panel central no repite la leyenda del izquierdo

    # ---- EL PEINE, en su propio panel y en escala logaritmica ------------------------------
    # Un primer intento lo dibujo ENCIMA del panel izquierdo, cuyo eje llega a 40.  Los dientes
    # van de 13 a 3200, asi que todos menos uno se salian por arriba y sus varillas quedaban
    # como rayas verticales sin sentido.  Tres ordenes de magnitud piden un eje logaritmico y un
    # panel propio.  [[a-figure-is-an-unaudited-document]]
    teeth = comb_teeth(ci, d["nlo_all"])
    ax = axes[2]
    xs = [t[2] for t in teeth]
    ys = [t[3] for t in teeth]
    ax.plot(xs, ys, color=AMBER, lw=1.0, alpha=0.5, zorder=4)
    ax.plot(xs, ys, ls="none", marker="D", ms=6.0, color=AMBER, zorder=9,
            markeredgecolor=STROKE, markeredgewidth=0.6)
    for k8, _A4, t, dc in teeth[:5]:
        ax.annotate("$8D=%d$" % k8, xy=(t, dc), xytext=(0, 9),
                    textcoords="offset points", ha="center", fontsize=6.6, color=STROKE)
    ax.axhline(3.84, color=STROKE, lw=1.0, alpha=0.7, zorder=3)
    ax.text(0.985, 3.84 * 1.13, T(r"$\Delta\chi^2=3.84$", r"$\Delta\chi^2=3.84$"),
            transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=7.2, color=STROKE)
    ax.axhspan(ax.get_ylim()[0] if False else 0.5, 3.84, color=GREEN, alpha=0.10, zorder=1)
    top = max(teeth, key=lambda r: r[2])
    # "cualquier contenido" es una afirmacion sobre TODA la clase, y esa generalidad depende
    # de la semilla: con la candidata el peine arranca en 8D = 2 y el diente mas alto es otro.
    # El pie ya lo dice; el texto DIBUJADO tambien tiene que decirlo, porque es donde se mira.
    ax.annotate(T("relaxation ceiling of the least constrained rung, on the seed of [2]\n"
                  "$%.2f$ TeV --- every content on it sits at $\\Delta\\chi^2\\geq%.1f$"
                  % (top[2], top[3]),
                  "techo de relajación del peldaño menos restringido, semilla de [2]\n"
                  "$%.2f$ TeV --- todo contenido suyo cae en $\\Delta\\chi^2\\geq%.1f$"
                  % (top[2], top[3])),
                xy=(top[2], top[3]), xytext=(0.44, 0.10), textcoords="axes fraction",
                ha="center", va="bottom", fontsize=7.2, color=AMBER, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=AMBER, lw=0.8, alpha=0.75,
                                shrinkA=2, shrinkB=5),
                bbox=dict(fc=SURF, ec=AMBER, lw=0.5, alpha=0.94,
                          boxstyle="round,pad=0.32"))
    ax.set_yscale("log")
    ax.set_ylim(1.0, 8000.0)
    ax.set_xlim(2.6, 11.4)
    ax.set_title(T("Where the allowed models actually fall",
                   "Dónde caen de verdad los modelos permitidos"),
                 fontsize=9.6, color=STROKE, fontweight="bold", pad=7)
    ax.set_xlabel(T(r"$1/R_5$  [TeV]", r"$1/R_5$  [TeV]"), fontsize=8.6)
    ax.set_ylabel(T(r"$\Delta\chi^2$ against the SM", r"$\Delta\chi^2$ frente al SM"),
                  fontsize=8.6)
    ax.grid(alpha=0.15, lw=0.6, which="both")
    for s in ax.spines.values():
        s.set_color(STROKE)
        s.set_alpha(0.35)
    ax.tick_params(labelsize=7.8, colors=STROKE)

    print("  el peine contra la verosimilitud:")
    print("    %6s %8s %12s %12s" % ("8D", "A4", "1/R5 [TeV]", "Delta chi2"))
    for k8, A4, t, dc in teeth:
        print("    %6d %8d %12.3f %12.1f" % (k8, A4, t, dc))
    print("    el minimo Delta chi2 sobre TODO el peine : %.1f, en 8D = %d"
          % (min(r[3] for r in teeth), top[0]))

    # ---- DE QUE DEPENDE ESE MINIMO, que es la pregunta que nadie hace ---------------------
    # Un solo diente decide, el de 8D = 1, porque el techo por peldano cae con k y 1 es el
    # peldano mas pequeno POSIBLE -- y lo es porque el Teorema de los octavos impares dice que 8D
    # es un entero impar.  Ese teorema lo carga un solo medio del sector gauge.  Asi que el margen
    # de la conclusion de colisionador no es 12.0 contra 3.84: es la INTEGRALIDAD de 8D, y esto
    # lo mide.  M^2 va como (6mu+A_4)/k, luego a A_4 fijo un cuanto q veces menor sube 1/R_5 por
    # 1/sqrt(q).  Es cota inferior del efecto: con A_4 libre subiria mas.
    # [[the-cost-is-in-the-step-i-dont-control]]
    sf = np.array(d["nlo_all"]["scan_full"])
    print("\n  DE QUE DEPENDE: si el cuanto de 8D no fuera 1")
    print("    %-14s %12s %14s" % ("cuanto de 8D", "1/R5 [TeV]", "Delta chi2"))
    for q, nm in ((1.0, "1 (el real)"), (0.5, "1/2"), (0.25, "1/4")):
        M = top[2] / q ** 0.5
        print("    %-14s %12.3f %14.1f" % (nm, M, float(np.interp(M, sf[:, 0], sf[:, 1]))))
    print("    umbral: 3.84.  Con el cuanto a la mitad la conclusion no se debilita: CAMBIA DE")
    print("    SIGNO, porque 14.19 TeV cae en la region que los datos levemente prefieren.")
    axes[0].set_ylabel(T(r"$\Delta\chi^2$", r"$\Delta\chi^2$"), fontsize=8.6)
    axes[0].legend(fontsize=7.2, loc="upper right", framealpha=0.92, edgecolor="none")

    # "envelope" y no "arithmetic comb": el peine aritmetico de la S9 vive DENTRO de un
    # peldanyo y tiene paso exacto en M^2.  Esto es otra cosa: el techo de cada peldanyo,
    # uno por peldanyo.  Dos objetos distintos no pueden llevar el mismo nombre.
    fig.suptitle(T("The likelihood, and the per-rung ceiling envelope laid on it",
                   "La verosimilitud, y la envolvente de techos por peldaño encima"),
                 fontsize=11.8, color=STROKE, fontweight="bold", y=0.985)
    fig.text(0.5, 0.005,
             T("Left: the data prefer a finite tower over none by $\\Delta\\chi^2=%.2f$. "
               "Middle: drop the $3.0$--$3.6$ TeV bin --- the one carrying CMS's own "
               "$2.0\\sigma$ local excess --- and the preference collapses to $%.2f$, so it is "
               "theirs and not ours.  Right: the question this paper can ask and a continuous "
               "lower limit cannot --- the relaxation ceiling of each reachable rung, "
               "against the same likelihood.  Contents sit at or below theirs, so by "
               "monotonicity these are conservative lower bounds and not attained points."
               % (a[1], b[1]),
               "Izquierda: los datos prefieren una torre finita a ninguna por "
               "$\\Delta\\chi^2=%.2f$. Centro: al quitar el bin de $3.0$--$3.6$ TeV ---el que "
               "lleva el exceso local de $2.0\\sigma$ de CMS--- la preferencia se desploma a "
               "$%.2f$, así que es suya y no nuestra.  Derecha: la pregunta que este artículo "
               "puede hacer y una cota inferior continua no --- el techo de relajación "
               "de cada peldaño alcanzable, contra la misma verosimilitud.  Los contenidos se "
               "sientan igual o por debajo del suyo, de modo que por monotonía son cotas "
               "inferiores conservadoras y no puntos alcanzados." % (a[1], b[1])),
             ha="center", va="bottom", fontsize=7.4, color=STROKE, wrap=True)

    fig.tight_layout(rect=(0, 0.045, 1, 0.955))
    name = "fig_chi2_es.pdf" if ES else "fig_chi2.pdf"
    for d2 in (HERE, HERE / "paper"):
        fig.savefig(d2 / name, facecolor=SURF)
    fig.savefig(HERE / name.replace(".pdf", ".png"), dpi=145, facecolor=SURF)
    print("  %s   minimo %.1f TeV, preferencia %.2f  |  CI: %.1f TeV, %.2f"
          % (name, a[0], a[1], b[0], b[1]))
    return 0


def both():
    """LAS DOS EDICIONES EN UNA INVOCACION.

    La salida archivada se hizo corriendo el guion dos veces, `>` y luego `--es >>`, y
    check_reproduces.py lo corre UNA. El archivo tenia por tanto el doble de lineas que la
    corrida que lo comprueba, y la compuerta lo marco como que difiere -- con razon: un archivo
    que sale de una invocacion distinta de la que se comprueba no comprueba nada.  Es la misma
    leccion que scan_mkk.both() y make_figures_vii.  [[a-patch-on-disk-is-not-a-patch-that-ran]]"""
    global ES
    rc = 0
    for es in (False, True):
        ES = es
        rc |= main()
    return rc


if __name__ == "__main__":
    sys.exit(both())
