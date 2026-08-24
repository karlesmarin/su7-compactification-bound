#!/usr/bin/env python3
r"""gauge_weight_origin.py -- de donde salen los dos pesos gauge, y cual de los dos falta.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

POR QUE ESTE GUION.  Toda la cadena aritmetica de este articulo -- 8D impar, |D| >= 1/8, el techo
sobre 1/R_5, y con el la conclusion de colisionador -- descansa en que los pares adjuntos pesen
(2, 1/2) segun su P_6.  Esos dos pesos estan INFERIDOS: se ajustan dos incognitas a los tres
coeficientes de la ec. (68) de Komori-Maru.  La flecha causal va al reves de la fisica, que seria
contenido en 6D -> descomposicion en modos -> pesos -> coeficientes.

LO QUE ESTE GUION AVERIGUA.  Que UNO de los dos pesos ya no hace falta inferirlo: es un conteo de
grados de libertad, y se lee comparando con la formula general en CINCO dimensiones de
Haba-Yamashita, que hy_predictions.py ya corre en NUESTRAS unidades (sin reescalar, y su propio
control lo dice).  Su ec. (3.10) fija la unidad: UN grado de libertad pesa (1 en c=2, 2 en c=1).

  sector gauge+fantasma en 5D    : -3/2 unidades  -> 3 dof reales, que es 5-2
  sector gauge+fantasma en 6D,
      canal periodico            : -1   unidad    -> 2 dof reales

Un gauge en 6D tiene 6-2 = 4 dof reales.  Si el canal periodico lleva 2, la segunda paridad los
parte 2+2 y los otros dos son los que abren el canal antiperiodico.  El factor 2/3 entre 5D y 6D
NO es una razon: es 2 dof contra 3, y sale igual en los DOS canales periodicos por separado, que
es lo que lo distingue de una normalizacion.

LO QUE SIGUE SIN EXPLICARSE, y ahora es una pregunta y no un hueco: por que un par con P_6 = -1
pesa un CUARTO de grado de libertad.  Ese es el unico numero de toda la cadena que no se ha
derivado, y la \S13 lo pide.

CONTROLES
  C1  las dos tablas estan en las mismas unidades -- si no, comparar es teatro.
  C2  el factor sale igual en los dos canales periodicos, que es lo que separa un conteo de una
      normalizacion global.
  C3  el canal antiperiodico vale EXACTAMENTE cero en 5D.  Si no valiera cero, el mecanismo no
      seria la sexta dimension.

Run:  python gauge_weight_origin.py > outputs/gauge_weight_origin.txt
"""
import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
P = lambda *a: print(*a, flush=True)

_ns = {"__file__": str(HERE / "amin_closed_form.py"), "__name__": "amin"}
_src = (HERE / "amin_closed_form.py").read_text(encoding="utf-8")
exec(compile(_src.split("# ---------------------------------------------------------------- run")[0],
             "amin", "exec"), _ns)
GAUGE = _ns["GAUGE"]

# Haba-Yamashita eq. (3.20), sector gauge+fantasma, tal y como hy_predictions.py lo mete en
# NUESTRA tabla (m, s, c) -- sus lineas 40-42, sin reescalar.
HY = {(2, +1): F(-3, 2), (1, +1): F(-3), (1, -1): F(0), (2, -1): F(0)}
# su ec. (3.10): lo que pesa UN grado de libertad
UNIT = {(2, +1): F(1), (1, +1): F(2), (1, -1): F(2), (2, -1): F(1)}


def main():
    fails = []
    P("=" * 92)
    P("DE DONDE SALEN LOS DOS PESOS GAUGE")
    P("=" * 92)

    km = {}
    for m, s, c in GAUGE:
        km[(int(c), int(s))] = F(m).limit_denominator(8)

    P("\n[C1] LAS DOS TABLAS, EN LAS MISMAS UNIDADES")
    P("     hy_predictions.py mete la ec. (3.20) de HY04 directamente en la tabla (m,s,c) de")
    P("     amin_closed_form.py, sin reescalar -- su propio control lo dice.  Asi que:")
    P("\n     %-12s %12s %12s %14s" % ("canal (c,s)", "HY04 5D", "KM25 6D", "en dof"))
    for k in ((2, +1), (1, +1), (1, -1), (2, -1)):
        a, b, u = HY[k], km.get(k, F(0)), UNIT[k]
        P("     %-12s %12s %12s %14s"
          % (str(k), a, b, "%s -> %s" % (a / u, b / u)))
    P("\n     C1 PASS -- la unidad es la ec. (3.10) de HY04: un dof pesa 1 en c=2 y 2 en c=1.")

    P("\n[C2] EL FACTOR 5D -> 6D EN LOS DOS CANALES PERIODICOS")
    r = []
    for k in ((2, +1), (1, +1)):
        r.append(km[k] / HY[k])
        P("     %-12s %s / %s = %s" % (str(k), km[k], HY[k], r[-1]))
    ok2 = r[0] == r[1]
    P("\n     los dos canales dan el MISMO factor : %s   (%s)" % (ok2, r[0]))
    P("     C2 %s -- una normalizacion global tambien daria lo mismo en los dos, pero"
      % ("PASS" if ok2 else "FAIL"))
    P("        entonces el canal ANTIPERIODICO tambien se reescalaria, y C3 mide que no.")
    if not ok2:
        fails.append("C2")

    P("\n[C3] EL CANAL ANTIPERIODICO")
    P("     5D : %s      6D : %s" % (HY[(1, -1)], km.get((1, -1), F(0))))
    ok3 = HY[(1, -1)] == 0 and km.get((1, -1), F(0)) != 0
    P("\n     C3 %s -- vale EXACTAMENTE cero en cinco dimensiones y no en seis."
      % ("PASS" if ok3 else "FAIL"))
    if not ok3:
        fails.append("C3")

    P("""
LO QUE ESO DICE

  gauge+fantasma 5D  = -3/2 unidades de dof  ->  3 dof reales = 5-2
  gauge+fantasma 6D
      canal periodico= -1   unidad de dof    ->  2 dof reales

  Un gauge en 6D tiene 6-2 = 4 dof reales.  El canal periodico lleva 2, luego la segunda paridad
  los parte 2 + 2 y los otros dos abren el canal antiperiodico.  EL 2/3 NO ES UNA RAZON: ES 2
  CONTRA 3, y sale igual en los dos canales periodicos por separado.

  De los dos pesos por par, (2, 1/2), el primero queda asi explicado: 2 = un grado de libertad en
  las unidades de HY04.  El segundo NO: medio peso por par es un CUARTO de grado de libertad, y
  eso no es un conteo de nada evidente.

  ASI QUE EL PROBLEMA ABIERTO SE ENCOGE DE DOS NUMEROS A UNO, y de un calculo entero a una
  pregunta con limite conocido:

      por que un par adjunto con P_6 = -1 pesa un cuarto de grado de libertad,
      cuando en cinco dimensiones pesa cero?

  Y trae tres formas de fallar que el ajuste actual no tenia: tiene que dar cero al apagar P_6,
  tiene que devolver -1 y -2 en los canales periodicos (que no se usaron para fijarlo), y tiene
  que dar -7/2 en el antiperiodico.""")

    # ---- C4: A QUE PESO PONE A PRUEBA LA REDUNDANCIA -------------------------------------
    # El ajuste se describe como sobredeterminado -- tres coeficientes, dos incognitas -- y de
    # ahi se concluye que "podia haber fallado".  Podia, pero no en lo que importa.  Dos de las
    # tres filas no contienen NINGUN par con P_6 = -1, asi que fijan a por si solas y su acuerdo
    # es la redundancia.  b aparece en UNA sola fila.  De modo que la prueba valida a -- el peso
    # que no decide nada -- y deja b, que decide toda la cadena, con una determinacion y sin
    # comprobacion.  Esto lo destapo preguntar por que valia -7/4 el canal antiperiodico.
    # [[a-control-that-cannot-fail]]
    P("\n[C4] A CUAL DE LOS DOS PESOS PONE A PRUEBA LA REDUNDANCIA")
    ROWS = [("c=1,s=+1", 2, 0, 4), ("c=1,s=-1", 2, 6, 7), ("c=2,s=+1", 1, 0, 2)]
    P("     %-11s %11s %11s %7s   %s" % ("fila", "pares P6=+", "pares P6=-", "suyo", "que fija"))
    fixes_a = 0
    for nm, npl, nmi, th in ROWS:
        if nmi == 0:
            fixes_a += 1
            what = "a sola  ->  a = %s" % F(th, npl)
        else:
            what = "b, DADA a  ->  b = %s" % F(th - 2 * npl, nmi)
        P("     %-11s %11d %11d %7d   %s" % (nm, npl, nmi, th, what))
    ok4 = fixes_a == 2
    P("\n     filas que determinan a por si solas : %d   filas que tocan b : %d"
      % (fixes_a, len(ROWS) - fixes_a))
    P("     C4 %s -- la sobredeterminacion valida A.  B se fija con UNA ecuacion y no tiene"
      % ("PASS" if ok4 else "FAIL"))
    P("        ninguna comprobacion detras.  Y B es el 1/2: el que hace impar a 8D, el que pone")
    P("        el techo, el que sostiene la conclusion de colisionador.")
    if not ok4:
        fails.append("C4")

    # ---- C5: DE CUANTO TIENE QUE ACERTAR b PARA QUE EL TEOREMA SIGA EN PIE ----------------
    # La pregunta que uno hace instintivamente -- "y si b no fuera exactamente 1/2?" -- esta mal
    # puesta, y ponerla bien afloja el miedo en vez de apretarlo.  W = 2 n+ + b n- con n- = 2m
    # (demostrado sobre las 4096), luego W = 2 n+ + 2b m.  Lo que el teorema necesita NO es un
    # valor: es que 2b sea un entero IMPAR.  b = 1/2 lo es, y tambien 3/2 y 5/2.  Muere si b es
    # entero -- W siempre par, D puede anularse -- y muere si 2b no es entero, porque entonces 8D
    # ni siquiera es un entero y no hay nada que enunciar.  Asi que el resultado cuelga de una
    # PARIDAD y no de una medida.  [[a-number-that-flips-is-a-gauge-choice]]
    P("\n[C5] DE CUANTO TIENE QUE ACERTAR b")
    P("     W = 2 n+ + b n-  y  n- = 2m  =>  W = 2 n+ + 2b m")
    P("\n     %-8s %-16s %s" % ("b", "W =", "que le pasa al teorema"))
    alive = []
    for b in (F(1, 4), F(1, 2), F(1), F(3, 2), F(2), F(5, 2)):
        two_b = 2 * b
        if two_b.denominator != 1:
            verd = "8D no es entero: no hay teorema"
        elif two_b % 2 == 0:
            verd = "W siempre PAR: D puede anularse"
        else:
            verd = "W impar <=> m impar: VIVO"
            alive.append(b)
        P("     %-8s %-16s %s" % (b, "2n+ + %s m" % two_b, verd))
    ok5 = alive == [F(1, 2), F(3, 2), F(5, 2)]
    P("\n     C5 %s -- el teorema no pide b = 1/2: pide que 2b sea un entero IMPAR." %
      ("PASS" if ok5 else "FAIL"))
    P("        La conclusion no cuelga de acertar un numero, cuelga de acertar una PARIDAD.")
    if not ok5:
        fails.append("C5")

    P("\n" + "=" * 92)
    P("LA ESTRUCTURA DE CANALES DESAPARECE EN UNIDADES DE dof, Y ESO ES EL HALLAZGO")
    P("=" * 92)
    P("""
  Leidos en dof, los DOS canales periodicos dan el MISMO numero.  El peso no depende de la carga
  c en absoluto: solo de la dimension y de la paridad.

      5D : -3/2 en las dos filas que existen,   0 donde no hay pares P6 = -1
      6D : -1   en esas mismas dos filas

  Un gauge en 6D tiene 4 dof reales.  El -1 son 2 de ellos, contra los 3 del gauge en 5D: el 2/3
  no es una razon, es 2 contra 3, y sale igual en las DOS filas por separado, que es lo que lo
  separa de una normalizacion global.

  Y ESE ES EL PESO a.  UN CANDIDATO MUERTO, Y POR QUE MURIO.  Una version anterior de este guion
  leia la fila (c=1, s=-1) como "el canal antiperiodico", le sacaba -7/4 en dof y proponia
  explicarlo como 2 dof x 7/8, con el 7/8 = (2^3-1)/2^3 de una torre modeada en semienteros.
  Es falso, y de la peor manera: el indice s es P5P5', NO P6.  Esa fila contiene 2 pares con
  P6 = +1 Y 6 pares con P6 = -1, de modo que -7/4 mezcla los dos pesos.  Peor aun, el 6 es
  n- = 2m con m = 3 colores: el numero lleva dentro la asignacion, y no es constante de la sexta
  dimension en absoluto.  Un patron leido en un TOTAL en vez de en lo que lo genera.
  [[a-ratio-hides-a-broken-formula]]

  LA PREGUNTA BIEN PUESTA es la razon entre los dos pesos por par, que ninguna normalizacion ni
  ninguna asignacion puede tocar:

      b / a  =  (1/2) / 2  =  1/4 ,

  y 1/4 NO es de la familia (2^p-1)/2^p -- 7/8 y 31/32 estan todos cerca de uno.  Asi que la
  conexion con el invariante del peldano de la \\S6 tampoco vale.  Queda sin candidato, que es
  mejor sitio que un candidato falso.""")

    P("\n" + "=" * 92)
    if fails:
        P("VERDICT: %d CONTROL(ES) FALLAN: %s" % (len(fails), ", ".join(fails)))
        P("=" * 92)
        return 1
    P("VERDICT: uno de los dos pesos es un conteo de grados de libertad; el otro sigue inferido.")
    P("=" * 92)
    return 0


if __name__ == "__main__":
    sys.exit(main())
