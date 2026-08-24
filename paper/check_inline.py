#!/usr/bin/env python3
"""check_inline.py - la matematica EN LINEA dice lo mismo en las dos ediciones.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

POR QUE EXISTE.  `check_formulas.py` empareja las formulas de DISPLAY de las dos ediciones y las
compara caracter a caracter.  Todo lo demas -- las celdas de tabla y las formulas dentro de la
prosa -- queda fuera, y ahi cabe una errata que ninguna otra compuerta ve.

El 24-ago-2026, al leer la edicion castellana de cabo a rabo, aparecio esta, en la tabla de la
escalera de la S4:

    EN   $(S_2+\\Delta_2)/2=\\mathcal{A}_2=A_4$   --- periodic only
    ES   $(S_2+\\Delta_2)/2=A_2$                 --- solo la periodica

es decir, la castellana decia que el CUARTO momento es $A_2$, que es el SEGUNDO.  Se contradecia
con su propia S8, que escribe $\\mathcal{A}_2=A_4$ tres paginas mas abajo.  Llevaba ahi desde que
se escribio la tabla.  [[a-table-row-is-one-claim]]

COMO FUNCIONA.  Extrae toda la matematica `$...$` de las dos ediciones (fuera displays, fuera
`\\texttt`), normaliza el espaciado y la prosa traducida de dentro de `\\text{}`, y pregunta que
formula existe en una edicion y NO EXISTE EN ABSOLUTO en la otra.  No compara cuentas ni orden:
que una formula salga tres veces en un idioma y dos en el otro es redaccion, porque la misma
frase se parte distinto al traducir.  Lo que delata una errata es la AUSENCIA.

LO QUE NO PUEDE HACER.  No sabe traducir.  Si una edicion escribe `$M_{jj}=5$--$8$~TeV` y la otra
`$M_{jj}$ entre $5$ y $8$~TeV`, lo marca, y es un falso positivo.  Por eso los casos comprobados
a mano viven en ACEPTADOS, con su razon escrita: una lista blanca vacia deja la compuerta
inservible, y una lista blanca sin razones deja de ser una lista blanca.  Tampoco mira simbolos
sueltos ---$G$, $m_h$, $-18$---, que son puro ruido de redaccion.

Uso:  python check_inline.py     (desde part_vii/paper/)
"""
import collections
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
EN, ES = "su7_hierarchy.tex", "su7_hierarchy_es.tex"
P = lambda *a: print(*a, flush=True)

# Diferencias comprobadas a mano el 24-ago-2026: la frase se parte distinto al traducir y la
# matematica dice lo mismo.  Cada una con su razon.
ACEPTADOS = {
    r"-27\to-18": "la ES escribe «pasa de $-27$ a $-18$», con los dos numeros sueltos",
    r"-39\to-30": "idem, «$2U$ de $-39$ a $-30$»",
    r"M_{jj}=5": "la ES escribe «$M_{jj}$ entre $5$ y $8$~TeV»",
    r"(TEXT,TEXT)": "«(representation, parity)» va en matematicas en la EN y en prosa en la ES",
}

# La errata real que hizo falta esta compuerta, guardada para que no se pueda volver a colar:
# tiene que salir marcada.  [[a-falsification-suite-beats-a-passing-control]]
FALSIFY = (r"(S_2+\Delta_2)/2=A_2", r"(S_2+\Delta_2)/2=\mathcal{A}_2=A_4")


def inline(path):
    s = io.open(path, encoding="utf-8").read()
    s = re.sub(r"(?<!\\)%.*", "", s)
    s = s.split(r"\begin{thebibliography}")[0]
    s = re.sub(r"\\begin\{equation\}.*?\\end\{equation\}", " ", s, flags=re.S)
    # OJO con el lookbehind: sin el, el `\\[3pt]` de los nodos TikZ casa como apertura de display
    # y el `.*?` se traga en silencio la tabla que viene detras.  Asi fue como la primera version
    # de este barrido dijo «sin diferencias» sobre la celda que si difiere.
    # [[falsify-the-instrument-first]]
    s = re.sub(r"(?<!\\)\\\[.*?(?<!\\)\\\]", " ", s, flags=re.S)
    s = re.sub(r"\\texttt\{(?:[^{}]|\{[^{}]*\})*\}", " ", s)
    out = []
    for m in re.finditer(r"(?<!\\)\$([^$]+)\$", s):
        f = " ".join(m.group(1).split())
        # la prosa de dentro va traducida a proposito: no es matematica
        f = re.sub(r"\\(text|mbox|mathrm)\{(?:[^{}]|\{[^{}]*\})*\}", "TEXT", f)
        for junk in (r"\,", r"\;", r"\ ", r"\!", "~"):
            f = f.replace(junk, "")
        f = re.sub(r"\s+", "", f)          # el espaciado dentro de $...$ no es contenido
        f = f.replace(r"\alpha_{\min}", r"\amin")   # la macro y su expansion son lo mismo
        # Un simbolo suelto -- $G$, $m_h$, $\ln2$, $-18$ -- aparece un numero distinto de veces
        # en cada idioma por pura redaccion, y compararlo solo produce ruido.  Lo que puede
        # esconder una errata es una formula con estructura.
        if len(f) >= 6:
            out.append(f)
    return collections.Counter(out)


def main():
    en = inline(os.path.join(HERE, EN))
    es = inline(os.path.join(HERE, ES))
    P("=" * 96)
    P("LA MATEMATICA EN LINEA, EDICION CONTRA EDICION")
    P("  EN %d formulas en linea (%d distintas)   ES %d (%d distintas)"
      % (sum(en.values()), len(en), sum(es.values()), len(es)))
    P("=" * 96)

    # AUSENTE, no «menos veces».  Que una formula salga tres veces en una edicion y dos en la
    # otra es redaccion: la misma frase se parte distinto.  Lo que delata una errata es que la
    # formula no exista EN ABSOLUTO al otro lado -- asi estaba el `(S_2+\Delta_2)/2=A_2` de la
    # tabla de la escalera, que en ingles no aparece ni una vez.
    bad = []
    for label, mine, other in (("solo en la inglesa", en, es), ("solo en la castellana", es, en)):
        rest = sorted(f for f in mine if f not in other and f not in ACEPTADOS)
        P("\n  %-24s %d formula(s) que no existen al otro lado" % (label, len(rest)))
        for f in rest:
            P("     %2dx  %s" % (mine[f], f[:86]))
            bad.append((label, f))

    P("\n  lista blanca (%d), comprobadas a mano:" % len(ACEPTADOS))
    for f, why in ACEPTADOS.items():
        P("     %-22s %s" % (f[:22], why))

    # falsacion: la errata del 24-ago tiene que seguir siendo visible para esta compuerta
    wrong, right = FALSIFY


    seen_wrong = (wrong in en) or (wrong in es)
    P("\n  falsacion: la errata del cuarto momento")
    P("     la version correcta   %-38s %s" % (right[:38], "presente" if (right in en and right in es) else "AUSENTE en una edicion"))
    P("     la version erronea    %-38s %s" % (wrong[:38], "ausente, bien" if not seen_wrong else "PRESENTE, ha vuelto"))
    if seen_wrong or not (right in en and right in es):
        bad.append(("falsacion", wrong))

    P("\n" + "=" * 96)
    if bad:
        P("MATEMATICA EN LINEA QUE NO CASA: %d" % len(bad))
        P("Una celda de tabla es una afirmacion y nadie la compara. Miralas una a una: o es una")
        P("errata en una de las dos ediciones, o es un reparto distinto de la frase y va a la")
        P("lista blanca CON SU RAZON.")
        P("=" * 96)
        return 1
    P("la matematica en linea de las dos ediciones dice lo mismo.")
    P("=" * 96)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
