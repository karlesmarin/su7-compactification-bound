#!/usr/bin/env python3
"""check_narrative.py -- does the paper agree with itself where it summarises itself?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

check_numbers.py asks whether every printed number has an archived run behind it.  That is a
different question from this one.  A number can be perfectly archived and still be the WRONG one
here: the abstract and the introduction are written first and edited last, so they are where a
superseded value survives after the body has moved on.  Nothing else in the pipeline compares the
summary against the thing summarised.

Three checks, and each of them caught something real the first time it was run:

  N1  every number the abstract prints appears somewhere in the body.  A number that appears in
      the abstract and NOWHERE else is either a typo or a value the body has since changed.
      (Found: the abstract rounding the closed-form median to 0.21 where the body and the archive
      both say 0.209.)
  N2  the same for the introduction, which carries the list of what is new and therefore the
      largest concentration of forward promises in the paper.
  N3  every \\S the abstract or the introduction points at exists.  A promise to a section that
      was renamed is a dangling promise, and \\ref prints a bare ?? that a reader reads as sloppy
      rather than as broken.

WHAT THIS DOES NOT CHECK, and it is the harder half: whether the abstract's PROSE contradicts the
body's.  That found the worst defect of the 23 August pass -- a paragraph reading "the recast is
now done, so the word may be exclusion" six lines above "why that is not yet the word excluded" --
and it found it by being read, not by being run.  This gate is a floor, not a substitute.

Numbers that legitimately live only in the abstract go in ALLOW with a reason.

Run:  python check_narrative.py     (from part_vii/paper/)
"""
import pathlib
import re
import sys

TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]

# keep the decimal point and match on token boundaries, the lesson check_numbers.py carries:
# a substring match passes `1.816` on the strength of `816`.
NUM = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)(?!\w)(?!\.\d)")
REF = re.compile(r"\\ref\{(sec:[a-z]+)\}")
LAB = re.compile(r"\\label\{(sec:[a-z]+)\}")

ALLOW = {
    # TikZ geometry inside the chain figure, which sits inside the introduction: millimetres and
    # per-cent colour mixes, not claims.  ONLY the two that the raw run actually reports go here.
    # The first version of this set also carried "2", "4", "7" and "52" "to be safe", which is a
    # hole and not safety: those pass on their own, and an allowlist wide enough to be comfortable
    # is wide enough to swallow a real claim.  [[a-control-that-cannot-fail]]
    "35", "60",
}


def sections(s):
    """(abstract, introduction, everything else) -- the two summaries against the whole."""
    a0, a1 = s.index(r"\begin{abstract}"), s.index(r"\end{abstract}")
    i0, i1 = s.index(r"\label{sec:intro}"), s.index(r"\label{sec:setting}")
    body = s[:a0] + s[a1:i0] + s[i1:]
    return s[a0:a1], s[i0:i1], body


def report(name, part, body, fails):
    have = set(NUM.findall(body))
    miss = sorted({t for t in NUM.findall(part) if t not in have and t not in ALLOW},
                  key=lambda t: (len(t), t))
    print("  %-14s %d distinct numbers, %d absent from the rest of the paper"
          % (name, len(set(NUM.findall(part))), len(miss)))
    for t in miss:
        j = part.find(t)
        print("     %-9s <-- %s" % (t, " ".join(part[max(0, j - 64):j + 36].split())))
        fails.append("%s:%s" % (name, t))


# --------------------------------------------------------------------------- N6
# NUMEROS ESCRITOS CON LETRA.  N5 compara solo los que llevan punto decimal, y lo dice: los
# enteros sueltos entran en giros de lengua.  Pero un numero escrito CON LETRA Y PEGADO A UN
# SUSTANTIVO TECNICO no es un giro, es una afirmacion, y ninguna compuerta lo miraba.
#
# El 24-ago-2026 la S13 INGLESA decia "$D=0$ is reached by TWO multiplets" donde su propio
# resumen, su propia S7 y las tres apariciones castellanas dicen TRES ---y el testigo de la S7
# los enumera: un adjunto mas dos fundamentales---.  Tres repasos externos y doce compuertas
# pasaron por encima, porque "two" no es un decimal y no lo miraba nadie.
# [[the-conclusion-line-is-a-prediction]] [[no-claim-lives-in-one-place]]
#
# Para no hacer ruido, solo se comparan los pares <numeral, sustantivo tecnico>: "two multiplets"
# contra "tres multipletes".
#
# Y HACE FALTA MAS QUE ESO, porque la primera version devolvia 39 diferencias y ninguna era real:
#   - "one" y "un/una" son ARTICULOS antes que numerales, y arrasaban con el 1;
#   - "the two branches" y "las dos ramas" son "both", no un recuento, y arrasaban con el 2;
#   - "twenty-three multiplets" casa `\b(three)[ -]multiplets` por la frontera del guion, mientras
#     que "veintitres multipletes" va en una palabra y no casa: los numerales COMPUESTOS metian
#     una diferencia sistematica en cada valor.
# Una compuerta que devuelve 39 falsos positivos no la mira nadie, y una lista blanca lo bastante
# ancha para acallarla se traga tambien el fallo real.
#
# LA CURA FUE PEOR UN MOMENTO.  El primer arreglo subio el minimo a 3 para matar el 1 y el 2 --- y
# con eso dejo fuera el "TWO multiplets" que es la razon de existir de esta comprobacion.  Lo
# canto el bloque de falsacion de abajo, no la lectura.  Lo que de verdad quita el ruido no es
# subir el minimo: es comparar CONJUNTOS DE VALORES por parrafo en vez de CUENTAS por documento,
# porque "the two branches" y "las dos ramas" aportan el mismo valor a las dos ediciones y se
# cancelan solos.  Con eso el minimo puede bajar a 2, que es donde vivia el fallo.
# El 1 se queda fuera: "one"/"un"/"una" son articulos antes que numerales.
# [[a-gate-can-forbid-its-own-remedy]] [[a-control-that-cannot-fail]]
MIN_VAL = 2
WORDS_EN = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8,
    "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14,
    "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
    "thirty": 30, "forty": 40, "fifty": 50, "sixty": 60, "seventy": 70, "eighty": 80,
    "ninety": 90, "hundred": 100,
}
WORDS_ES = {
    "un": 1, "uno": 1, "una": 1, "dos": 2, "tres": 3, "cuatro": 4, "cinco": 5, "seis": 6,
    "siete": 7, "ocho": 8, "nueve": 9, "diez": 10, "once": 11, "doce": 12, "trece": 13,
    "catorce": 14, "quince": 15, "dieciséis": 16, "diecisiete": 17, "dieciocho": 18,
    "diecinueve": 19, "veinte": 20, "treinta": 30, "cuarenta": 40, "cincuenta": 50,
    "sesenta": 60, "setenta": 70, "ochenta": 80, "noventa": 90, "cien": 100, "ciento": 100,
}
# concepto -> (sustantivo ingles, sustantivo castellano).  Solo objetos que se CUENTAN y cuyo
# recuento es una afirmacion del articulo.
NOUNS = [
    ("multiplete", r"multiplets?", r"multipletes?"),
    ("peldano", r"rungs?", r"pelda(?:ñ|n)os?"),
    ("semilla", r"seeds?", r"semillas?"),
    ("coordenada", r"coordinates?", r"coordenadas?"),
    ("congruencia", r"congruences?", r"congruencias?"),
    ("dimension", r"dimensions?", r"dimensiones?|dimensi(?:ó|o)n"),
    ("momento", r"moments?", r"momentos?"),
    ("generador", r"generators?", r"generadores?|generador"),
    ("paridad", r"parities|parity", r"paridades?"),
    ("rama", r"branch(?:es)?", r"ramas?"),
    ("techo", r"ceilings?", r"techos?"),
    ("corte", r"cuts?", r"cortes?|corte"),
    ("lazo", r"loops?", r"lazos?"),
    ("color", r"colours?", r"colores?|color"),
    ("vertice", r"vertices|vertex", r"v(?:é|e)rtices?"),
    ("modo", r"modes?", r"modos?"),
    ("torre", r"towers?", r"torres?"),
    ("condicion", r"conditions?", r"condiciones?|condici(?:ó|o)n"),
    ("identidad", r"identities|identity", r"identidades?|identidad"),
    ("teorema", r"theorems?", r"teoremas?"),
    ("prueba de cierre", r"closure tests?", r"pruebas? de cierre"),
]


# Diferencias comprobadas a mano que NO son defectos.  Cada una con su razon escrita: una lista
# blanca sin razones deja de ser una lista blanca.  [[a-control-that-cannot-fail]]
N6_ACEPTADOS = {
    # El ingles escribe "in five dimensions --- where they differ" y el castellano elide el
    # sustantivo, "y en cinco ---donde difieren---", porque la frase anterior acaba de decir
    # "En seis dimensiones".  La elipsis es correcta en castellano y no hay nada que arreglar.
    ("dimension", "en cinco"),
}


def _paras(path):
    """Los parrafos de prosa, en orden, sin flotantes ni comentarios.

    Los comentarios se BORRAN enteros en vez de dejar la linea en blanco: si se deja el hueco,
    un bloque de notas dentro de un parrafo lo parte en dos y el emparejamiento se descuadra a
    partir de ahi.
    """
    s = open(path, encoding="utf-8").read().split(r"\begin{thebibliography}")[0]
    s = "\n".join(ln for ln in s.split("\n") if not ln.lstrip().startswith("%"))
    s = re.sub(r"(?<!\\)%.*", "", s)
    for env in ("figure", "table", "tikzpicture", "longtable", "tabular"):
        s = re.sub(r"\\begin\{" + env + r"\*?\}.*?\\end\{" + env + r"\*?\}", " ", s, flags=re.S)
    s = re.sub(r"\\texttt\{[^}]*\}", " ", s)
    return [p.strip().lower() for p in re.split(r"\n\s*\n", s) if len(p.strip()) > 120]


def _pairs(text, words, noun_key):
    """{concepto: {valores mencionados}} en un trozo de texto."""
    usables = {w: v for w, v in words.items() if v >= MIN_VAL}
    alt = "|".join(sorted(usables, key=len, reverse=True))
    out = {}
    for row in NOUNS:
        concept, noun = row[0], row[noun_key]
        # (?<![\w-]) mata los compuestos: "twenty-three multiplets" no cuenta como "three",
        # porque en castellano "veintitres" va en una palabra y no habria con que emparejarlo.
        # hasta DOS palabras cortas entre el numeral y el sustantivo: el ingles pone el
        # adjetivo delante ("seventy odd rungs") y el castellano detras ("setenta peldanyos
        # impares"), y sin este hueco la misma frase cuenta en una edicion y no en la otra.
        pat = (r"(?<![\w-])(" + alt + r")[\s-]+(?:[a-z]{1,9}[\s-]+){0,2}(?:" + noun + r")\b")
        for m in re.finditer(pat, text):
            out.setdefault(concept, set()).add(usables[m.group(1)])
    return out


def _falsify_n6():
    """La frase REAL que se colo el 24-ago-2026, .la ve N6?

    La S13 inglesa decia "reached by two multiplets" donde la castellana decia "tres
    multipletes".  Si esto deja de saltar, N6 ha vuelto a ser un adorno.
    """
    en = ("the protection needs six dimensions: on $s^1/\\zz_2$ every rung is even and $d=0$ "
          "is reached by two multiplets, while here the rungs below the fourth moment are odd.")
    es = ("la proteccion necesita seis dimensiones: en $s^1/\\zz_2$ todos los peldanos son pares "
          "y $d=0$ se alcanza con tres multipletes, mientras que aqui los peldanos por debajo "
          "del cuarto momento son impares.")
    ca, cb = _pairs(en, WORDS_EN, 1), _pairs(es, WORDS_ES, 2)
    ok = ca.get("multiplete") == {2} and cb.get("multiplete") == {3}
    print("  falsacion: la frase de la S13 del 24-ago %s"
          % ("salta, EN {2} contra ES {3}" if ok
             else "NO SALTA -- N6 esta rota (EN %s, ES %s)"
                  % (ca.get("multiplete"), cb.get("multiplete"))))
    return 0 if ok else 1


def spelled(texs):
    """N6 -- mismo parrafo, mismo objeto, .mismo valor?

    Comparar CUENTAS a nivel de documento no sirve, y se probo: el castellano escribe "seis
    dimensiones" donde el ingles escribe "six-dimensional", y "setenta peldanyos" donde el
    ingles imprime "$70$", de modo que las cuentas difieren en todas partes por idioma y en
    ninguna por error.  Aquella version devolvia nueve diferencias y las nueve eran ruido.

    Lo que NO puede diferir por idioma es el VALOR: si el parrafo ingles dice DOS multipletes
    y su gemelo castellano dice TRES, uno de los dos esta mal.  Eso es lo que se compara.
    """
    print("\nN6  los numeros con letra: mismo parrafo, mismo objeto, .mismo valor?")
    pe, ps = _paras(texs[0]), _paras(texs[1])
    if len(pe) != len(ps):
        print("  *** %d parrafos en la inglesa y %d en la castellana: no emparejan, y esa"
              % (len(pe), len(ps)))
        print("      diferencia es en si misma el hallazgo.  N6 no puede correr. ***")
        return ["letra:parrafos:%d-%d" % (len(pe), len(ps))]
    out, n = [], 0
    for i, (a, b) in enumerate(zip(pe, ps)):
        ca, cb = _pairs(a, WORDS_EN, 1), _pairs(b, WORDS_ES, 2)
        for concept in sorted(set(ca) & set(cb)):
            n += 1
            if ca[concept] != cb[concept] and not any(
                    c == concept and marca in b for c, marca in N6_ACEPTADOS):
                print("     parrafo %-3d %-12s EN %-12s ES %s"
                      % (i, concept, sorted(ca[concept]), sorted(cb[concept])))
                print("        EN: %s" % " ".join(a.split())[:86])
                print("        ES: %s" % " ".join(b.split())[:86])
                out.append("letra:%d:%s" % (i, concept))
    print("  %d parrafos emparejados, %d objetos contados en las dos, %d con valor distinto"
          % (len(pe), n, len(out)))
    if _falsify_n6():
        out.append("letra:falsacion")
    return out


def main():
    fails = []
    for path in TEXS:
        s = open(path, encoding="utf-8").read()
        abst, intro, body = sections(s)
        print("\n%s" % path)
        report("N1 abstract", abst, body, fails)
        report("N2 intro", intro, s[:s.index(r"\label{sec:intro}")] + s[s.index(r"\label{sec:setting}"):], fails)

        labs = set(LAB.findall(s))
        dangling = sorted((set(REF.findall(abst)) | set(REF.findall(intro))) - labs)
        print("  %-14s %s" % ("N3 sections",
                              "every section the summaries point at exists" if not dangling
                              else "*** dangling: %s ***" % ", ".join(dangling)))
        fails += ["%s:%s" % (path, d) for d in dangling]

        # N4 -- ETIQUETAS DUPLICADAS.  LaTeX solo AVISA de esto: compila, no falla, y \eqref
        # apunta silenciosamente a la ultima definicion.  Reutilizar eq:comb para dos ecuaciones
        # distintas hizo que una referencia de la seccion 9 apuntase a una ecuacion de la 10 sin
        # que nada se quejase.  Un aviso que nadie lee no es un aviso.
        # [[an-overloaded-symbol-becomes-a-false-claim]]
        body = s.split(r"\begin{thebibliography}")[0]
        seen = {}
        for m in re.finditer(r"\\label\{([^}]+)\}", body):
            seen.setdefault(m.group(1), []).append(m.start())
        dup = sorted(k for k, v in seen.items() if len(v) > 1)
        print("  %-14s %s" % ("N4 labels",
                              "%d labels, none defined twice" % len(seen) if not dup
                              else "*** defined twice: %s ***" % ", ".join(dup)))
        fails += ["%s:dup:%s" % (path, d) for d in dup]

        # y lo mismo desde el otro lado: lo que el propio pdflatex reporto en su .log
        log = pathlib.Path(path.replace(".tex", ".log"))
        if not log.exists():
            print("  %-14s *** no hay .log: compila antes, o esta comprobacion no existe ***"
                  % "N4b build log")
            fails.append("%s:log:missing" % path)
        else:
            txt = log.read_text(encoding="utf-8", errors="replace")
            bad = [w for w in ("multiply defined",
                               "There were undefined references",
                               "Citation `",              # pdflatex escribe esto solo si falta
                               "Reference `") if w in txt]
            print("  %-14s %s" % ("N4b build log",
                                  "clean" if not bad else "*** %s ***"
                                  % ", ".join(w.strip(" `") for w in bad)))
            fails += ["%s:log:%s" % (path, w.strip(" `")) for w in bad]

    # ---- N5: LAS DOS EDICIONES TIENEN QUE CITAR LOS MISMOS NUMEROS -------------------------
    # Ninguna otra compuerta compara una edicion con la otra.  check_parity cuenta estructuras --
    # secciones, figuras, entornos tabular -- y check_numbers comprueba archivos; si el ingles
    # dice 10.86 y el espanol se quedo en 11.2, las dos pasan.  La primera vez que esto corrio
    # encontro cuatro cosas: la tabla de criterios con cuatro columnas solo en ingles, el parrafo
    # entero sobre la procedencia de la matriz de HEPData ausente en espanol, una fila del libro
    # de cuentas una version por detras, y -- la peor -- el espanol llamando todavia "el punto del
    # 95 %" a un umbral que el ingles ya habia dejado de llamar asi.
    # [[translating-finds-what-review-missed]]
    print("\nN5  las dos ediciones, numero a numero")
    cnt = {}
    for path in TEXS:
        s = open(path, encoding="utf-8").read().split(r"\begin{thebibliography}")[0]
        s = re.sub(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", " ", s, flags=re.S)
        s = re.sub(r"(?<!\\)%.*", "", s)          # los comentarios no los lee nadie
        c = {}
        for t in NUM.findall(s):
            c[t] = c.get(t, 0) + 1
        cnt[path] = c
    en, es = cnt[TEXS[0]], cnt[TEXS[1]]
    keys = sorted(set(en) | set(es))
    # solo los numeros con punto decimal: los enteros pequenos entran en giros de lengua
    # ("a factor of three") y compararlos produciria ruido en vez de senal.
    diff = [k for k in keys if "." in k and en.get(k, 0) != es.get(k, 0)]
    print("  %d numeros decimales distintos; %d con conteo distinto entre ediciones"
          % (len([k for k in keys if "." in k]), len(diff)))
    for k in diff:
        print("     %-10s EN %2d   ES %2d" % (k, en.get(k, 0), es.get(k, 0)))
        fails.append("edicion:%s" % k)

    fails += spelled(TEXS)

    print("")
    if fails:
        print("PROBLEMAS: %d" % len(fails))
        return 1
    print("OK -- the abstract and the introduction say what the body says,")
    print("     and the two editions quote the same numbers, in digits and in words.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
