#!/usr/bin/env python3
"""check_prose.py - decima compuerta: una frase que se quedo huerfana al mover un parrafo.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

La edicion castellana imprimia, en Data availability:

    ... un verificador lo impone y ha saltado dos veces. Los artefactos principales son
    Un instrumento interactivo de esta serie corre en karlesmarin.github.io/ghu-explorer ...

El parrafo del explorador se inserto ENTRE la entradilla y su lista, y la entradilla se quedo
colgando delante de el.  Ocho compuertas verdes y una lectura completa no lo vieron: ninguna
lee prosa.  Lo vio la comparacion linea a linea de los dos fuentes.

La firma es estrecha y por eso sirve de compuerta: un verbo copulativo seguido de una palabra
capitalizada SIN punto en medio.  Sobre las dos ediciones da tres aciertos legitimos, todos
nombres propios (`is Sym2`, `is Weinberg`, `es Sym2`), asi que la lista blanca es corta y
explicita -- si crece, hay que mirar por que, no ampliarla a ciegas.

  python check_prose.py             comprueba las dos ediciones
  python check_prose.py --falsify   vuelve a meter el fallo real y comprueba que salta
"""
import re
import sys
import pathlib

import fitz

# The snippets this gate prints come out of the PDF's own text layer, so they carry whatever
# glyphs the paper does -- arrows, primes, Greek.  On Windows the default console encoding is
# cp1252 and print() raised UnicodeEncodeError *while reporting a finding*: the gate died with a
# traceback and a non-zero exit that looked like a crash rather than like a hit.  A gate that
# cannot say what it found is a gate that can be misread as broken.  [[edits-must-fail-loudly]]
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except Exception:
    pass

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent
PDFS = ["su7_hierarchy.pdf", "su7_hierarchy_es.pdf"]

COPULA = {"son", "es", "are", "is", "fue", "era", "eran", "were", "was"}
# nombres propios que SI van detras de una copula, uno a uno y con su motivo
ALLOWED = {
    # uno a uno y con su motivo; si esta lista crece hay que mirar POR QUE, no ampliarla
    "Part": "D is Part VI's curvature...; its role here is Part IV's -- el nombre de una parte",
    "Fermat": "the mod-six law follows from c^4 = c^2 (mod 3), which is Fermat -- un teorema",
    "Sk": "the rung-k coefficient is Sk + (2^p-1)Dk / 2^p -- un simbolo compuesto en el PDF",
    "Sl": "stationarity is Sl4(theta)+Sl4(2theta)=0 -- la funcion de Clausen, un simbolo",
    "Theorem": "which is Theorem 1 / Theorem 2 -- referencia a un teorema numerado",
    # entro el 24-ago con el pie de la Figura 6, al marcar alli la semilla: es la MISMA clase
    # que "Theorem", una referencia numerada que el PDF imprime capitalizada.  No se amplia a
    # ciegas -- se anyade porque el motivo es el de la linea de arriba, palabra por palabra.
    "Corollary": "which is Corollary 1 -- referencia a un corolario numerado, como Theorem",
    "Representation": "el TITULO de la Parte III es '...Is Representation-Dependent'",
    "Euler": "the resummation is Euler's identity -- el apellido de quien la escribio",
}
# la cola tiene que ser {1,}: el fallo real era "son Un instrumento", y con {2,} la compuerta
# no lo veia.  Se descubrio con --falsify, que es exactamente para lo que esta.
PAT = re.compile(r"\b([a-zA-Zaeiouñüà-ÿ]{2,6}) ([A-ZÀ-Ý][a-zà-ÿ]{1,})")


def scan(text):
    text = re.sub(r"-\n", "", text)
    text = re.sub(r"\s+", " ", text)
    out = []
    for m in PAT.finditer(text):
        if m.group(1).lower() not in COPULA:
            continue
        out.append((m.group(2), text[max(0, m.start() - 60):m.end() + 40]))
    return out


def main():
    falsify = "--falsify" in sys.argv
    bad = 0
    for f in PDFS:
        d = fitz.open(HERE / f)
        text = "\n".join(d.load_page(i).get_text() for i in range(d.page_count))
        hits = scan(text)
        unknown = [(w, c) for w, c in hits if w not in ALLOWED]
        P("  %-24s %d copula+mayuscula, %d en la lista blanca, %d sin explicar"
          % (f, len(hits), len(hits) - len(unknown), len(unknown)))
        for w, c in unknown:
            P("       *** %s ***  ...%s" % (w, c))
        bad += len(unknown)

    if falsify:
        P("")
        P("  --falsify: se reinyecta el fallo real que motivo esta compuerta.")
        broken = ("un verificador lo impone y ha saltado dos veces. Los artefactos "
                  "principales son Un instrumento interactivo de esta serie corre en")
        hits = [(w, c) for w, c in scan(broken) if w not in ALLOWED]
        P("     %s" % ("DETECTADO: %s" % hits[0][0] if hits else "*** NO DETECTADO ***"))
        if not hits:
            return 1

    P("")
    if bad:
        P("PROSA HUERFANA: %d sitio(s) donde una frase entra en otra sin punto." % bad)
        return 1
    P("prosa: ninguna entradilla colgando en ninguna de las dos ediciones")
    return 0


if __name__ == "__main__":
    sys.exit(main())
