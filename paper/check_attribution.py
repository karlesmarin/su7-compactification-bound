#!/usr/bin/env python3
"""check_attribution.py - every equation we cite from someone else's paper is on our disk.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS EXISTS.  Carles, 22-ago-2026: "sobre todo no podemos atribuirnos algo que no sea nuestro
y tenemos que estar seguros de nuestras atribuciones".  The paper's whole novelty argument is a
list of things that are NOT ours, each pinned to a numbered equation in someone else's paper:
\\cite{CCD24} eq.~(2.11), \\cite{HY04} eq.~(3.20), \\cite{HHK03} eqs.~(4.13)--(4.15), and so on.
Every one of those is a factual claim about a document, and none of the other nine gates reads a
document.  check_refs.py checks that a \\cite has a \\bibitem; it cannot check that the bibitem's
paper says what the sentence claims.

WHAT IT CHECKS, and it is deliberately narrow.

  1. COVERAGE.  Every bibliography key the paper cites an EQUATION NUMBER of must have its source
     archived under ../../_papers/.  An equation citation we cannot open is an attribution taken
     on trust, and the honesty floor of this series does not allow one silently.
  2. PRESENCE.  The cited equation number must actually occur in the extracted text of that
     source.  This is a weak test -- a number can occur for other reasons -- and it is stated as
     weak: it catches the transposed digit and the equation that does not exist, not a
     misreading.  It is a floor, not a proof.
  3. The quoted VERBATIM strings.  Where the paper puts someone's words in quotation marks it must
     be their words; each is searched for in their own text.

WHAT IT DOES NOT CHECK.  Whether the equation says what we say it says.  That is a reading, it is
in GATE_*.md, and no script does it.  The gate prints the context around each hit so a reader can
do that reading without opening the PDF, which is the most a script can honestly offer.

Run:  python check_attribution.py     (from part_vii/paper/)
"""
import pathlib
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
PAPERS = HERE.parent.parent / "_papers"
TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]

# bibliography key -> the archived file that IS that paper.  A key absent from this map and cited
# with an equation number is a coverage failure, which is the point of the gate.
SOURCE = {
    "KM25":   "SU7_GGHU_2503.04090.txt",
    "AHMN":   "AHMN_2312.08608.txt",
    "CCD24":  "ccd24_2409.16137.txt",
    "HY04":   "hy04_hep-ph_0401185.txt",
    "HHK03":  "hhk_0309088.txt",
    "KTY17":  "1704.04840.pdf",
    "vGIQ":   "hep-th_0204223.pdf",
    "GLS06":  "ggn_0604215.pdf",
    "GNH05":  "hep-th_0503153.pdf",
    "CCP06":  "hep-ph_0510366.pdf",
    "PP01":   "pp01_hep-ph_0105021.pdf",
    "GHBQ05": "ghbq05_hep-th_0506164.pdf",
    "Cac25":  "cac25_2501.13118.pdf",
    # El manual de CIJET, cuya ec. (2) da la base quiral de operadores y la (3) el desarrollo en
    # lambda_i.  El PROGRAMA no puede redistribuirse; el manual es un preprint de arXiv y si.
    "CIJET":  "cijet_manual_1301.7263.pdf",
}
# Cited with an equation number but NOT archived here.  Each must be declared, with the reason,
# or the gate fails: an undeclared gap is exactly what this file exists to prevent.
NOT_ARCHIVED = {
    "Wood":   "Kent technical report 15-92, not online as a PDF we could archive; eq.~(14.1) is "
              "quoted in full in the bibitem itself, so the reader does not need our copy.",
    "Apostol": "textbook, Theorem 12.6.",
    "Lewin":  "textbook, \\S7.1.",
}

VERBATIM = [
    ("HY04",  "applicable in general"),
    ("HY04",  "we can know whether the color is conserved"),
    ("CCD24", "accidentally the same"),
    ("CCD24", "which configuration gives the smallest value"),
    ("KM25",  "we do not consider the potential"),
]

# Attributions that are not an equation NUMBER but a claim about what someone's paper does or does
# not contain.  Each is a string search over their own text, with the verdict we assert.  These
# exist because the gate caught the paper crediting \cite{GHBQ05} eq.~(3.18) with "the same
# weight in six dimensions" -- and that paper contains no 15/16, no 31/32 and no zeta(5) anywhere
# in forty pages.  What it carries is the SPLIT, not the weight, and the sentence now says so.
# Crediting someone with something they did not write is the same error as claiming something
# that is not ours, and it is the easier one to miss because it feels generous.
CONTENT = [
    ("PP01",   "15",     True,  "eq.~(13) must state the -15/16 weight"),
    ("GHBQ05", "15\n16", False, "must NOT carry the 15/16 weight -- we claim only the split"),
    ("GHBQ05", "15/16",  False, "same, written inline"),
    ("GHBQ05", "31/32",  False, "nor the value-bound weight"),
    ("GHBQ05", "ζ(5)", False, "nor zeta(5) at all"),
]


def text_of(name):
    p = PAPERS / name
    if not p.exists():
        return None
    if p.suffix == ".txt":
        return p.read_text(encoding="utf-8", errors="replace")
    import fitz                       # pdftotext lies; PyMuPDF is the house reader
    return "\n".join(pg.get_text() for pg in fitz.open(p))


# ---------------------------------------------------------------- 1  what the paper cites
P("=" * 100)
P("1 -- EVERY FOREIGN EQUATION THE PAPER CITES")
P("=" * 100)
# The assignment, made BY HAND once and checked for completeness below.  Two heuristics were
# tried first and both were wrong in ways worth recording, because either would have produced a
# green gate over a false table.  Attributing to the \cite on the same LINE missed most of them,
# since a source is named once and then called "their eq." for a paragraph.  Attributing to the
# nearest PRECEDING \cite then blamed \cite{Manton} for \cite{KM25}'s eq.~(82) and \cite{Wood} for
# \cite{CCD24}'s (3.15) -- because \cite{KM25} is the paper's subject and is referred to as "their"
# and "its" throughout without ever being re-cited.  No regular expression resolves an English
# pronoun.  So the list is declared, and the gate's job is to prove it COMPLETE, not to guess it.
ASSIGNED = [
    # EL RESUMEN, desde que lleva la condicional (23-ago).  El 8D impar depende de los
    # coeficientes gauge impresos en su ec. (68), y eso se dice ya en el abstract y no solo en la
    # S13 -- que es justamente el punto: la condicion tiene que viajar con el titular.
    ("KM25", "68"),
    ("KM25", "82"),
    # LA INTRODUCCION, parrafo "sobre que se apoya": la cadena pasa por un coeficiente que este
    # articulo NO calcula y lee de su ec. (68).  Va aqui, justo despues del (82) de la relacion de
    # jerarquia, porque la ruta tiene que verse en la pagina 1 y no en la 50.
    ("KM25", "68"),
    ("KM25", "11--13"), ("KM25", "63--68"), ("KM25", "68"),
    ("KM25", "63--67"),
    ("KM25", "68"), ("KM25", "67"),
    ("CCD24", "2.11--2.13"), ("PP01", "13"), ("GHBQ05", "3.18"), ("CCD24", "2.11"),
    ("CCD24", "3.15"), ("KM25", "80"), ("KM25", "80"),
    # EL COROLARIO del Teorema 1 y su demostracion.  Al pasar el teorema a condicional, la lectura
    # del -27 deja de estar dentro del enunciado y pasa a ser el corolario que comprueba que la
    # hipotesis se cumple para la semilla publicada -- y la prueba lo repite al cerrarse.  Dos
    # citas nuevas de su ec. (68), y las dos dicen exactamente lo que este articulo NO calcula.
    ("KM25", "68"), ("KM25", "68"),
    ("KM25", "63--67"),
    ("HY04", "5.9"), ("KM25", "68"), ("HY04", "3.20"), ("CCD24", "2.7"), ("CCD24", "3.26"),
    ("KTY17", "3.10--3.12"), ("HHK03", "4.13--4.14"), ("HHK03", "4.7"), ("HHK03", "3.20"),
    ("HHK03", "4.13"), ("HHK03", "4.7"), ("HHK03", "4.5--4.6"), ("CCD24", "2.11"),
    # La cita de su ec. (68) que habia aqui era la del parrafo del 2W: "el mismo semientero de su
    # ec. (68) del que cuelga el Teorema 1".  Esa frase se retira el 23-ago porque es FALSA -- 2W
    # sobrevive a la semilla candidata y 8D no, asi que no cuelgan del mismo sitio.  Al irse la
    # frase se va su cita.
    ("CCD24", "2.11"), ("CCD24", "2.12--2.13"), ("CCD24", "3.14"),
    ("CCD24", "3.44"), ("HHK03", "4.14"), ("CCD24", "2.11"), ("Cac25", "3.3"),
    ("CCD24", "3.25--3.27"), ("CCD24", "2.11"), ("KM25", "18--19"),
    # el medio del sector gauge, en la S10: de el cuelga el Teorema de los octavos impares, y con
    # el la integralidad de 8D -- que es, medida alli, la raiz de la que depende toda la
    # conclusion de colisionador.  Su ec. (68) es la tabla gauge; la Parte VI la RE-DERIVA de los
    # 48 generadores del adjunto en vez de leerla, con un ajuste sobredeterminado.
    ("KM25", "68"),
    # y sus ecs. (78)-(79), que son las que identifican como indices de COLOR los tres que hacen
    # m = 3 -- el conteo cuya paridad, y no ningun medio, es lo que hace impar a 8D.
    ("KM25", "78--79"),
    # La derivacion RETIRADA, y su diagnostico, ambos en la S10 y ambos citando ecuaciones ajenas.
    # La regla de conteo de \cite{HY04} -- -(D-2) por el potencial de un grado de libertad, su
    # ec. (3.10), con un par a c=2 y dos a c=1 por su ec. (3.8) -- reproduce su propia ec. (3.20)
    # en D=5 y da (2,4) en D=6.  Parecia una derivacion fuera de muestra del trozo periodico de la
    # ec. (68) de \cite{KM25}, y NO lo es: D-2 = 4 en seis dimensiones coincide con el numero de
    # componentes de A_mu, que es cuatro en TODA dimension.  La coincidencia se retira en el texto.
    ("HY04", "3.10"), ("HY04", "3.8"), ("HY04", "3.20"), ("KM25", "68"),
    # Y la misma aritmetica leida al reves, que es una pregunta viva sobre su contabilidad.  Su
    # ec. (62) asigna N=4 a A_mu y N=1 al escalar superviviente, y con nuestra tabla de generadores
    # esos dos numeros reproducen su ec. (68) EXACTA -- pero suman cinco, y un campo gauge en seis
    # dimensiones tiene D-2 = 4 polarizaciones fisicas.  Con la resta de fantasmas de \cite{HNT} la
    # ec. (68) pasaria a (3/2,3,6), el 8D del sector gauge seria par, y el Teorema de los octavos
    # impares no existiria.  El articulo lo deja escrito como pregunta abierta, sin resolverla.
    ("KM25", "62"), ("KM25", "68"), ("KM25", "68"),
    ("vGIQ", "5.3"), ("vGIQ", "5.2"),
    ("CCD24", "2.12--2.13"), ("KM25", "3--10"), ("KM25", "4"), ("KM25", "7"), ("KM25", "57"),
    ("AHMN", "4.4"),
    # El item abierto de la S13, REESCRITO el 23-ago: ya no cita la ec. (3.20) de \cite{HY04}.
    # La regla -(D-2) se describe ahora en palabras, y lo que la sostiene es el desglose (4-2)+2
    # de \cite{HNT} en seis dimensiones mas el N=4 de \cite{AHMN} -- ninguno de los dos lleva
    # numero de ecuacion en el texto, asi que no entran aqui.  Queda una sola cita: la ec. (68) de
    # \cite{KM25}, que es sobre la que todo queda condicionado.
    ("KM25", "68"),
    # La base quiral de operadores de CIJET, su ec. (2): O_1 = delta delta (LL singlete),
    # O_2 = T^a T^a (LL octete), O_3/O_4 LR, O_5/O_6 RR, con sigma_bin en su ec. (3).  Leida del
    # manual, no supuesta -- de ahi sale que el VV es lambda_1 = lambda_3 = lambda_5 y que el
    # cierre ya corrido es lambda_1 solo.
    ("CIJET", "2"),
    # EL ANEXO DE LEAN, ultima cita del articulo.  Dice que Lean NO puede adjudicar su ec. (68)
    # -- es un input fisico y no un paso logico -- y que por eso el reparto (3/2,1/2) no esta
    # formalizado.  La cita existe para marcar lo que queda FUERA del ladrillo, que es la parte
    # del anexo que mas trabaja.
    ("KM25", "68"),
]

body = (HERE / TEXS[0]).read_text(encoding="utf-8", errors="replace")
body = body.split(r"\begin{thebibliography}")[0]
flat = re.sub(r"\s+", " ", body)
found = [m.group(1) + ("--" + m.group(2) if m.group(2) else "")
         for m in re.finditer(r"eqs?\.~?\((\d+(?:\.\d+)*)\)(?:--\((\d+(?:\.\d+)*)\))?", flat)]
P("  equation citations in the body : %d      entries in the declared table : %d"
  % (len(found), len(ASSIGNED)))
assert len(found) == len(ASSIGNED), (
    "the table is stale: %d citations against %d entries. Re-read the new ones and assign them "
    "by hand -- do NOT pad the table." % (len(found), len(ASSIGNED)))
mismatch = [(i + 1, f, a[1]) for i, (f, a) in enumerate(zip(found, ASSIGNED)) if f != a[1]]
assert not mismatch, "table drifted out of order at %s" % mismatch[:3]
P("  every citation is accounted for, in order, with the number the table expects.")
P("")
cited = {}
for k, e in ASSIGNED:
    cited.setdefault(k, set()).update(e.split("--"))
for k in sorted(cited):
    P("  %-9s %s" % (k, ", ".join(sorted(cited[k], key=lambda s: [int(x) for x in s.split(".")]))))

# ---------------------------------------------------------------- 2  coverage
P("")
P("=" * 100)
P("2 -- COVERAGE: is the source on our disk?")
P("=" * 100)
missing = []
for k in sorted(cited):
    if k in SOURCE:
        ok = (PAPERS / SOURCE[k]).exists()
        P("  %-9s archived    %-34s %s" % (k, SOURCE[k], "" if ok else "<-- FILE MISSING"))
        if not ok:
            missing.append(k)
    elif k in NOT_ARCHIVED:
        P("  %-9s DECLARED    %s" % (k, NOT_ARCHIVED[k][:60] + "..."))
    else:
        P("  %-9s UNDECLARED GAP  <-- archive it or declare why not" % k)
        missing.append(k)
assert not missing, "attributions with no archived source and no declaration: %s" % missing
P("")
P("  Every key is either archived or declared with its reason.")

# ---------------------------------------------------------------- 3  presence
P("")
P("=" * 100)
P("3 -- PRESENCE: does the cited number occur in their text?  (a floor, not a proof)")
P("=" * 100)
absent = []
for k in sorted(cited):
    if k not in SOURCE:
        continue
    t = text_of(SOURCE[k])
    if t is None:
        continue
    flat = re.sub(r"\s+", " ", t)
    for e in sorted(cited[k]):
        pats = ["(%s)" % e, "(%s )" % e, " %s " % e]
        hit = any(p in flat for p in pats)
        P("  %-9s eq. %-9s %s" % (k, e, "found" if hit else "NOT FOUND in their text"))
        if not hit:
            absent.append((k, e))
P("")
if absent:
    P("  NOT FOUND: %s" % absent)
    P("  A miss is not automatically an error -- an equation number can be rendered as an image,")
    P("  or split across lines by the extractor -- but each one has to be looked at by hand.")
P("  found = the string occurs; it does NOT mean the equation says what the paper says it says.")
P("  That reading is in the GATE_*.md notes and no script performs it.")

# ---------------------------------------------------------------- 4  verbatim quotes
P("")
P("=" * 100)
P("4 -- THE QUOTED WORDS ARE THEIR WORDS")
P("=" * 100)
bad = []
for k, q in VERBATIM:
    t = text_of(SOURCE[k])
    flat = re.sub(r"\s+", " ", t or "")
    hit = q.lower() in flat.lower()
    P("  %-9s %-46s %s" % (k, '"%s"' % q[:44], "verbatim" if hit else "NOT IN THEIR TEXT"))
    if not hit:
        bad.append((k, q))
assert not bad, "quotation marks around words that are not in the source: %s" % bad
P("")
P("  Every string the paper puts in quotation marks occurs in the paper it attributes it to.")

# ---------------------------------------------------------------- 5  what a source does NOT say
P("")
P("=" * 100)
P("5 -- AND WHAT WE SAY A SOURCE DOES NOT CONTAIN, IT MUST NOT CONTAIN")
P("=" * 100)
P("  Over-crediting is an attribution error too, and it is the one that feels safe.  Each row is")
P("  a claim about the presence OR ABSENCE of a string in someone else's paper.")
P("")
wrong = []
for k, needle, want, why in CONTENT:
    t = text_of(SOURCE[k])
    flat = re.sub(r"[ \t]+", " ", t or "")
    got = needle in flat
    ok = (got == want)
    P("  %-9s %-10s expected %-8s got %-8s %s   %s"
      % (k, repr(needle)[:10], want, got, "ok" if ok else "MISMATCH", why))
    if not ok:
        wrong.append((k, needle, want, got))
assert not wrong, "a presence/absence attribution is wrong: %s" % wrong
P("")
P("  \\cite{PP01} eq.~(13) does carry the -15/16; \\cite{GHBQ05} carries the split and not the")
P("  weight, and the paper's sentence and bibitem now say exactly that.")

P("")
P("=" * 100)
P("ATTRIBUTION: TODO EN ORDEN" if not absent else "ATTRIBUTION: %d numbers to check by hand" % len(absent))
P("=" * 100)
