#!/usr/bin/env python3
"""check_parity.py - las dos ediciones de la Parte VII dicen estructuralmente lo mismo.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

Adaptado del verificador homonimo de part_vi/paper. La convencion de la casa es
traduccion COMPLETA, no resumen. Pero una edicion se edita despues de la otra y
nada avisa si un bloque se queda sin traducir: un parrafo que existe en una
edicion y no en la otra no rompe ningun build, y tampoco lo rompe una fila de
tabla de mas o de menos --- eso ultimo lo aprendimos en la Parte VI, con 37
filas en ingles y 38 en castellano.

Cuenta los elementos estructurales de cada edicion, compara etiquetas y citas, y
compara las filas de cada longtable una a una. Ademas comprueba lo especifico de
esta traduccion: que la edicion castellana apunte a las figuras _es y que no
quede ningun "\\,\\%" en ella (babel-spanish lo rechaza en modo matematico).

Uso:  python check_parity.py     (desde part_vii/paper/)
"""
import re
import sys

PAR = [("secciones", r"\\section\{"), ("puentes", r"\\bridge\{"),
       ("keyeq", r"\\begin\{keyeq\}"), ("figuras", r"\\begin\{figure\}"),
       ("theorem", r"\\begin\{theorem\}"), ("proof", r"\\begin\{proof\}"),
       ("paragraph", r"\\paragraph\{"), ("bibitem", r"\\bibitem\{"),
       ("longtable", r"\\begin\{longtable\}"), ("tabular", r"\\begin\{tabular\}"),
       ("itemize", r"\\begin\{itemize\}"), ("item", r"\\item\b")]


def carga(p):
    s = open(p, encoding="utf-8").read()
    return re.sub(r"(?<!\\)%[^\n]*", "", s)


en, es = carga("su7_hierarchy.tex"), carga("su7_hierarchy_es.tex")

print("  %-14s %-6s %-6s %s" % ("elemento", "EN", "ES", ""))
malos = []
for nombre, pat in PAR:
    a, b = len(re.findall(pat, en)), len(re.findall(pat, es))
    if a != b:
        malos.append("%s: EN %d, ES %d" % (nombre, a, b))
    print("  %-14s %-6d %-6d %s" % (nombre, a, b, "" if a == b else "*** DESIGUAL ***"))

for nombre, pat in (("label", r"\\label\{([^}]*)\}"), ("cite", r"\\cite\{([^}]*)\}"),
                    ("eqref", r"\\eqref\{([^}]*)\}")):
    A, B = set(re.findall(pat, en)), set(re.findall(pat, es))
    if A != B:
        malos.append("%s solo en EN: %s | solo en ES: %s" % (nombre, sorted(A - B), sorted(B - A)))
    print("  %-14s %-6d %-6d %s" % (nombre + "s", len(A), len(B),
                                    "" if A == B else "*** DISTINTAS ***"))


def filas_longtable(s):
    out, pos = [], 0
    while True:
        i = s.find("\\begin{longtable}", pos)
        if i < 0:
            return out
        j = s.find("\\end{longtable}", i)
        cuerpo = s[i:j]
        out.append(len([r for r in cuerpo.split("\\\\") if "&" in r]))
        pos = j + 1


fa, fb = filas_longtable(en), filas_longtable(es)
print()
print("  %-14s %-6s %-6s %s" % ("longtable", "EN", "ES", "filas de cada una"))
for n, (a, b) in enumerate(zip(fa, fb)):
    if a != b:
        malos.append("filas de la longtable %d: EN %d, ES %d" % (n, a, b))
    print("  %-14s %-6d %-6d %s" % ("  tabla %d" % n, a, b, "" if a == b else "*** DESIGUAL ***"))

# --- lo especifico de esta traduccion.
figs_en = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]*)\}", en)
figs_es = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]*)\}", es)
esperado = [f.replace(".pdf", "_es.pdf") for f in figs_en]
print()
print("  figuras: EN %d, ES %d" % (len(figs_en), len(figs_es)))
if figs_es != esperado:
    malos.append("la edicion ES no apunta a las figuras _es: %s" % figs_es)
else:
    print("  todas las figuras de la edicion ES son las variantes _es")

n_pct = len(re.findall(r"\\,\\%", es))
if n_pct:
    malos.append("quedan %d '\\,\\%%' en la edicion ES (babel-spanish los rechaza)" % n_pct)
else:
    print("  ningun '\\,\\%' en la edicion ES")

# --- el cuerpo esta realmente traducido: ninguna palabra funcional inglesa suelta.
cuerpo = es[es.find("\\begin{document}"):es.find("\\begin{thebibliography}")]
cuerpo = re.sub(r"\\bibitem\{[^}]*\}", "", cuerpo)
# UNA CITA TEXTUAL NO ES UN RESIDUO DE TRADUCCION.  El articulo cita el registro de HEPData con
# sus propias palabras -- "the correlation matrix of the maximum likelihood estimators..." -- y
# traducir una cita documental para contentar a una compuerta seria falsearla: lo que se afirma
# es precisamente que el registro dice ESO.  Asi que lo entrecomillado con ``...'' se descuenta,
# igual que ya se descuentan los \bibitem.  [[quote-the-number-from-the-aux]]
n_citas = len(re.findall(r"``.*?''", cuerpo, flags=re.S))
cuerpo = re.sub(r"``.*?''", " ", cuerpo, flags=re.S)
print("  %d citas textuales descontadas antes de buscar palabras inglesas" % n_citas)
ingles = re.findall(r"(?<![\\\w])(the|which|whose|there|these|those|from|that|when|where|what)"
                    r"(?![\w])", cuerpo)
if ingles:
    malos.append("quedan %d palabras funcionales inglesas en el cuerpo ES: %s"
                 % (len(ingles), sorted(set(ingles))))
else:
    print("  el cuerpo de la edicion ES no lleva palabras funcionales inglesas")

print()
if malos:
    print("PROBLEMAS:")
    for m in malos:
        print("   " + m)
else:
    print("paridad estructural: las dos ediciones llevan los mismos elementos")
sys.exit(1 if malos else 0)
