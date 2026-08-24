#!/usr/bin/env python3
"""check_keyeq.py - las cajas enmarcadas dicen lo mismo en las dos ediciones.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

POR QUE EXISTE.  El entorno `keyeq` esta definido en el preambulo como "a framed display for the
statements the paper stands on -- and only those".  Es lo primero que mira quien hojea, y por eso
es el peor sitio posible para una frase rancia.

Y ha fallado dos veces:

  24-ago-2026.  La caja de |F(0)| >= zeta(5)/32 se leia incondicional en las dos ediciones, cuando
  en la rama candidata la demostracion se cae.  Lo encontro un repaso a mano.

  24-ago-2026, mismo dia, tercer repaso.  La caja del final de la S9 castellana decia que "una
  maquina cuyo alcance en masa cubra 10-20 TeV ZANJA LA CLASE".  La inglesa habia retirado esa
  frase el 23-ago (commit d5b5f58, "cuatro afirmaciones falsas mias, TODAS A FAVOR") y la habia
  sustituido por la correcta y mas debil: un analisis con sensibilidad demostrada PODRIA PONER A
  PRUEBA la clase.  La castellana se quedo con la version fuerte durante un dia entero, dentro
  del marco.  [[no-claim-lives-in-one-place]] [[the-abstract-must-carry-the-hypothesis]]

Ninguna de las doce compuertas anteriores mira dentro de una caja como caja: `check_formulas`
compara los DISPLAYS, `check_inline` la matematica EN LINEA, `check_narrative` los NUMEROS del
resumen.  La prosa enmarcada no la miraba nadie.

COMO FUNCIONA.  Empareja las cajas de las dos ediciones por posicion y compara lo que una
traduccion NO cambia: la etiqueta de la ecuacion que llevan dentro, los numeros, las referencias
cruzadas, las citas, y el numero de enfasis.  No compara palabras: no sabe traducir.

LO QUE NO PUEDE HACER.  Si las dos ediciones dicen lo mismo de forma equivocada, pasa.  Es una
compuerta de PARIDAD, no de verdad.  Contra eso solo hay leerlas.

Uso:  python check_keyeq.py     (desde part_vii/paper/)
"""
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
EN, ES = "su7_hierarchy.tex", "su7_hierarchy_es.tex"
P = lambda *a: print(*a, flush=True)


def strip_comments(s):
    return re.sub(r"(?<!\\)%.*", "", s)


def boxes(path):
    body = strip_comments(open(path, encoding="utf-8", errors="replace").read())
    out = []
    for m in re.finditer(r"\\begin\{keyeq\}(.*?)\\end\{keyeq\}", body, re.S):
        line = body[:m.start()].count("\n") + 1
        out.append((line, m.group(1)))
    return out


def sig(b):
    """Lo que una traduccion no mueve."""
    return {
        "label": sorted(re.findall(r"\\label\{(.+?)\}", b)),
        "ref": sorted(re.findall(r"\\(?:eq)?ref\{(.+?)\}", b)),
        "cite": sorted(re.findall(r"\\cite\{(.+?)\}", b)),
        # numeros con decimales y enteros de tres cifras o mas: los que son AFIRMACIONES.
        # Los sueltos (2, 5, 8) son ruido de redaccion y de subindices.
        "num": sorted(set(re.findall(r"(?<![\w.])\d+\.\d+", b))
                      | set(re.findall(r"(?<![\w.\\{])\d{3,}(?![\d}])", b))),
        "emph": len(re.findall(r"\\emph\{", b)),
        "bf": len(re.findall(r"\\textbf\{", b)),
    }


# --------------------------------------------------------------------------- falsacion
# Las dos cajas REALES que se colaron.  Si esta compuerta deja de verlas, se ha roto.
FALSIFY = [
    ("la maquina que zanjaba la clase, 24-ago-2026",
     r"programa entero. Un analisis con sensibilidad demostrada",   # EN (sin acentos: es el ES)
     r"programa entero. Una maquina cuyo \emph{alcance en masa de resonancia} ---no cuya energia "
     r"de haz--- cubra $10$--$20$~TeV zanja la clase"),
]


def falsify():
    P("\n" + "=" * 92)
    P("FALSACION -- la caja rancia del 24-ago, .la ve?")
    P("=" * 92)
    fails = 0
    for name, a, b in FALSIFY:
        sa, sb = sig(a), sig(b)
        # la castellana rancia lleva un numero (10, 20) que la inglesa no lleva: la firma difiere
        differs = any(sa[k] != sb[k] for k in ("num", "ref", "cite", "emph"))
        fails += 0 if differs else 1
        P("  %-46s %s" % (name, "la firma difiere, salta" if differs
                          else "NO LA VE -- la compuerta esta rota"))
    if fails:
        P("\n%d caja(s) rancia(s) que esta compuerta dejaria pasar." % fails)
    else:
        P("\nla caja rancia salta.")
    return fails


def main():
    be, bs = boxes(os.path.join(HERE, EN)), boxes(os.path.join(HERE, ES))
    P("=" * 92)
    P("LAS CAJAS ENMARCADAS, .DICEN LO MISMO EN LAS DOS EDICIONES?")
    P("=" * 92)
    P("  %d cajas en la inglesa, %d en la castellana" % (len(be), len(bs)))
    bad = 0
    if len(be) != len(bs):
        P("\n  EL NUMERO DE CAJAS NO COINCIDE.  Una caja de mas o de menos es un enunciado")
        P("  de mas o de menos en el sitio donde mas se lee.")
        bad += 1
    P("")
    for i in range(min(len(be), len(bs))):
        (le, te), (ls, ts) = be[i], bs[i]
        se, ss = sig(te), sig(ts)
        diffs = {k: (se[k], ss[k]) for k in se if se[k] != ss[k]}
        lab = (se["label"] or ss["label"] or ["(sin label)"])[0]
        if not diffs:
            P("  [%2d] %-22s  EN l.%-5d ES l.%-5d  igual" % (i, lab[:22], le, ls))
            continue
        bad += 1
        P("  [%2d] %-22s  EN l.%-5d ES l.%-5d  DIFIERE" % (i, lab[:22], le, ls))
        for k, (a, b) in diffs.items():
            P("       %-6s EN %-34s ES %s" % (k, str(a)[:34], str(b)[:34]))
        P("       EN: %s" % " ".join(te.split())[:96])
        P("       ES: %s" % " ".join(ts.split())[:96])
    P("\n" + "=" * 92)
    if bad:
        P("CAJAS QUE NO DICEN LO MISMO: %d" % bad)
        P("Una caja keyeq es, por definicion del propio articulo, un enunciado sobre el que se")
        P("apoya.  Si las dos ediciones enmarcan cosas distintas, una de las dos esta rancia.")
        P("=" * 92)
        falsify()
        return 1
    P("las cajas enmarcadas llevan los mismos numeros, refs y citas en las dos ediciones.")
    P("=" * 92)
    return 1 if falsify() else 0


if __name__ == "__main__":
    raise SystemExit(main())
