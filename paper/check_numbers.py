#!/usr/bin/env python3
"""check_numbers.py - is every number printed in the Part VII paper greppable in an archived run?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Ported from part_vi/paper/check_numbers.py, with the two lessons burnt into that one kept
verbatim: keep the decimal point and match on token boundaries (a gate that matches substrings
passed `1.816` on the strength of `816`), and treat an ESCAPED percent as content, because the
naive r"%.*" silently swallows every number to the right of a printed percentage.

The Data availability section CLAIMS every displayed number regenerates from the ancillary
scripts and sits in the archived stdout beside them.  This is what makes that claim true rather
than decorative.

Run:  python check_numbers.py     (from part_vii/paper/)
"""
import os
import re

TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]
# Part VII quotes Part VI throughout -- its five rows, its 1.94/1.20 residual, the cost of the
# donation -- so the archive it is checked against has to include Part VI's runs.  A gate that
# cannot see the run behind a number reports that number as unbacked when it is not.
OUT = ["../outputs", "../../part_vi/outputs"]

ALLOW = {
    # --- their own equation numbers, quoted throughout
    "11", "12", "13", "63", "64", "65", "66", "67", "68", "79", "80", "82",
    "2.5", "2.9", "2.11", "2.13", "3.15", "5.2", "5.3", "240",
    # --- identifiers: arXiv, journals, volumes, pages, years
    "2503.04090", "2409.16137", "0204223", "1409.6539", "1305.6846",
    "1882", "1976", "2004", "2006", "2024", "27", "86", "12.6",
    # anyos del parrafo de historia de la S13: 1979 (Manton), 1983 (Hosotani) y 2007 (el
    # articulo por el que \cite{KM25} SI cita a Takenaga).  Son fechas de publicacion leidas de
    # la bibliografia, no medidas nuestras -- 1989 y 2004 ya estaban permitidos.
    "1979", "1983", "2007",
    # --- structural counts stated in words (dimensions, ranks, index labels)
    "1", "0", "2", "3", "4", "5", "6", "7", "8", "9", "10", "14",
    "28", "48", "84", "16", "31", "32", "45",
    # --- quoted from elsewhere, not measured by us
    "125", "127", "0.63", "0.653", "80.4",
    # --- numeros de ECUACION de otros articulos, uno a uno y con su fuente
    "4.4",      # AHMN, su eq. (4.4): g4 = 2 m_W / v
    "4.43",     # vGIQ, sus eqs. (4.42)-(4.43): la traza torcida
    "3.25", "3.26", "3.27",   # CCD24, sus ecs. (3.25)-(3.27): dos modelos, un solo potencial
    "2.7", "2.9",     # CCD24, sus ecs. (2.7)-(2.9): el potencial como forma lineal en el contenido
    "3.20",           # Haba-Yamashita, su ec. (3.20): la misma forma lineal, y su counting rule
    "4.13", "4.15",   # Haba-Hosotani-Kawamura, sus ecs. (4.13)-(4.15): dos funcionales enteros
    "3.10", "3.12",   # Kojima-Takenaga-Yamashita, sus ecs. (3.10)-(3.12): 21 + 28 = 48 + singlete
    "3.3",      # Cacciapaglia 2501.13118, su ec. (3.3): repite su (2.11)
    "14.1",     # Wood, Tech. Rep. 15-92, su ec. (14.1): la formula de duplicacion
    "3.18",     # Ghilencea-Hoover-Burgess-Quevedo, su ec. (3.18): el -15/16 en 6D
    "25.12",    # NIST DLMF, su seccion 25.12
    # --- anos
    "1981",     # Gross-Pisarski-Yaffe, Rev. Mod. Phys. 53 (1981) 43
    "2001",     # Ponton-Poppitz, JHEP 0106 (2001) 019 -- el -15/16, en print desde entonces
    # --- publicado por ELLOS, no medido por nosotros
    "0.71",     # vGIQ: sus dos minimos degenerados en su(3) adjunto, 0.29 y 0.71
}

# Exact rationals are PRINTED AS RATIONALS by the runs, so they are checked as the literal
# string "a/b" rather than through a decimal expansion no run ever emits.
FRACTIONS = ["29/8", "21/8", "15/8", "11/8", "9/8", "7/8", "5/8", "3/8", "1/8",
             "25/12", "15/16", "3/4", "1/4", "1/2", "9/2", "5/4", "1/3", "2/3",
             "4/3", "1/16", "5/16", "3/8"]


def tokens(s):
    # La cola de Part VI era (?![\w.]), y con ella un numero AL FINAL DE UNA FRASE no existe:
    # el archivo de ceiling_ilp dice "the escape divides the ceiling by 2.53." y el punto final
    # hacia fallar el lookahead, asi que el gate declaraba 2.53 sin respaldo estando impreso.
    # Lo que hay que rechazar no es un punto cualquiera: es un punto SEGUIDO DE DIGITO, que
    # significa que estamos partiendo un decimal mas largo (la leccion del 1.816 / 816).
    return set(re.findall(r"(?<![\w.])\d+(?:\.\d+)?(?!\w)(?!\.\d)", s))


HAVE_NUM = []          # every archived token as a float, filled by main()


def rounds_from(t):
    """A paper number printed to FEWER decimals than the run prints it is still backed --
    but only if some archived value actually rounds to it at that precision.  This is a
    check, not a whitelist: it cannot admit a number nothing computed.  `-2.45` in the
    prose is `F'' = -2.4488` in outputs/mh_closed_form.txt, and a list of exceptions would
    have hidden that instead of proving it."""
    if "." not in t:
        return None
    d = len(t.split(".")[1])
    # DOS decimales como minimo.  Con uno solo la regla blanquea: `4.4` es la ecuacion (4.4)
    # de AHMN y quedaba "respaldada" por un 4.41 de otro sitio.  Un numero de una cifra
    # decimal casa con demasiadas cosas; si es legitimo va a ALLOW, con su motivo escrito.
    if d < 2:
        return None
    p = float(t)
    for v in HAVE_NUM:
        if round(v, d) == p and v != p:
            return "%.*f" % (min(d + 4, 10), v)
    return None


def audit(tex, have, HAVE_RAW):
    s = open(tex, encoding="utf-8").read()
    # (?<!\\): an ESCAPED percent is content, not a comment
    s = re.sub(r"(?<!\\)%.*", "", s)
    s = re.sub(r"\\begin\{thebibliography\}.*", "", s, flags=re.S)
    # layout is not data
    s = re.sub(r"\\definecolor\{[^}]*\}\{RGB\}\{[^}]*\}", "", s)
    s = re.sub(r"[\d.]+\\textwidth", "", s)
    s = re.sub(r"\\(geometry|includegraphics|usepackage|documentclass)\[[^\]]*\]", "", s)
    s = re.sub(r"p\{[\d.]+\\textwidth\}|\[[\d.]+pt\]|\{[\d.]+cm\}", "", s)
    s = re.sub(r"\\setlength\{\\[a-zA-Z]+\}\{[^}]*\}", "", s)
    s = re.sub(r"\\rowcolor\{[^}]*\}|\\cellcolor\{[^}]*\}|\\rowcolors\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", s)
    s = re.sub(r"[a-zA-Z]+!\d+", "", s)                    # colour mixes: hdrblue!12
    s = re.sub(r"\\\\\[\d+pt\]", "", s)                    # \\[7pt]
    s = re.sub(r"\\fboxrule|\\fboxsep", "", s)
    # the TikZ chain of Figure 1 is a drawing: distances, seps, widths and arrow lengths
    s = re.sub(r"(node distance|inner sep|minimum height|text width|length)=[\d.]+(mm|pt)", "", s)
    # a thin space inside a numeral is typesetting: 18\,648 IS the number 18648
    s = re.sub(r"(?<=\d)\\,(?=\d)", "", s)

    found, missing, allowed, rounded, frac = [], [], [], [], []
    for t in sorted(tokens(s), key=lambda x: (-len(x), x)):
        if len(t) < 2 and t not in ALLOW:
            continue
        if t in ALLOW:
            allowed.append(t)
        elif t in have:
            found.append(t)
        elif rounds_from(t):
            # backed, but printed to fewer decimals than the run prints it
            rounded.append("%s <- %s" % (t, rounds_from(t)))
        else:
            missing.append(t)
    for f in FRACTIONS:
        if f in s:
            (frac if f in HAVE_RAW else missing).append("the rational %s" % f)

    print("%s: every printed number against %s" % (tex, ", ".join(d + "/*" for d in OUT)))
    print("  greppable in an archived run : %d" % len(found))
    print("  backed, rounded when printed : %d" % len(rounded))
    for r in rounded:
        print("     %s" % r)
    print("  exact rationals, as printed  : %d" % len(frac))
    print("  declared non-measurements    : %d" % len(allowed))
    print("  NOT FOUND                    : %d" % len(missing))
    for t in missing:
        print("     %-14s <-- archive its run, or remove it from the paper" % t)
    return missing


def main():
    corpus = ""
    for d in OUT:
        for fn in sorted(os.listdir(d)):
            p = os.path.join(d, fn)
            if os.path.isfile(p):
                corpus += open(p, encoding="utf-8", errors="ignore").read() + "\n"
    have = tokens(corpus)
    # a run that printed 2.73e+00 DID archive 2.73
    have |= set(re.findall(r"(?<![\w.])(\d+(?:\.\d+)?)e[+-]?\d+", corpus))
    HAVE_NUM.extend(sorted({float(x) for x in have}))
    bad = []
    for tex in TEXS:
        bad += audit(tex, have, corpus)
        print()
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
