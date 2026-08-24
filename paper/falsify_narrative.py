#!/usr/bin/env python3
"""falsify_narrative.py -- can check_narrative.py fail?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

A gate that has only ever been seen to pass is not evidence.  This plants each of the three
defects check_narrative.py claims to catch, one at a time, into the real English source, and
demands a non-zero exit for every one of them -- then puts the file back and re-runs the gate to
prove the restoration is clean.  [[a-falsification-suite-beats-a-passing-control]]

Run:  python falsify_narrative.py     (from part_vii/paper/)
"""
import pathlib
import shutil
import subprocess
import sys

P = pathlib.Path(__file__).resolve().parent
TEX = P / "su7_hierarchy.tex"
BAK = P / "_falsify_backup.tex"


def run():
    r = subprocess.run([sys.executable, "check_narrative.py"], cwd=P,
                       capture_output=True, text=True)
    return r.returncode, r.stdout


CASES = [
    ("N1 a number only the abstract knows",
     "so the Wilson-line vacuum \\emph{is} the electroweak hierarchy.",
     "so the Wilson-line vacuum \\emph{is} the electroweak hierarchy, at $77.7731$~TeV."),
    ("N2 a number only the intro knows",
     "so $\\amin$ \\emph{is} the electroweak hierarchy. A small phase",
     "so $\\amin$ \\emph{is} the electroweak hierarchy, all $88.8842$ of it. A small phase"),
    ("N3 a dangling section promise",
     "\\item A second, independent anchor (\\S\\ref{sec:anchor})",
     "\\item A second, independent anchor (\\S\\ref{sec:nosuchthing})"),
    # N4 -- la que costo una tarde: reutilizar una etiqueta.  LaTeX compila, avisa en el .log y
    # \eqref apunta a la ultima definicion.  El guion NO recompila, asi que esta planta la
    # duplicacion en la FUENTE y comprueba que N4 la ve sin ayuda del registro.
    ("N4 a label defined twice",
     "\\begin{equation}\\label{eq:tulimit}",
     "\\begin{equation}\\label{eq:chisq}"),
    # N5 -- una edicion que se queda atras.  Es EXACTAMENTE lo que paso cuatro veces el 23 de
    # agosto: se corrige el ingles, se olvida el espanol, y las nueve compuertas siguen verdes
    # porque ninguna compara una edicion con la otra.
    ("N5 the editions disagree on a number",
     "  \\frac{1}{R_5}\\;>\\;10.86\\ \\text{TeV}\\,,",
     "  \\frac{1}{R_5}\\;>\\;11.20\\ \\text{TeV}\\,,"),
]

shutil.copy(TEX, BAK)
orig = TEX.read_text(encoding="utf-8")
try:
    code, out = run()
    print("baseline: exit=%d  %s" % (code, "PASS" if code == 0 else "*** already failing ***"))
    if code != 0:
        sys.exit("cannot falsify from a failing baseline")
    bad = 0
    for name, old, new in CASES:
        if old not in orig:
            print("  %-38s *** anchor text not found; test is broken ***" % name)
            bad += 1
            continue
        TEX.write_text(orig.replace(old, new, 1), encoding="utf-8")
        code, out = run()
        ok = code != 0
        bad += not ok
        print("  %-38s exit=%d  %s" % (name, code, "CAUGHT" if ok else "*** MISSED ***"))
finally:
    TEX.write_text(orig, encoding="utf-8")
    BAK.unlink()

code, _ = run()
print("\nrestored: exit=%d" % code)
sys.exit(1 if bad or code else 0)
