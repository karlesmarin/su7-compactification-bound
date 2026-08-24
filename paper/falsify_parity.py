#!/usr/bin/env python3
"""Descontar las citas textuales, no puede haber abierto un agujero.

Se planta una palabra funcional inglesa FUERA de comillas y se exige que check_parity la cace.
Y una DENTRO de comillas, que debe seguir pasando.
"""
import pathlib
import shutil
import subprocess
import sys

P = pathlib.Path(r"E:\proyectos\Curiosity\research\smeft_formalization\part_vii\paper")
TEX = P / "su7_hierarchy_es.tex"
orig = TEX.read_text(encoding="utf-8")


def run():
    r = subprocess.run([sys.executable, "check_parity.py"], cwd=P, capture_output=True, text=True)
    return r.returncode


ANCLA = "La contribución, dicha sin rodeos.}"
assert ANCLA in orig

CASES = [
    ("una palabra inglesa FUERA de comillas", ANCLA + " This is the contribution.", True),
    ("la misma DENTRO de comillas", ANCLA + " ``This is the contribution.''", False),
]

try:
    code = run()
    print("baseline: exit=%d  %s" % (code, "PASS" if code == 0 else "*** ya falla ***"))
    if code:
        sys.exit(1)
    bad = 0
    for name, new, want_fail in CASES:
        TEX.write_text(orig.replace(ANCLA, new, 1), encoding="utf-8")
        c = run()
        ok = (c != 0) == want_fail
        bad += not ok
        print("  %-40s exit=%d  %s" % (name, c, "correcto" if ok else "*** MAL ***"))
finally:
    TEX.write_text(orig, encoding="utf-8")

print("\nrestaurado: exit=%d" % run())
sys.exit(1 if bad else 0)
