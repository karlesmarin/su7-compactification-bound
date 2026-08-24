#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lanza check_reproduces a terminar, sobreviva o no la sesion, sin dejar un parcial enganyoso.

Escribe en <base>.EN_CURSO.txt/.err mientras corre y solo lo renombra a <base>.txt/.err cuando
el proceso termina, dejando en la ultima linea el codigo de salida.  Asi, si esto se interrumpe,
lo que queda en outputs/ lleva EN_CURSO en el nombre y nadie lo confunde con una corrida entera.

Uso:  python _correr_reproduces.py 2026-08-25
"""
from __future__ import print_function
import os, subprocess, sys, io

AQUI = os.path.dirname(os.path.abspath(__file__))
PAPEL = os.path.join(AQUI, "paper")
SALIDAS = os.path.join(AQUI, "outputs")

fecha = sys.argv[1] if len(sys.argv) > 1 else "sin_fecha"
base = os.path.join(SALIDAS, "_g_check_reproduces_%s" % fecha)

tmp_out, tmp_err = base + ".EN_CURSO.txt", base + ".EN_CURSO.err"
fin_out, fin_err = base + ".txt", base + ".err"

for f in (tmp_out, tmp_err, fin_out, fin_err):
    if os.path.exists(f):
        os.remove(f)

env = dict(os.environ)
env["PYTHONIOENCODING"] = "utf-8"
env["PYTHONUNBUFFERED"] = "1"          # que el fichero no vaya media hora por detras

with io.open(tmp_out, "wb") as o, io.open(tmp_err, "wb") as e:
    codigo = subprocess.call([sys.executable, "check_reproduces.py"],
                             cwd=PAPEL, stdout=o, stderr=e, env=env)

with io.open(tmp_out, "a", encoding="utf-8", newline="\n") as o:
    o.write(u"\n--- codigo de salida: %d ---\n" % codigo)

os.rename(tmp_out, fin_out)
os.rename(tmp_err, fin_err)
print("terminado, codigo %d" % codigo)
sys.exit(codigo)
