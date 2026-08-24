#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EL MARCADOR: mira que ha cambiado y corre SOLO las compuertas que ese cambio puede romper.

Por que existe.  Hasta el 24-ago-2026 la rutina era "corro las catorce" despues de cualquier
cambio.  Ese dia se edito CUATRO FRASES de prosa y a continuacion se lanzo un guion de treinta
minutos que verifica que 72 programas siguen reproduciendo sus salidas byte a byte.  No podia
encontrar nada: no se habia tocado ni un programa.  Se lanzo porque tocaba.  Un programa que se
corre porque toca ha dejado de ser una herramienta y es un horario.

La correspondencia de abajo NO esta supuesta: sale de leer que ficheros abre cada compuerta
(`grep` sobre los check_*.py, 24-ago-2026).  Si una compuerta cambia lo que lee, esta tabla
miente, y por eso hay un control que lo detecta ---ver comprobar_el_mapa().

Uso:
    python comprobar.py                 # segun lo que hay sin comprometer (git status)
    python comprobar.py --desde HEAD~1  # segun lo que cambio desde ese punto
    python comprobar.py --todo          # las catorce, cuando de verdad quieras las catorce
    python comprobar.py --solo-decir    # dice que correria y no corre nada
"""
from __future__ import print_function
import argparse, io, os, re, subprocess, sys, time

sys.stdout.reconfigure(encoding="utf-8")

AQUI = os.path.dirname(os.path.abspath(__file__))
PAPEL = os.path.join(AQUI, "paper")
RAIZ_GIT = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))

# ---------------------------------------------------------------- que lee cada compuerta
# derivado de los propios ficheros, no supuesto.  cada entrada: que clases de cambio la afectan.
LEE = {
    "check_attribution": {"tex", "papers"},
    "check_branch":      {"tex", "figura", "log"},
    "check_figures":     {"figura", "pdf", "salida", "log"},
    "check_formulas":    {"tex", "salida", "log"},
    "check_inline":      {"tex"},
    "check_keyeq":       {"tex"},
    "check_layout":      {"pdf", "log"},
    "check_narrative":   {"tex", "log"},
    "check_numbers":     {"tex", "salida"},
    "check_parity":      {"tex", "figura"},
    "check_prose":       {"pdf"},
    "check_refs":        {"tex"},
    "check_reproduces":  {"tex_guiones", "guion", "salida", "log"},
    "check_scripts":     {"tex_guiones", "guion", "salida"},
    "check_sync":        {"tex"},
}

# la unica cara: cuanto cuesta cada una.  el marcador no puede ahorrarte tiempo si no sabe donde esta.
CARA = {"check_reproduces": u"~30 min, corre 72 programas"}


def clasificar(ruta):
    """De una ruta de git a una CLASE de cambio.  None = no afecta a ninguna compuerta."""
    r = ruta.replace("\\", "/")
    if "/part_vii/" not in r and not r.startswith("research/smeft_formalization/"):
        return None
    base = r.rsplit("/", 1)[-1]
    if "/_papers/" in r:
        return "papers"
    if re.match(r"su7_hierarchy(_es)?\.tex$", base):
        return "tex"
    if re.match(r"su7_hierarchy(_es)?\.pdf$", base):
        return "pdf"
    if re.match(r"su7_hierarchy(_es)?\.(log|aux)$", base):
        return "log"
    if base.startswith("fig_") and base.endswith(".pdf"):
        return "figura"
    if "/outputs/" in r:
        return None if base.startswith("_g_") else "salida"   # los _g_ son archivo de compuertas
    if base.startswith("check_") and base.endswith(".py"):
        return "compuerta"
    if base.endswith(".py"):
        # un .py solo es un GUION DE CALCULO si el articulo lo cita.  comprobar.py,
        # lectura_citas.py y demas herramientas no pueden romper la reproduccion de nada,
        # y si se cuelan aqui disparan el reproductor de media hora para nada.
        return "guion" if base in _citados_por_el_articulo() else None
    return None


_CITADOS = None


def _citados_por_el_articulo():
    """Los nombres de .py que aparecen en el texto de las dos ediciones.

    Es la misma fuente que usa check_reproduces para saber que tiene que correr, asi que las
    dos listas no pueden divergir."""
    global _CITADOS
    if _CITADOS is None:
        _CITADOS = set()
        for f in ("su7_hierarchy.tex", "su7_hierarchy_es.tex"):
            try:
                t = io.open(os.path.join(PAPEL, f), encoding="utf-8").read()
            except IOError:
                continue
            _CITADOS |= set(re.findall(r"[A-Za-z0-9_]+\.py", t))
    return _CITADOS


def lista_de_guiones_cambio(desde):
    """Solo el .tex? Entonces el reproductor solo importa si cambio QUE guiones cita el articulo.

    check_reproduces lee el .tex para una sola cosa: saber que programas tiene que correr.  Si
    editas prosa y la lista de programas citados es identica, ese guion de media hora no puede
    encontrar nada.  El 24-ago-2026 se lanzo igualmente, tras editar cuatro frases.
    """
    ref = desde or "HEAD"
    for f in ("su7_hierarchy.tex", "su7_hierarchy_es.tex"):
        ruta = "research/smeft_formalization/part_vii/paper/" + f
        try:
            antes = subprocess.check_output(["git", "show", "%s:%s" % (ref, ruta)],
                                            cwd=RAIZ_GIT, stderr=subprocess.DEVNULL)
            antes = antes.decode("utf-8", "replace")
        except subprocess.CalledProcessError:
            return True                      # no puedo comparar: no me la juego
        ahora = io.open(os.path.join(PAPEL, f), encoding="utf-8").read()
        rx = r"[A-Za-z0-9_]+\.py"
        if set(re.findall(rx, antes)) != set(re.findall(rx, ahora)):
            return True
    return False


def cambiado(desde):
    if desde:
        cmd = ["git", "diff", "--name-only", desde, "--", "research/smeft_formalization/part_vii"]
    else:
        cmd = ["git", "status", "--porcelain", "--untracked-files=all",
               "--", "research/smeft_formalization/part_vii"]
    out = subprocess.check_output(cmd, cwd=RAIZ_GIT).decode("utf-8", "replace")
    rutas = []
    for l in out.splitlines():
        l = l.rstrip()
        if not l:
            continue
        rutas.append(l[3:].strip().strip('"') if not desde else l.strip())
    return rutas


def comprobar_el_mapa():
    """Control que puede fallar: que la tabla LEE nombre exactamente las compuertas que existen.

    Si manyana aparece un check_nuevo.py y nadie toca esta tabla, el marcador lo saltaria en
    silencio ---y un saltado en silencio es peor que no tener marcador---.  Asi que revienta.
    """
    reales = {f[:-3] for f in os.listdir(PAPEL)
              if f.startswith("check_") and f.endswith(".py")}
    faltan, sobran = reales - set(LEE), set(LEE) - reales
    if faltan or sobran:
        print(u"\n  EL MAPA MIENTE.")
        if faltan:
            print(u"  compuertas que existen y no estan en la tabla: %s" % u", ".join(sorted(faltan)))
        if sobran:
            print(u"  entradas de la tabla sin compuerta detras:      %s" % u", ".join(sorted(sobran)))
        print(u"  arregla comprobar.py antes de fiarte de el.\n")
        sys.exit(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--desde", default=None)
    ap.add_argument("--todo", action="store_true")
    ap.add_argument("--solo-decir", action="store_true")
    a = ap.parse_args()

    comprobar_el_mapa()

    if a.todo:
        clases, rutas = set(LEE and {"tex", "tex_guiones", "pdf", "log", "figura", "salida", "guion", "papers"}), []
        print(u"\n  LAS CATORCE, a mano.\n")
    else:
        rutas = cambiado(a.desde)
        clases = set()
        for r in rutas:
            c = clasificar(r)
            if c:
                clases.add(c)
        print(u"\n  %d fichero(s) cambiado(s)%s" % (len(rutas), u" desde %s" % a.desde if a.desde else u" sin comprometer"))
        for r in sorted(rutas)[:14]:
            c = clasificar(r)
            print(u"     %-9s %s" % (u"[%s]" % c if c else u"[ - ]", r.replace("\\", "/").split("part_vii/")[-1]))
        if len(rutas) > 14:
            print(u"     ... y %d mas" % (len(rutas) - 14))

    _sin_guiones_nuevos = False
    if "tex" in clases and not a.todo:
        if lista_de_guiones_cambio(a.desde):
            clases.add("tex_guiones")
        else:
            _sin_guiones_nuevos = True

    if "compuerta" in clases:
        print(u"\n  Has tocado una COMPUERTA. Eso no lo decide el marcador: corre --todo.")
        clases |= {"tex", "tex_guiones", "pdf", "log", "figura", "salida", "guion", "papers"}

    toca = sorted(g for g, necesita in LEE.items() if necesita & clases)
    if not toca:
        print(u"\n  Ningun cambio afecta a ninguna compuerta. No hay nada que correr.\n")
        return 0

    print(u"\n  clases de cambio: %s" % u", ".join(sorted(clases)))
    print(u"  correrian %d de %d:" % (len(toca), len(LEE)))
    for g in toca:
        print(u"     %-20s %s" % (g, CARA.get(g, u"")))
    if _sin_guiones_nuevos and "check_reproduces" not in toca:
        print(u"  (el .tex cambio, pero cita los mismos guiones: el reproductor no puede")
        print(u"   encontrar nada, y por eso no esta en la lista)")
    saltadas = sorted(set(LEE) - set(toca))
    if saltadas:
        print(u"  se saltan %d: %s" % (len(saltadas), u", ".join(s.replace("check_", "") for s in saltadas)))

    if a.solo_decir:
        print()
        return 0

    print()
    fallos = []
    for g in toca:
        t0 = time.time()
        p = subprocess.run([sys.executable, g + ".py"], cwd=PAPEL,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        ok = p.returncode == 0
        print(u"  %-20s %s   %5.1fs" % (g, u"verde" if ok else u"FALLA", dt))
        if not ok:
            fallos.append((g, p.stdout.decode("utf-8", "replace")))

    if fallos:
        for g, salida in fallos:
            print(u"\n%s\n  %s\n%s" % (u"=" * 78, g, u"=" * 78))
            print(u"\n".join(salida.splitlines()[-25:]))
        print(u"\n  %d compuerta(s) en rojo.\n" % len(fallos))
        return 1
    print(u"\n  todas las que tocaban, verdes.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
