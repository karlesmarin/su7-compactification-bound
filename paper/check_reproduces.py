#!/usr/bin/env python3
"""check_reproduces.py - novena compuerta: los guiones citados VUELVEN A DAR lo archivado.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

Las otras ocho leen el .tex, el .log y los PDF.  Ninguna ejecuta nada.  `check_numbers.py`
comprueba que todo numero impreso es greppable en `../outputs/`, y `check_scripts.py` que
todo guion citado existe y tiene salida archivada -- pero las dos creen al archivo.  Si el
guion deja de producirlo, las dos siguen verdes y la frase de Data availability ---*every
displayed number regenerates from the ancillary scripts*--- pasa a ser falsa sin que nada
suene.

Eso es exactamente lo que paso, y es el motivo de esta compuerta.  `su7_vacuum.py` se edito
despues de escribirse cuatro guiones que rehusan sus funciones haciendo `exec` de sus
primeras 102 lineas; a partir de esa edicion la linea 102 caia DENTRO de `minimise()`, que
compilaba sin su `return` y devolvia `None` a todo el mundo.  `su7_sixth_row.py` -- la fila
sexta pre-registrada de la seccion 7 -- y `su7_wedge_direction.py` -- la figura 5 --
reventaban.  Los numeros del articulo estaban bien; lo roto era poder rehacerlos.

Cada guion se ejecuta a un fichero NUEVO y se compara con el archivo; el archivo no se
toca nunca -- es contra lo que `check_numbers.py` greppa, y reescribirlo con la corrida de
hoy convertiria la comprobacion en una tautologia.

  python check_reproduces.py             ejecuta todo (lento: ~10 min)
  python check_reproduces.py --falsify   demuestra que la compuerta PUEDE fallar
  python check_reproduces.py NOMBRE ...  solo esos guiones

Lo que NO cubre, dicho en voz alta: los `.sage` (necesitan el contenedor) y los guiones que
el articulo cita sin salida archivada.  Los dos casos se imprimen, uno a uno, con su motivo.
"""
import difflib
import os
import pathlib
import re
import subprocess
import sys
import tempfile

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent          # part_vii/paper
PART = HERE.parent                                       # part_vii
OUT = PART / "outputs"
TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]

# LOS GUIONES NO ESTAN TODOS EN part_vii/.  Los nueve del recast de colisionador viven en
# part_vii/tu_limit/, y esta compuerta solo miraba en part_vii, de modo que los daba por
# inexistentes: `no comprobado  scan_mkk.py  el guion no existe`.  Un guion que la compuerta no
# encuentra se lee igual que uno que no hace falta comprobar, y esos nueve son justo los que mas
# se han movido.  Cada uno se ejecuta desde SU PROPIO directorio, que es como estan escritos.
# [[a-gate-outside-the-list-does-not-exist]]
SEARCH = [PART, PART / "tu_limit"]

# Y tres de ellos leen un ingrediente que este repositorio NO puede llevar: la salida de CIJET,
# cuya licencia permite el uso academico y prohibe la redistribucion.  Sin el fichero, el guion
# no falla por estar roto -- falla por no tener con que correr, y esta compuerta tiene que
# distinguir las dos cosas o convierte una restriccion legal en un defecto de codigo.
CIJET_GRID = pathlib.Path("E:/proyectos/Curiosity/tools/cijet/CIJET/data/grid_fitALresults.dat")
NEEDS = {
    "cijet_control.py": (CIJET_GRID, "necesita la salida de CIJET, que su licencia prohibe redistribuir"),
    "alphas_diagnosis.py": (CIJET_GRID, "necesita la salida de CIJET, que su licencia prohibe redistribuir"),
    "cms_ci_closure.py": (CIJET_GRID, "necesita la salida de CIJET, que su licencia prohibe redistribuir"),
}


# --------------------------------------------------------------------------- variantes
# Guiones que el articulo usa CON BANDERA.  Sin esto la variante que se publica no se
# comprueba nunca: se comprueba la otra, y pasa.
#
#   make_fig_collapse.py            -> familia sintetica, cargas 1-3   -> make_fig_collapse.txt
#   make_fig_collapse.py --lattice  -> los ocho multipletes de KM      -> ..._lattice.txt
#                                      ^^ ESTA es fig_collapse_lattice.pdf, la del articulo
#
# La clave del diccionario es el nombre del guion; el valor, la lista de (argumentos, archivo).
VARIANTES = {
    "make_fig_collapse.py": [([], "make_fig_collapse.txt"),
                             (["--lattice"], "make_fig_collapse_lattice.txt")],
}


def locate(name):
    """(fichero, directorio desde el que ejecutarlo) o (None, None) si no esta en ningun sitio."""
    for d in SEARCH:
        p = d / name
        if p.exists():
            return p, d
    return None, None


def cited():
    """los .py y .sage que el articulo nombra en \\texttt{} -- la lectura de check_scripts.py."""
    out = set()
    for t in TEXS:
        s = (HERE / t).read_text(encoding="utf-8")
        for m in re.findall(r"\\texttt\{([^}]*?\.(?:py|sage))\}", s):
            n = "".join(m.replace("\\allowbreak", "").replace("\\_", "_").split())
            if "/" not in n and not n.startswith("."):
                out.add(n)
    return sorted(out)


def run(name, args=()):
    """devuelve (returncode, stdout) ejecutando el guion desde SU directorio, como esta escrito."""
    src, cwd = locate(name)
    p = subprocess.run([sys.executable, str(src)] + list(args), cwd=str(cwd),
                       capture_output=True, timeout=3600)
    return p.returncode, p.stdout.decode("utf-8", "replace").replace("\r\n", "\n")


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    falsify = "--falsify" in sys.argv
    names = argv or cited()

    P("guiones citados por el articulo: %d" % len(cited()))
    P("")

    same, diff, crash, skipped = [], [], [], []
    # cada guion se corre una vez por variante: sin bandera si no tiene, y una por cada
    # bandera declarada arriba.  `etiqueta` es lo que se imprime, para que se vea cual paso.
    trabajos = []
    for n in names:
        for args, archivo in VARIANTES.get(n, [([], os.path.splitext(n)[0] + ".txt")]):
            trabajos.append((n, tuple(args), archivo,
                             n if not args else "%s %s" % (n, " ".join(args))))

    for n, args, archivo, etiqueta in trabajos:
        src, _cwd = locate(n)
        arch = OUT / archivo
        if src is None:
            skipped.append((etiqueta, "el guion no existe en %s" % " ni ".join(d.name for d in SEARCH)))
            continue
        if n.endswith(".sage"):
            skipped.append((etiqueta, "Sage: necesita el contenedor, fuera del alcance de esta compuerta"))
            continue
        if n in NEEDS and not NEEDS[n][0].exists():
            skipped.append((etiqueta, NEEDS[n][1]))
            continue
        if not arch.exists():
            skipped.append((etiqueta, "no tiene salida archivada en outputs/"))
            continue
        try:
            rc, got = run(n, args)
        except subprocess.TimeoutExpired:
            crash.append((etiqueta, "TIMEOUT a los 3600 s"))
            P("  %-34s TIMEOUT" % etiqueta)
            continue
        want = arch.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
        # UN CODIGO DE SALIDA NO CERO NO ES UN FALLO DE REPRODUCIBILIDAD.  Varios guiones de
        # esta serie salen con 1 cuando un control suyo cae, y eso es lo que deben hacer: el
        # veredicto honesto de holonomy_skeleton.py es que C1 y C2 fallan.  Lo que esta
        # compuerta mide es si el guion vuelve a producir SU ARCHIVO, no si le gusta el
        # resultado.  Antes marcaba REVIENTA sin llegar a comparar la salida, de modo que un
        # negativo bien medido se leia como codigo roto.
        if rc and got != want:
            crash.append((etiqueta, "salio con codigo %d y ademas difiere" % rc))
            P("  %-34s REVIENTA (codigo %d)" % (etiqueta, rc))
            continue
        if got == want:
            if rc:
                P("  %-34s reproduce su archivo (control rojo, codigo %d: es su veredicto)"
                  % (etiqueta, rc))
                same.append(etiqueta)
                continue
            same.append(etiqueta)
            P("  %-34s reproduce su archivo" % etiqueta)
        else:
            diff.append(etiqueta)
            P("  %-34s *** DIFIERE DEL ARCHIVO ***" % etiqueta)
            for line in list(difflib.unified_diff(want.splitlines(), got.splitlines(),
                                                  "archivo", "corrida de hoy",
                                                  lineterm="", n=1))[:20]:
                P("        %s" % line)

    P("")
    for n, why in skipped:
        P("  no comprobado  %-30s %s" % (n, why))
    if skipped:
        P("")
    P("  reproducen %d | difieren %d | revientan %d | no comprobados %d"
      % (len(same), len(diff), len(crash), len(skipped)))
    for n, why in crash:
        P("     REVIENTA  %-28s %s" % (n, why))

    if falsify:
        P("")
        P("  --falsify: la compuerta tiene que poder fallar, asi que se le miente.")
        victim = same[0] if same else None
        if victim is None:
            P("     nada verde que estropear; falsificacion no concluyente")
            return 1
        arch = OUT / (os.path.splitext(victim)[0] + ".txt")
        want = arch.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
        rc, got = run(victim)
        mutated = want.replace("0", "9", 1)
        P("     %s contra un archivo con UN digito cambiado: %s"
          % (victim, "DETECTADO" if got != mutated else "*** NO DETECTADO ***"))
        if got == mutated:
            return 1

    if diff or crash:
        P("")
        P("REPRODUCIBILIDAD ROTA: %d guion(es) ya no dan lo que el articulo cita."
          % (len(diff) + len(crash)))
        return 1
    P("")
    P("todo guion citado que puede ejecutarse aqui vuelve a dar su salida archivada, byte a byte")
    return 0


if __name__ == "__main__":
    sys.exit(main())
