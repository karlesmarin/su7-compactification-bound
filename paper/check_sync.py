#!/usr/bin/env python3
"""check_sync.py - una correccion que se hace en una edicion, .llego a la otra?

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

POR QUE EXISTE.  Es el modo de fallo que mas veces ha aparecido en este articulo, cuatro dias
seguidos, y el unico que ninguna compuerta cazaba:

    22-ago  `613a507` anyade DOS parrafos a la S4 y la S5 inglesas.  Su propio mensaje dice
            "las dos ligaduras estan escritas en las dos ediciones, en los dos sentidos".
            No lo estaban: +29 lineas en la inglesa contra +12 en la castellana.
    23-ago  `0bee2ae` baja la condicionalidad al resto del texto.  Un parrafo entero de la S13
            ---"cerrar esto pide un calculo, y mientras tanto todo lo que necesita que el
            desplazamiento gauge sea impar queda condicionado a la ec. (68)"--- no cruzo.
    23-ago  `d5b5f58` retira "una maquina que cubra 10-20 TeV zanja la clase", una de las cuatro
            afirmaciones falsas *a nuestro favor* que ese commit destapo.  La castellana la
            conservo un dia entero, dentro de una caja `keyeq`.
    24-ago  y al reves: `b065ce5` arregla la castellana y deja la inglesa con la misma frase.

`check_parity` compara ELEMENTOS ---secciones, figuras, tabulares---, `check_narrative` NUMEROS,
`check_formulas` FORMULAS, `check_inline` la matematica en linea, `check_keyeq` las cajas.  Ninguna
compara ARGUMENTOS, y una retractacion tiene los mismos elementos, los mismos numeros y las mismas
formulas que la afirmacion que retracta.  [[translating-finds-what-review-missed]]

COMO FUNCIONA.  Dos medidas, y ninguna entiende el texto.

  S1  LA FIRMA POR PARRAFO EMPAREJADO.  Empareja los parrafos de prosa de las dos ediciones por
      posicion dentro de cada seccion y compara lo que una traduccion NO mueve: cuantos `\\emph`,
      cuantos `\\textbf`, cuantos `\\eqref`, cuantas `\\cite`, cuantos `\\texttt`, cuantos numeros
      con decimales.  Una traduccion mueve las palabras; no mueve un `\\eqref` ni un numero.  Si la
      castellana lleva dos citas donde la inglesa lleva cuatro, ahi falta texto.
      Encontro el parrafo de la S13 y el `gate_km_ghosts.py` perdido de la S11.

  S2  LA ASIMETRIA POR COMMIT.  Para cada commit que toca alguna de las dos ediciones, cuantas
      lineas anyadio a cada una.  No es prueba ---el castellano corre un 4 % mas largo y las lineas
      se parten distinto--- es una LISTA DE SITIOS DONDE MIRAR, ordenada por cuanto se desviaron.
      Localizo `613a507` en un minuto.  No hace fallar: informa.

LO QUE NO PUEDE HACER.  Si las dos ediciones dicen lo mismo y lo mismo esta mal, pasa.  Y si una
correccion cambia palabras sin tocar ninguna cita, ningun numero y ninguna enfasis ---que es
exactamente lo que fue el "no se sostiene" contra "la hipotesis no se cumple" de la S11--- S1 no la
ve.  Contra eso solo hay leer las dos en paralelo.  Esta compuerta es un suelo, no un sustituto.

Uso:  python check_sync.py     (desde part_vii/paper/)
"""
import os
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
EN, ES = "su7_hierarchy.tex", "su7_hierarchy_es.tex"
P = lambda *a: print(*a, flush=True)

FLOAT = re.compile(r"\\begin\{(figure|table|tikzpicture|longtable|tabular)\*?\}.*?"
                   r"\\end\{\1\*?\}", re.S)

# Diferencias comprobadas A MANO que NO son defectos.  Cada una con su razon escrita: una lista
# blanca sin razones deja de ser una lista blanca.  [[a-control-that-cannot-fail]]
ACEPTADOS = {
    # La inglesa escribe "the potential of \cite{KM25} already reduced in that limit" y la
    # castellana "su potencial ya reducido en ese limite": el posesivo sustituye a la cita, que
    # esta dos frases antes en el mismo parrafo.  Correcto en castellano.
    ("sec:setting", 0, "cite"),
    # La inglesa entrecomilla con ``...'' y la castellana con \emph{«...»}, de modo que cada cita
    # textual traducida anyade dos \emph que la inglesa no lleva.  Es el convenio tipografico de
    # la edicion, no contenido.
    ("sec:ceiling", 13, "emph"),
    # La castellana escribe `\emph{nuisance}` dos veces donde la inglesa escribe "nuisance
    # parameters" a secas: es un prestamo, y el convenio de la edicion es marcarlo en cursiva
    # ---lo hace ya en otros cuatro sitios---.  El ingles no toma prestado nada.
    ("sec:collider", 39, "emph"),
}


def strip_comments(s):
    """Fuera los comentarios.  Las lineas que son SOLO comentario se borran ENTERAS, no se dejan
    en blanco: un bloque de notas dentro de un parrafo lo partiria en dos y el emparejamiento se
    descuadraria a partir de ahi.  Las dos ediciones no llevan las notas en los mismos sitios, de
    modo que el descuadre no seria ni siquiera simetrico."""
    s = "\n".join(ln for ln in s.split("\n") if not ln.lstrip().startswith("%"))
    return re.sub(r"(?<!\\)%.*", "", s)


def sections(path):
    """[(etiqueta, cuerpo)] por seccion, con la etiqueta \\label{sec:...} como nombre."""
    s = strip_comments(open(os.path.join(HERE, path), encoding="utf-8").read())
    s = s.split(r"\begin{thebibliography}")[0]
    out, name, buf = [], "portada", []
    for ln in s.split("\n"):
        m = re.match(r"\\section\*?\{.*?\}(?:\\label\{(.+?)\})?", ln)
        if m:
            out.append((name, "\n".join(buf)))
            name, buf = m.group(1) or "sin-label", []
            continue
        buf.append(ln)
    out.append((name, "\n".join(buf)))
    return out


def paragraphs(body):
    body = FLOAT.sub(" ", body)
    # los comentarios ya no estan, pero dejaron lineas en blanco que partirian un parrafo en dos
    body = re.sub(r"\n{3,}", "\n\n", body)
    return [p.strip() for p in re.split(r"\n\s*\n", body) if len(p.strip()) > 120]


def sig(p):
    return {
        "emph": len(re.findall(r"\\emph\{", p)),
        "bf": len(re.findall(r"\\textbf\{", p)),
        "eqref": len(re.findall(r"\\eqref\{", p)),
        "cite": len(re.findall(r"\\cite\{", p)),
        "tt": len(re.findall(r"\\texttt\{", p)),
        "num": len(set(re.findall(r"(?<![\w.])\d+\.\d+", p))),
    }


# los campos DUROS: una diferencia de uno ya es senyal, porque no dependen del idioma.
# `emph` y `bf` son blandos: la puntuacion de una cita traducida los mueve, asi que ahi hace
# falta una diferencia de dos.
DUROS = ("eqref", "cite", "tt", "num")


def s1():
    P("=" * 92)
    P("S1  LA FIRMA POR PARRAFO EMPAREJADO -- lo que una traduccion no mueve")
    P("=" * 92)
    en, es = sections(EN), sections(ES)
    bad = 0
    if len(en) != len(es):
        P("  *** %d secciones en la inglesa y %d en la castellana ***" % (len(en), len(es)))
        return 1
    total = 0
    for (ne, be), (ns, bs) in zip(en, es):
        pe, ps = paragraphs(be), paragraphs(bs)
        if len(pe) != len(ps):
            P("  %-14s %d parrafos EN contra %d ES  <-- FALTA O SOBRA UN PARRAFO"
              % (ne[:14], len(pe), len(ps)))
            bad += 1
            continue
        total += len(pe)
        for i, (a, b) in enumerate(zip(pe, ps)):
            sa, sb = sig(a), sig(b)
            dif = {k: (sa[k], sb[k]) for k in sa if sa[k] != sb[k]
                   and (k in DUROS or abs(sa[k] - sb[k]) >= 2)
                   and (ne, i, k) not in ACEPTADOS}
            if dif:
                bad += 1
                P("  %-14s parrafo %-3d %s" % (ne[:14], i,
                  "  ".join("%s %d/%d" % (k, v[0], v[1]) for k, v in dif.items())))
                P("        EN: %s" % " ".join(a.split())[:82])
                P("        ES: %s" % " ".join(b.split())[:82])
    P("  %d parrafos emparejados en %d secciones, %d con la firma distinta"
      % (total, len(en), bad))
    return bad


def s2():
    P("")
    P("=" * 92)
    P("S2  LA ASIMETRIA POR COMMIT -- informativa, no hace fallar")
    P("=" * 92)
    rel = "research/smeft_formalization/part_vii/paper/"
    repo = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))

    def sh(*a):
        r = subprocess.run(["git"] + list(a), cwd=repo, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        return r.stdout

    log = sh("log", "--format=%H|%ad|%s", "--date=short", "--",
             rel + EN, rel + ES).strip().split("\n")
    if not log or not log[0]:
        P("  no hay historia de git aqui; S2 no puede correr")
        return
    P("  %-9s %-11s %5s %5s  %s" % ("commit", "fecha", "+EN", "+ES", "asunto"))
    n = 0
    for line in log[:40]:
        if "|" not in line:
            continue
        h, date, subj = line.split("|", 2)
        add = {}
        for l in sh("show", "--numstat", "--format=", h, "--",
                    rel + EN, rel + ES).strip().split("\n"):
            p = l.split("\t")
            if len(p) == 3:
                add[os.path.basename(p[2])] = int(p[0]) if p[0] != "-" else 0
        a, b = add.get(EN), add.get(ES)
        if a is None or b is None or a < 4:
            continue
        r = b / a
        marca = "  <-- MIRAR" if (r < 0.60 or r > 1.9) else ""
        if marca:
            n += 1
        P("  %-9s %-11s %5d %5d  %-46s%s" % (h[:7], date, a, b, subj[:46], marca))
    P("")
    P("  %d commit(s) con la castellana desviada mas de un 40 %% de la inglesa." % n)
    P("  No es prueba: el castellano corre un 4 % mas largo y las lineas se parten distinto.")
    P("  Es la lista de sitios donde mirar cuando algo no cuadra.")


# --------------------------------------------------------------------------- falsacion
# Los dos casos REALES del 24-ago-2026 que S1 tiene que ver.  Si dejan de saltar, la firma se ha
# quedado sin campos y esta compuerta ha vuelto a ser un adorno.
# [[a-falsification-suite-beats-a-passing-control]]
FALSIFY = [
    ("la S13, el parrafo de la condicionalidad que no cruzo",
     r"\emph{Closing this takes a calculation, not another check}: a background-gauge determinant "
     r"for one charged orbit, with both parities kept, returning the split and not just the total. "
     r"Until that exists, everything here that needs the gauge offset to be \emph{odd} is "
     r"conditional on \cite{KM25}'s eq.~(68) as published.",
     r"\emph{Pero si podemos decir lo que cuesta la otra rama, y cuesta menos de lo que parece.} "
     r"Correr el mismo programa entero con la semilla candidata es un cambio pequenyo."),
    ("la S11, el gate_km_ghosts.py que se perdio",
     r"Measured in \texttt{gauge\_weight\_origin.py} and \texttt{gate\_km\_ghosts.py}. What "
     r"\S\ref{sec:open} records is a different and narrower gap.",
     r"Medido en \texttt{gauge\_weight\_origin.py}. Lo que la \S\ref{sec:open} recoge es un hueco "
     r"distinto y mas estrecho."),
]


def falsify():
    P("")
    P("=" * 92)
    P("FALSACION -- los dos parrafos del 24-ago, .los ve la firma?")
    P("=" * 92)
    fails = 0
    for name, a, b in FALSIFY:
        sa, sb = sig(a), sig(b)
        dif = [k for k in DUROS if sa[k] != sb[k]]
        ok = bool(dif)
        fails += 0 if ok else 1
        P("  %-48s %s" % (name, ("salta, difieren en " + ",".join(dif)) if ok
                          else "NO LO VE -- la firma se ha quedado ciega"))
    if fails:
        P("\n%d caso(s) real(es) que esta compuerta dejaria pasar." % fails)
    else:
        P("\nlos dos parrafos rancios saltan.")
    return fails


def main():
    bad = s1()
    s2()
    bad += falsify()
    P("")
    P("=" * 92)
    if bad:
        P("PARRAFOS QUE NO DICEN LO MISMO: %d" % bad)
        P("Una correccion que se hace en una edicion y no en la otra deja el articulo diciendo")
        P("dos cosas distintas segun por donde se abra.  Y va en los dos sentidos.")
        P("=" * 92)
        return 1
    P("cada parrafo de la castellana lleva las mismas citas, ecuaciones y numeros que su gemelo.")
    P("=" * 92)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
