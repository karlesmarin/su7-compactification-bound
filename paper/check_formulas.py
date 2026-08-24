#!/usr/bin/env python3
"""check_formulas.py - octava compuerta: las FÓRMULAS, no los números.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

`check_numbers.py` comprueba que cada número impreso esté en una salida archivada, y
`check_parity.py` cuenta entornos. Ninguno mira el interior de una fórmula. Éste sí, y hace
tres cosas que ninguno de los dos hace:

  1. EMPAREJA las fórmulas de las dos ediciones en orden y las compara carácter a carácter
     tras normalizar el espaciado. Una fórmula que difiere entre ediciones es un error en una
     de las dos: o una errata, o una corrección aplicada sólo a un lado. Ésta es la razón de
     ser del fichero.
  2. Comprueba el EQUILIBRIO de delimitadores dentro de cada fórmula ---llaves, \\left/\\right,
     paréntesis--- y que ningún `$` quede impar en una línea de texto.
  3. Comprueba que toda macro `\\newcommand` del preámbulo usada en fórmulas exista en LA OTRA
     edición también, que es como se cuela un `\\rep` sin definir tras copiar un bloque.

Salida en ../outputs/check_formulas.txt. Código de salida 1 si algo falla.
"""
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "outputs" / "check_formulas.txt"
EN, ES = HERE / "su7_hierarchy.tex", HERE / "su7_hierarchy_es.tex"

# entornos de display que llevan matemáticas
ENVS = ("equation", "equation*", "align", "align*", "gather", "gather*", "multline")


def strip_comments(text):
    """quita comentarios TeX sin comerse un \\% escapado (ver la nota de memoria)."""
    out = []
    for line in text.splitlines():
        i, esc = None, False
        for k, ch in enumerate(line):
            if ch == "\\":
                esc = not esc
                continue
            if ch == "%" and not esc:
                i = k
                break
            esc = False
        out.append(line if i is None else line[:i])
    return "\n".join(out)


def norm(s):
    """normaliza espaciado y etiquetas para comparar dos ediciones."""
    s = re.sub(r"\\label\{[^}]*\}", "", s)
    # la prosa dentro va traducida.  UN NIVEL DE ANIDAMIENTO: el `[^{}]*` de Part VI no casa
    # \text{for \emph{arbitrary} bulk content...}, y la ec. (ceiling) salia como "difiere entre
    # ediciones" cuando lo unico que difiere es la prosa traducida de dentro.
    s = re.sub(r"\\(text|mbox|mathrm)\{(?:[^{}]|\{[^{}]*\})*\}", "TEXT", s)
    s = re.sub(r"\{\\(rm|it|bf)\s+[^{}]*\}", "TEXT", s)       # ídem con {\rm palabra}
    s = re.sub(r"\\\\\[[^\]]*\]", r"\\\\", s)
    return re.sub(r"\s+", "", s)


def displays(text):
    """toda fórmula de display, en orden de aparición."""
    found = []
    for env in ENVS:
        for m in re.finditer(r"\\begin\{" + re.escape(env) + r"\}(.*?)\\end\{" + re.escape(env) + r"\}",
                             text, re.S):
            found.append((m.start(), env, m.group(1)))
    for m in re.finditer(r"(?<!\\)\\\[(.*?)(?<!\\)\\\]", text, re.S):
        found.append((m.start(), r"\[\]", m.group(1)))
    found.sort()
    return [(env, body) for _, env, body in found]


def balance(body):
    """delimitadores desequilibrados dentro de una fórmula."""
    bad = []
    depth = 0
    esc = False
    for ch in body:
        if esc:
            esc = False
            continue
        if ch == "\\":
            esc = True
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth < 0:
                bad.append("llave } de más")
                depth = 0
    if depth:
        bad.append("%d llave(s) { sin cerrar" % depth)
    if body.count(r"\left") != body.count(r"\right"):
        bad.append(r"\left %d vs \right %d" % (body.count(r"\left"), body.count(r"\right")))
    return bad


def odd_dollars(text):
    """prosa con un número impar de $ ---mata el resto del párrafo en silencio.

    POR PÁRRAFO, no por línea.  Part VI contaba línea a línea, y una fórmula en linea
    partida por el ajuste de margen ---`$c^4\\equiv c^2\\ (\\mathrm{mod}` en una linea y
    `3)$` en la siguiente, que es legal y frecuente--- daba DOS falsos positivos por cada
    corte.  El párrafo es la unidad en la que el fallo real (un `$` sin cerrar) se nota.
    """
    bad = []
    n = 1
    for para in re.split(r"\n\s*\n", text):
        clean = re.sub(r"\\\$", "", para)
        if clean.count("$") % 2:
            first = next((l for l in para.splitlines() if l.strip()), "")
            bad.append((n, first.strip()[:90]))
        n += para.count("\n") + 2
    return bad


def macros(text):
    return set(re.findall(r"\\newcommand\{\\(\w+)\}", text))


def main():
    log, fail = [], 0
    src = {}
    for tag, path in (("EN", EN), ("ES", ES)):
        src[tag] = strip_comments(path.read_text(encoding="utf-8"))

    # --- 1. las dos ediciones, fórmula a fórmula -------------------------------------
    den, des = displays(src["EN"]), displays(src["ES"])
    log.append("fórmulas de display: EN %d, ES %d" % (len(den), len(des)))
    if len(den) != len(des):
        fail += 1
        log.append("  FALLO: las ediciones no llevan el mismo número de fórmulas")
    for i, ((e1, b1), (e2, b2)) in enumerate(zip(den, des), 1):
        if e1 != e2:
            fail += 1
            log.append("  FALLO fórmula %d: entorno %s (EN) vs %s (ES)" % (i, e1, e2))
        if norm(b1) != norm(b2):
            fail += 1
            log.append("  FALLO fórmula %d (%s): difiere entre ediciones" % (i, e1))
            log.append("      EN  %s" % re.sub(r"\s+", " ", b1.strip())[:200])
            log.append("      ES  %s" % re.sub(r"\s+", " ", b2.strip())[:200])
    if fail == 0:
        log.append("  las %d fórmulas de display son idénticas en las dos ediciones" % len(den))

    # --- 2. delimitadores ------------------------------------------------------------
    for tag in ("EN", "ES"):
        n_bad = 0
        for i, (env, body) in enumerate(displays(src[tag]), 1):
            for msg in balance(body):
                fail += 1
                n_bad += 1
                log.append("  FALLO %s fórmula %d (%s): %s" % (tag, i, env, msg))
        odd = odd_dollars(src[tag])
        for n, line in odd:
            fail += 1
            n_bad += 1
            log.append("  FALLO %s línea %d: número impar de $ -- %s" % (tag, n, line))
        log.append("%s: delimitadores y $ -- %s" % (tag, "ok" if not n_bad else "%d fallos" % n_bad))

    # --- 3. macros usadas y definidas ------------------------------------------------
    men, mes = macros(src["EN"]), macros(src["ES"])
    for tag, mine, other, othertag in (("EN", men, mes, "ES"), ("ES", mes, men, "EN")):
        gone = sorted(m for m in mine if m not in other)
        if gone:
            log.append("%s define macros que %s no: %s" % (tag, othertag, ", ".join(gone)))
    for tag in ("EN", "ES"):
        defined = macros(src[tag]) | {"rep", "su", "ZZ"}
        used = set(re.findall(r"\\(\w+)", src[tag]))
        # sólo interesan las que parecen nuestras: definidas en la otra edición y no aquí
        missing = sorted(u for u in used if u in (men | mes) and u not in macros(src[tag]))
        if missing:
            fail += 1
            log.append("  FALLO %s: usa macro no definida aquí: %s" % (tag, ", ".join(missing)))
    log.append("macros: %d en EN, %d en ES" % (len(men), len(mes)))

    # --- 4. el log del compilador --------------------------------------------------
    # Un glifo que se cae no es un error de compilación: pdflatex termina en 0 y lo dice
    # sólo aquí. Así se publicó el subíndice `\sum_{\rm términos}` de la edición castellana,
    # con la tilde ausente del PDF y las siete compuertas en verde.
    for tag, path in (("EN", EN), ("ES", ES)):
        logf = path.with_suffix(".log")
        if not logf.exists():
            fail += 1
            log.append("  FALLO %s: no hay .log -- compila antes de pasar esta compuerta" % tag)
            continue
        txt = logf.read_text(encoding="utf-8", errors="replace")
        miss = re.findall(r"Missing character: There is no (.+?) in font", txt)
        math = re.findall(r"Command (\\\S+) invalid in math mode on input line (\d+)", txt)
        if miss:
            fail += 1
            log.append("  FALLO %s: %d glifo(s) caído(s) del PDF: %s" % (tag, len(miss), ", ".join(miss[:5])))
        for cmd, line in math:
            fail += 1
            log.append("  FALLO %s: %s inválido en modo matemático, línea %s" % (tag, cmd, line))
        log.append("%s log: %s" % (tag, "sin glifos caídos ni comandos inválidos en matemáticas"
                                   if not (miss or math) else "%d aviso(s)" % (len(miss) + len(math))))

    body = "\n".join(log)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(body + "\n", encoding="utf-8")
    print(body)
    print()
    print("fórmulas: TODO EN ORDEN" if not fail else "fórmulas: %d FALLO(S)" % fail)
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
