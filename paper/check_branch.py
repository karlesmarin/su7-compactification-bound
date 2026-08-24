#!/usr/bin/env python3
"""check_branch.py - toda afirmacion que depende de la semilla gauge lo dice DONDE SE DICE.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

POR QUE EXISTE.  El 23-ago-2026 el articulo se bifurco: hay dos semillas gauge posibles y unos
resultados sobreviven a las dos y otros no.  Tres repasos externos seguidos encontraron el MISMO
modo de fallo, unas cuarenta veces: **la afirmacion se corrige donde uno esta trabajando y sus
otras apariciones se quedan como estaban**.  El resumen decia la version nueva y la S7 la vieja;
la S13 explicaba la bifurcacion y la Figura 1 dibujaba la cadena sin ella.

Eso no es un fallo de razonamiento, es de barrido, y por tanto es mecanizable.  Es el tercero de
la familia: `check_numbers.py` existe porque los numeros mienten, `check_attribution.py` porque las
citas se desincronizan, y este porque **una afirmacion no vive en un solo sitio**.
[[no-claim-lives-in-one-place]]

COMO FUNCIONA.  Cada afirmacion que carga se declara UNA vez aqui, con:

  - el patron que la reconoce en el texto (regex, sobre el .tex de las dos ediciones);
  - a que rama pertenece;
  - y, si es condicional, los marcadores que la absuelven -- las palabras que, cerca, indican que
    la condicion esta dicha.

La compuerta busca TODAS las apariciones y falla si alguna condicional aparece sin marcador cerca.
No juzga la redaccion: juzga que la marca este.

LO QUE ESTA COMPUERTA NO PUEDE HACER, y conviene saberlo.  No entiende el texto.  Si una frase
afirma la cadena con otras palabras que las del patron, no la ve.  Es una red de arrastre, no un
lector.  Ampliar los patrones cuando aparezca un fallo que se le escapo es parte de usarla.
[[a-keyword-sweep-misses-the-possessive]]

Uso:  python check_branch.py     (desde part_vii/paper/)
"""
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]
P = lambda *a: print(*a, flush=True)

# Cuantos caracteres alrededor de la aparicion cuentan como "cerca".  Generoso a proposito: la
# marca puede estar en la frase anterior o en la siguiente.
WINDOW = 700

# ...PERO NO DENTRO DEL LEDGER.  El 24-ago-2026 tres filas del ledger llevaban la version de una
# sola rama y esta compuerta las daba por marcadas: en una longtable de cien lineas, 700 caracteres
# cubren media docena de filas, y el "KM25" de la fila de al lado absolvia a la que no lo tenia.
# Una fila del ledger es UNA afirmacion y su ventana es su propia fila.
# [[a-table-row-is-one-claim]]
ROWSEP = "\\\\\n"

# --------------------------------------------------------------------------- las afirmaciones
# (etiqueta, patron, marcadores que la absuelven, por que carga)
CLAIMS = [
    ("8D impar",
     # el lookbehind no es cosmetico: sin el, "no devuelve ningun 8D impar" -- el resultado
     # NEGATIVO de un barrido -- se contaba como si el articulo afirmara la imparidad.  Un falso
     # positivo en una compuerta cuesta el doble: se arregla poniendo una marca donde no hacia
     # falta, y esa marca debilita un enunciado que estaba bien.
     # el `\$?` tampoco lo es: el articulo escribe `$8D$ is odd`, con el dolar de cierre en medio,
     # y sin el la rama inglesa devolvia 0 apariciones -- una compuerta que no encuentra nada no
     # puede fallar.  Lo caza el bloque FALSIFY del final.  [[a-control-that-cannot-fail]]
     # ...y `must be` tampoco: el pie de la Figura 6 decia «$8D$ must be odd» y el patron solo
     # conocia `is`.  Los PIES son donde peor mira esta compuerta, porque sus patrones estan
     # escritos contra la prosa del cuerpo. [[a-figure-is-an-unaudited-document]]
     r"8D\}?\$?\s*(?:is|es|sea|must\s+be|ha\s+de\s+ser|debe\s+ser)\s+(?:an\s+|un\s+)?"
     r"(?:entero\s+)?(?:odd|impar)"
     r"|(?<!ningún \$)8D\$?\s+impar|odd\s+integer.{0,40}8D",
     [r"conditional", r"condicional", r"hypothesis", r"hip[oó]tesis", r"seed", r"semilla",
      r"KM25", r"sec:open", r"Corollary", r"Corolario", r"either\s+branch", r"dos\s+ramas",
      r"either\s+gauge\s+seed", r"cualquiera\s+de\s+las\s+dos"],
     "la imparidad depende de la semilla; es lo que la bifurcacion pone en duda"),

    ("|D| >= 1/8",
     r"\|D\|\s*\\?ge\\?\s*\\tfrac1?8|\|D\|\s*\\ge\s*\\tfrac\{1\}\{8\}",
     [r"conditional", r"condicional", r"seed", r"semilla", r"sec:open", r"KM25",
      r"Corollary", r"Corolario", r"hypothesis", r"hip[oó]tesis"],
     "se sigue de 8D impar, luego hereda su condicionalidad"),

    ("cota del valor zeta(5)/32",
     r"\|F\(0\)\|\s*\\;?\\ge",
     [r"seed", r"semilla", r"KM25", r"sec:open", r"printed", r"imprimen", r"HHHK",
      r"not proved", r"sin demostrar"],
     "32F(0)/zeta(5) vale 9 con una semilla y 18 con la otra: la demostracion no sobrevive"),

    ("el techo en TeV",
     r"10\.03~?\\?TeV|10\.03\$?~TeV",
     # OJO CON LOS MARCADORES SUELTOS.  Hasta el 24-ago-2026 esta lista llevaba `rama` y `branch`
     # a pelo, y las dos absolvian por accidente:
     #   - `rama` casa dentro de **prog-rama**, y este articulo dice "programa entero" en cada
     #     pagina.  El 10.03 de la PORTADA castellana ---donde la inglesa si dice "on the seed
     #     KM25 print"--- quedaba absuelto por la palabra "programa" de cuatro lineas antes;
     #   - `branch` casa en "branch point", que es el punto de ramificacion alpha^4 log alpha y
     #     no tiene nada que ver con la semilla gauge.
     # Un marcador que casa dentro de una palabra corriente no marca nada: convierte la compuerta
     # en un sello de goma.  Van con frontera de palabra y con el sentido desambiguado.
     # [[the-guard-itself-can-be-the-liar]] [[a-guard-must-not-stand-on-a-public-name]]
     [r"seed", r"semilla", r"KM25", r"sec:open", r"7\.38",
      r"(?:candidate|either|other|both|two)\s+branch(?:es)?",
      r"\bbranch(?:es)?\b(?!\s*point)(?=[^.]{0,60}\bseed\b)",
      r"(?:las\s+dos|ambas|otra|la)\s+ramas?\b(?!\s*(?:de\s+ramificaci|analítica|rota))",
      r"\bramas?\s+candidata"],
     "el techo se calcula sobre el coset, y el coset lo fija la semilla"),

    # --- los cuatro que se le escaparon el 24-ago-2026, todos del mismo tipo: la frase nombra
    # una propiedad de UNA semilla como si fuera el objeto que decide.  El objeto es chi_gauge.
    ("los colores impares",
     r"number of colours is\s+odd|n[uú]mero de colores es impar",
     [r"seed", r"semilla", r"KM25", r"sec:open", r"chigauge", r"2,\\tfrac12"],
     "con (3/2,1/2) los tres colores siguen ahi y 8D sale par: la lectura es de una semilla"),

    ("un solo semientero",
     r"single half-integer|un solo medio|[uú]nico\s+\S+\s+semientero|[uú]nico semientero"
     # el patron mira tambien los PIES, donde la frase se comprime y pierde el sustantivo: la
     # Figura 14 decia «the half-integer everything hangs on», sin «single» ni «only», y por eso
     # no casaba con nada.  [[a-figure-is-an-unaudited-document]]
     r"|only half-integer|half-integer everything hangs on|semientero (?:del|en el) que todo",
     # OJO: AQUI NO VALEN LOS MARCADORES DE SEMILLA.  Hasta el 24-ago-2026 esta lista llevaba
     # `seed`, `semilla`, `KM25` y `sec:open`, copiados de las afirmaciones de rama --- y marcan
     # OTRA COSA.  Que una frase diga de que semilla habla no dice nada sobre si el medio es lo
     # que carga la imparidad, que es justo lo que esta afirmacion pone en duda.  Con ellos, la
     # frase de la S7 «the obstruction is carried by the single half-integer coefficient» quedaba
     # absuelta por el \cite{KM25} de la demostracion de dos lineas antes, y la compuerta pasaba
     # en verde mientras la S13 decia lo contrario: «it is not the half-integer ... what the whole
     # table turns on is the parity character of the affine gauge base point».
     # Un marcador que absuelve por la razon equivocada es peor que no tener marcador.
     # [[the-absolving-half-of-a-gate-is-never-tested]]
     [r"chigauge", r"not what carries it", r"no es eso lo que la carga",
      r"tempting", r"tentador", r"\bnot\b.{0,40}same half-integer",
      r"no es el mismo semientero", r"parity of the \\emph\{sum\}",
      r"paridad\s+de la \\emph\{suma\}", r"not a single half-integer",
      r"no un semientero suelto", r"2a\)\s*\+\s*\(2b\)", r"2a\+6b",
      # y la marca de RETRACTACION, cenyida a la formula exacta con que el articulo retracta.
      # Sin ella la compuerta saltaba sobre «An earlier version of this paragraph read ... the
      # single half-integer», que no es una afirmacion sino su entierro.  Una compuerta que no
      # distingue una tesis de su retractacion obliga a marcar las retractaciones, y eso las
      # debilita.  Cenyida a proposito: no vale un «earlier version» cualquiera.
      # `\s+` y no un espacio: el .tex va a 100 columnas y la frase se parte en dos lineas.
      # Un patron con espacios literales no casa nada en cuanto el parrafo se re-fluye, y una
      # compuerta que no encuentra nada no puede fallar.  [[a-control-that-cannot-fail]]
      r"[Aa]n\s+earlier\s+version\s+of\s+this\s+paragraph\s+read",
      r"Una\s+versi[oó]n\s+anterior\s+de\s+este\s+p[aá]rrafo\s+le[ií]a"],
     "el medio no es lo que decide: lo decide chi_gauge = (2a)+(2b) mod 2"),

    ("efecto de la sexta dimension",
     r"effects? of (?:the )?(?:\\emph\{)?six(?:th)?\}?\s+dimensions?"
     r"|efecto de la[s]? (?:\\emph\{)?se(?:is|xta)\}?( dimensi[oó]n(?:es)?)?"
     r"|protection is six-dimensional|protecci[oó]n es hexadimensional",
     [r"seed", r"semilla", r"KM25", r"sec:open", r"chigauge", r"either", r"cualquiera",
      r"together with the gauge seed", r"y la semilla"],
     "seis dimensiones abren el canal antiperiodico; en que coset cae el punto base es otra cosa"),

    ("ambas columnas de su Tabla 1",
     r"columns of their Table~1|columnas de su Tabla~1",
     [r"\\emph\{our\}", r"\\emph\{nuestros\}", r"\bour\b", r"nuestros", r"reconstruction",
      r"reconstrucci[oó]n", r"sec:setting"],
     "la S2 dice que nuestra alpha_min no reproduce su columna publicada: son NUESTROS valores"),
]

# Afirmaciones que NO deben llevar marca de rama porque valen en las dos.  Si alguien les pone
# un qualifier de semilla, es un error en el otro sentido: debilita un resultado robusto.
ROBUST = [
    ("ley mod 6", r"8D\s*\\;?\\equiv\\;?\s*2A_4\+3",
     "vale con las dos semillas: el punto base da 9 = 3 (mod 6) en las dos"),
    ("2W impar", r"2W\$?\s+(?:is\s+)?(?:an\s+)?odd|2W\$?\s+es\s+(?:un\s+entero\s+)?impar",
     "el gauge aporta -3 con las dos semillas"),
]


def strip_comments(s):
    """fuera los comentarios de LaTeX: una nota nuestra no es una afirmacion del articulo."""
    return re.sub(r"(?<!\\)%.*", "", s)


def longtable_spans(body):
    """(inicio, fin) de cada longtable.  El ledger es la unica del articulo, pero no se supone."""
    spans = []
    for m in re.finditer(r"\\begin\{longtable\}", body):
        e = body.find(r"\end{longtable}", m.end())
        spans.append((m.start(), e if e >= 0 else len(body)))
    return spans


def item_spans(body):
    """(inicio, fin) de cada \\item de una lista numerada o con vinyetas.

    Un \\item es una afirmacion completa, igual que una fila del ledger.  El 24-ago-2026 el
    punto 4 de la lista "que hay de nuevo aqui" de la edicion castellana decia
    "$1/R_5\\le 10.03$~TeV para contenido arbitrario" a secas, donde la inglesa dice "for
    arbitrary content ON THE GAUGE SEED KM25 PRINT" -- y esta compuerta lo absolvia con el
    \\S\\ref{sec:open} del punto 3, tres lineas mas arriba.  La marca de un item no cubre al de
    al lado.  [[a-table-row-is-one-claim]]
    """
    spans = []
    for env in ("enumerate", "itemize"):
        for m in re.finditer(r"\\begin\{" + env + r"\}", body):
            e = body.find(r"\end{" + env + "}", m.end())
            e = e if e >= 0 else len(body)
            starts = [i.start() for i in re.finditer(r"\\item\b", body[m.end():e])]
            starts = [s + m.end() for s in starts]
            for k, s in enumerate(starts):
                spans.append((s, starts[k + 1] if k + 1 < len(starts) else e))
    return spans


def window_around(body, spans, lo_hit, hi_hit, items=()):
    """La ventana de 'cerca'.  Dentro de una longtable es LA FILA, dentro de una lista es EL ITEM;
    fuera de las dos, WINDOW caracteres.

    Una fila del ledger es una afirmacion completa con su propia casilla de 'donde'.  Si se le
    deja mirar 700 caracteres, mira las filas vecinas y se absuelve con la marca de otra.  Un
    \\item de una lista es exactamente lo mismo y tardo un mes en verse.
    """
    for a, b in spans:
        if a <= lo_hit < b:
            start = body.rfind(ROWSEP, a, lo_hit)
            start = a if start < 0 else start + len(ROWSEP)
            end = body.find(ROWSEP, hi_hit, b)
            end = b if end < 0 else end
            return body[start:end]
    for a, b in items:
        if a <= lo_hit < b:
            return body[a:b]
    return body[max(0, lo_hit - WINDOW):min(len(body), hi_hit + WINDOW)]


# --------------------------------------------------------------------------- falsacion
# Las frases REALES que el articulo llevaba el 23-ago-2026 y que esta compuerta NO cazo.  Cada
# una tiene que salir marcada como sin marca; si alguna deja de saltar, el patron se ha roto y la
# compuerta ha vuelto a ser una red con un agujero del tamanyo de un fallo que ya ocurrio.
# Un conjunto de falsacion vale mas que un control que pasa.
# [[a-falsification-suite-beats-a-passing-control]]
FALSIFY = [
    ("los colores impares",
     r"\emph{$D\neq0$ because the number of colours is odd}, and the protection is an effect "
     r"of \emph{six} dimensions rather than of the group."),
    ("los colores impares",
     r"\emph{$D\neq0$ porque el número de colores es impar}, y la protección es un efecto de "
     r"las \emph{seis} dimensiones, no del grupo."),
    ("un solo semientero",
     r"The mechanism is a single half-integer in the gauge table, and a half-integer needs a "
     r"second orbifold parity."),
    ("un solo semientero",
     r"la imparidad del sector gauge la carga \emph{un solo medio}: el peso $\tfrac12$ de cada "
     r"par adjunto con $P_6=-1$."),
    ("efecto de la sexta dimension",
     r"Both are effects of the \emph{sixth} dimension --- the same computation in five is even "
     r"at every rung."),
    ("efecto de la sexta dimension",
     r"el peso gauge semientero lo prohíbe, y ese peso es un efecto de las seis dimensiones."),
    ("ambas columnas de su Tabla 1",
     r"the same expansion gives the Higgs mass, so \emph{both} columns of their Table~1 follow "
     r"from the moments with no minimisation anywhere."),
    ("ambas columnas de su Tabla 1",
     r"así que \emph{las dos} columnas de su Tabla~1 se siguen de los momentos sin minimizar."),
    ("8D impar",
     r"Electroweak breaking needs $D>0$, and $8D$ is odd by Theorem~1, so a content can afford "
     r"the escape iff $8D\ge11$."),
    ("8D impar",
     r"Right: why one in six --- $8D$ must be odd and $8D+A_4\equiv0\pmod3$, which is "
     r"Theorem~\ref{thm:mod6}."),
    ("8D impar",
     r"A la derecha, por qué uno de cada seis: $8D$ ha de ser impar y "
     r"$8D+A_4\equiv0\pmod3$, que es el Teorema~\ref{thm:mod6}."),
    # El epigrafe de PORTADA de la edicion castellana, 24-ago-2026.  La inglesa dice "from
    # $10.03$~TeV there --- on the seed \cite{KM25} print ---"; la castellana no lo decia, y esta
    # compuerta la absolvia porque la frase de al lado contiene "programa".  Es la primera pagina.
    ("el techo en TeV",
     r"el techo por peldaño cae \emph{demostradamente}, de los $10.03$~TeV de allí a una "
     r"asíntota de Lambert-$W$ de $2.678$~TeV. Maximizar sobre el programa entero."),
]

BY_LABEL = {c[0]: c for c in CLAIMS}


def falsify():
    """Cada frase rancia tiene que dar: la caza el patron Y no la absuelve ningun marcador."""
    P("\n" + "=" * 96)
    P("FALSACION -- las frases que se colaron el 23-ago, .las caza ahora?")
    P("=" * 96)
    fails = 0
    for label, snippet in FALSIFY:
        _, pat, marks, _ = BY_LABEL[label]
        caught = re.search(pat, snippet) is not None
        absolved = any(re.search(k, snippet, re.I) for k in marks)
        ok = caught and not absolved
        fails += 0 if ok else 1
        P("  %-30s %s %s   %s"
          % (label,
             "caza" if caught else "NO CAZA",
             "sin absolver" if not absolved else "ABSUELTA POR UN MARCADOR",
             snippet[:56].replace("\n", " ")))
    if fails:
        P("\n%d frase(s) rancia(s) que esta compuerta seguiria dejando pasar." % fails)
    else:
        P("\nlas %d frases rancias saltan, y ninguna se absuelve sola." % len(FALSIFY))
    return fails


# ===========================================================================================
# LAS FIGURAS, POR DENTRO.  La tercera capa de la misma leccion, y la ultima que quedaba.
#
#   24-ago, manyana : el CUERPO decia la version condicionada y los PIES no.  Arreglados.
#   24-ago, noche   : el TEXTO DIBUJADO dentro de la figura tampoco.  Esto.
#
# Solo `fig_coset3d` ---la figura que trata DE las dos ramas--- llevaba la marca dentro. Las
# demas afirmaban en su propio dibujo:
#     fig_ceiling, fig_cone : el eje rotulado «$8D$ (odd, by the theorem)»
#     fig_lift              : la leyenda «admissible: $8D$ odd, ...»
#     fig_chi2              : «the best any bulk content can do»
#     fig_cone              : «$8D$ odd and $8D+A_4\equiv0$ (mod 3)», juntando en una frase
#                             lo que depende de la semilla con lo que no
#
# Un eje es donde se dice la afirmacion tanto como un pie, y el lector lo mira ANTES. Lo que
# no cabe ahi es una subordinada; lo que si cabe es una referencia.
#
# Y LOS PATRONES SON OTROS.  Dentro de una figura no hay verbo: el cuerpo escribe «$8D$ is
# odd» y el eje escribe «8D odd» a secas, de modo que los patrones del .tex ---que piden
# is/es/must be/ha de ser--- no encuentran nada. Una compuerta que no encuentra nada no puede
# fallar. [[a-control-that-cannot-fail]] [[a-figure-is-an-unaudited-document]]
FIG_CLAIMS = [
    # `[\s(]+` y no `\s+`: el eje de fig_surface rotulaba «$8D$ (odd, log scale)» y el
    # parentesis rompia la adyacencia, de modo que el barrido no lo veia.  Comprimir un rotulo
    # quita justo las palabras ---y ahora tambien los espacios--- por los que la red busca.
    ("8D impar, dibujado", r"8\s*D\s*[\s(]\s*(?:is\s+)?odd|8\s*D\s*[\s(]\s*impar"
                           r"|odd,\s*by the theorem|impar,\s*por el teorema"),
    # cenyido a la forma que AFIRMA, no a la que describe.  La primera version casaba tambien
    # el pie de fig_chi2, «cada 1/R5 que un contenido del bulk puede alcanzar», que dice que el
    # panel dibuja el conjunto alcanzable ---no que ningun contenido pueda superarlo--- y el
    # marcador estaba en la anotacion de al lado, a mas de setenta caracteres.  Ensanchar la
    # ventana habria sido el arreglo comodo y el equivocado: dentro de una figura de tres
    # paneles, una marca en el primero no condiciona una afirmacion del tercero.
    ("todo contenido, dibujado", r"best any bulk content can do"
                                 r"|lo mejor que puede hacer un contenido"),
]
# Dentro de una figura la cita sale numerada, no como \cite{KM25}: los marcadores son otros.
FIG_MARCAS = [r"\bseed\b", r"semilla", r"published seed", r"candidate seed", r"\[2\]",
              r"Komori", r"conditional", r"condicional"]


def figuras_por_dentro():
    """Devuelve el numero de afirmaciones de rama SIN MARCA dibujadas dentro de una figura."""
    import glob
    try:
        import fitz
    except ImportError:
        P("\n  *** PyMuPDF no esta: el texto de las figuras no se ha mirado. ***")
        P("  Una comprobacion que no corre no es una comprobacion que pasa.")
        return 1

    P("")
    P("=" * 96)
    P("EL TEXTO DIBUJADO DENTRO DE LAS FIGURAS -- donde el lector mira antes que al pie")
    P("=" * 96)
    pdfs = sorted(glob.glob(os.path.join(HERE, "fig_*.pdf")))
    if not pdfs:
        P("  *** no hay ninguna fig_*.pdf: no hay nada que mirar, y eso es un fallo ***")
        return 1
    malas = 0
    mirados = 0
    for p in pdfs:
        nombre = os.path.basename(p)
        try:
            txt = " ".join(" ".join(pg.get_text() for pg in fitz.open(p)).split())
        except Exception as e:                       # noqa: BLE001
            P("  %-26s *** no se puede leer: %s ***" % (nombre, e))
            malas += 1
            continue
        mirados += 1
        for etiqueta, pat in FIG_CLAIMS:
            for m in re.finditer(pat, txt, re.I):
                frag = txt[max(0, m.start() - 70):m.end() + 70]
                if not any(re.search(k, frag, re.I) for k in FIG_MARCAS):
                    malas += 1
                    P("  %-26s %-22s SIN MARCA" % (nombre, etiqueta))
                    P("        ...%s..." % frag[:84])
    P("  %d figuras leidas, %d afirmacion(es) de rama sin marca dentro del dibujo"
      % (mirados, malas))

    # ...Y EL NUMERO DE LA CITA, que las figuras llevan A MANO.  Dentro de un dibujo no se
    # puede escribir \cite{KM25}: se escribe «[2]».  Ese 2 lo decide el orden de la
    # bibliografia, de modo que anyadir una referencia antes de Komori-Maru convertiria todas
    # esas marcas en una cita a otro articulo, en silencio y en las dos ediciones.
    # [[cross-document-citations-are-unguarded]]
    P("")
    P("  el numero con que las figuras citan a Komori-Maru, contra el .aux de cada edicion:")
    for tex in TEXS:
        aux = os.path.join(HERE, tex.replace(".tex", ".aux"))
        if not os.path.exists(aux):
            P("    %-24s *** no hay .aux: compila antes ***" % tex)
            malas += 1
            continue
        m = re.search(r"\\bibcite\{KM25\}\{(\d+)\}",
                      open(aux, encoding="utf-8", errors="replace").read())
        n = m.group(1) if m else None
        ok = (n == "2")
        P("    %-24s \\cite{KM25} -> [%s]   %s"
          % (tex, n, "ok, es el que dibujan las figuras" if ok
             else "*** las figuras dibujan [2] y esto NO es 2 ***"))
        malas += 0 if ok else 1
    return malas


# --------------------------------------------------------------------------- falsacion (figuras)
# Los rotulos REALES que llevaban las figuras el 24-ago-2026, antes de condicionarlos.
FIG_FALSIFY = [
    ("el eje de fig_ceiling y fig_cone", "8D (odd, by the theorem) 3 4 5 6 7 8 9 10"),
    ("la leyenda de fig_lift", "admissible: 8D odd, A4 + 8D 0 (mod 3)"),
    ("la anotacion de fig_chi2", "the best any bulk content can do 10.03 TeV"),
]


def falsify_figuras():
    P("")
    P("=" * 96)
    P("FALSACION -- los rotulos que las figuras llevaban dibujados, .los ve?")
    P("=" * 96)
    fails = 0
    for nombre, texto in FIG_FALSIFY:
        caza = any(re.search(pat, texto, re.I) for _e, pat in FIG_CLAIMS)
        absuelto = any(re.search(k, texto, re.I) for k in FIG_MARCAS)
        ok = caza and not absuelto
        fails += 0 if ok else 1
        P("  %-36s %s" % (nombre, "caza sin absolver" if ok
                          else ("NO LA CAZA" if not caza else "LA ABSUELVE SOLA")))
    if fails:
        P("\n%d rotulo(s) real(es) que este bloque dejaria pasar." % fails)
    else:
        P("\nlos tres rotulos rancios saltan.")
    return fails


def main():
    bad = []
    P("=" * 96)
    P("TODA AFIRMACION QUE DEPENDE DE LA SEMILLA GAUGE, .LO DICE DONDE SE DICE?")
    P("=" * 96)
    for tex in TEXS:
        path = os.path.join(HERE, tex)
        if not os.path.exists(path):
            P("  FALTA %s" % tex)
            bad.append((tex, "no existe"))
            continue
        body = strip_comments(open(path, encoding="utf-8", errors="replace").read())
        body = body.split(r"\begin{thebibliography}")[0]
        spans = longtable_spans(body)
        items = item_spans(body)
        P("\n%s" % tex)
        for label, pat, marks, why in CLAIMS:
            hits = list(re.finditer(pat, body))
            unmarked = []
            for m in hits:
                near = window_around(body, spans, m.start(), m.end(), items)
                if not any(re.search(k, near, re.I) for k in marks):
                    ln = body[:m.start()].count("\n") + 1
                    unmarked.append(ln)
            flag = "sin marca en %d" % len(unmarked) if unmarked else "todas marcadas"
            P("  %-26s %3d aparicion(es)   %s" % (label, len(hits), flag))
            if unmarked:
                for ln in unmarked[:6]:
                    frag = body.split("\n")[ln - 1].strip()[:88]
                    P("        linea %-5d %s" % (ln, frag))
                bad.append((tex, label, unmarked))
        for label, pat, why in ROBUST:
            n = len(re.findall(pat, body))
            P("  %-26s %3d aparicion(es)   robusta, no necesita marca" % (label, n))

    P("\n" + "=" * 96)
    if bad:
        P("AFIRMACIONES CONDICIONADAS SIN DECIRLO: %d" % len(bad))
        for b in bad:
            P("   %s" % (b,))
        P("")
        P("Una afirmacion no vive en un solo sitio.  Si se corrige donde uno esta trabajando y no")
        P("donde ademas aparece, el articulo dice dos cosas distintas segun por donde se abra.")
        P("=" * 96)
        falsify()
        figuras_por_dentro()
        falsify_figuras()
        return 1
    P("cada afirmacion condicionada lleva su marca cerca, en las dos ediciones.")
    P("=" * 96)
    mal_fig = figuras_por_dentro() + falsify_figuras()
    return 1 if (falsify() or mal_fig) else 0


if __name__ == "__main__":
    raise SystemExit(main())
