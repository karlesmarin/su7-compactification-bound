#!/usr/bin/env python3
"""check_layout.py - no page of the SU(7) paper carries an internal blank band.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

A tall unbreakable box (the ledger, the verification table, a figure) can push a
page break early and leave a hole in the MIDDLE of a page.  Slack at the FOOT is
what \\raggedbottom is for and is fine; slack in the middle is not.  This measures
the largest gap between consecutive pieces of inked content on every page.

Run:  python check_layout.py     (from part_vi/paper/)
"""
import re
import sys

import fitz

# La compuerta imprime el TEXTO del bloque que se sale, y ese texto puede llevar cualquier glifo
# del PDF.  Sin esto, un '<=' en la linea desbordada mata al guion con UnicodeEncodeError y la
# compuerta muere en vez de fallar: sale un traceback donde tendria que salir un veredicto.
# Paso el 2026-08-22, con el recuadro del techo riguroso.  Ver [[edits-must-fail-loudly]].
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PDFS = ["su7_hierarchy.pdf", "su7_hierarchy_es.pdf"]
LIMIT = 120.0          # pt; a band larger than this is worth looking at
# The page number sits in the footer, below the text block (A4 at margin 2.4cm
# ends at ~774pt).  Counting it as content makes the gap between the last line
# and the folio look like a band, which flags the LAST page of every document --
# a guard that always fires measures nothing.  So the footer is not content.
FOOTER = 780.0
# ...but excluding the footer must not blind the guard to slack ABOVE it. On any
# page but the last, a large gap between the last line and the bottom of the text
# block means an unbreakable box (a table, a figure) did not fit and was pushed
# over, usually leaving a stranded section heading behind. That is a defect and
# the first version of this file could not see it.
TEXT_BOTTOM = 774.0        # A4 at margin 2.4cm
FOOT_LIMIT = 150.0
# ...and the footer filter above had a second, worse blind spot, found by Carles reading the
# PDF on 2026-08-20: it drops EVERYTHING below y = FOOTER, so ink that runs off the BOTTOM of
# the page is dropped too. Figure 1's caption was doing exactly that, in both editions -- the
# last line sat at y = 848.5 on a 841.9 pt page, i.e. printed past the paper edge, and the
# gate said "both editions clean". A float page can overrun without LaTeX warning at all.
# The folio sits at ~785, so anything below BLEED is not a folio: it is spill.
BLEED = 35.0               # pt of true bottom margin that must stay free of ink
# ...and a fourth blind spot, found by Carles reading the Spanish PDF on 2026-08-21: everything
# above measures the page VERTICALLY.  Eq. (17)'s boxed line ran off the RIGHT of its frame -- the
# border cut through the D of D = 1/8 -- and this file said "both editions clean".  A translated
# line is routinely ten characters longer than its original, so horizontal overflow is a failure
# mode of the ES edition specifically, and nothing was watching it.  A4 at margin 2.4cm gives a
# text block from 68 to 527 pt; RIGHT is where ink stops being text and starts being spill.
LEFT_EDGE, RIGHT_EDGE = 68.0, 528.0
SIDE_TOL = 2.0             # pt of slack for glyph bearings and rules that sit on the margin

ALLBAD = []
for PDF in PDFS:
 doc = fitz.open(PDF)
 print("%s: %d pages" % (PDF, doc.page_count))
 worst = []
 for i, page in enumerate(doc, 1):
     h = page.rect.height
     rows = set()
     for b in page.get_text("blocks"):
         y0, y1 = b[1], b[3]
         for y in range(int(y0), int(y1) + 1):
             rows.add(y)
     for d in page.get_drawings():
         r = d["rect"]
         for y in range(int(r.y0), int(r.y1) + 1):
             rows.add(y)
     for b in page.get_images(full=True):
         for r in page.get_image_rects(b[0]):
             for y in range(int(r.y0), int(r.y1) + 1):
                 rows.add(y)
     # --- spill off the SIDES.  Text blocks only: get_drawings() returns the page rules and
     #     the figure frames, which legitimately touch the margin, so they are not evidence.
     for b in page.get_text("blocks"):
         x0, x1, txt = b[0], b[2], b[4]
         if not txt.strip():
             continue
         if x1 > RIGHT_EDGE + SIDE_TOL or x0 < LEFT_EDGE - SIDE_TOL:
             over = max(x1 - RIGHT_EDGE, LEFT_EDGE - x0)
             print("  p%-3d TEXT RUNS OFF THE SIDE by %5.1f pt  [%.0f..%.0f]  %r"
                   % (i, over, x0, x1, " ".join(txt.split())[:60]))
             ALLBAD.append((PDF, i, "side overflow %.1f pt" % over))

     # --- spill off the bottom of the sheet, measured BEFORE the footer filter
     if rows and max(rows) > h - BLEED:
         print("  p%-3d ink runs to %6.1f on a %6.1f pt page   SPILLS OFF THE BOTTOM by "
               "%6.1f pt" % (i, max(rows), h, max(rows) - (h - BLEED)))
         ALLBAD.append((PDF, i, "spill %.1f pt past the bottom margin" % (max(rows) - (h - BLEED))))

     rows = {y for y in rows if y < FOOTER}
     if rows and i < doc.page_count:
         slack = TEXT_BOTTOM - max(rows)
         if slack > FOOT_LIMIT:
             # ...but slack at the foot is only a defect if something was STRANDED. If the very
             # next page opens with a full-width float, the slack is just that float not fitting
             # in what was left, which is what \raggedbottom is for and reads correctly. This is
             # measured, not assumed: it asks whether the next page carries an image high up.
             # A \includegraphics of a PDF is a Form XObject, so get_images() is blind to it;
             # what is always there is the caption. A caption in the top half of the next page
             # means that page opens with the float.
             nxt = doc.load_page(i)          # 0-based: the page after this one
             floated = any(re.match(r"\s*(Figure|Figura|Table|Cuadro|Tabla)\s+\d+", b[4])
                           and b[1] < TEXT_BOTTOM * 0.55
                           for b in nxt.get_text("blocks"))
             # ...and the same question for an unbreakable INLINE box, which the float test
             # above cannot see because it has no caption. The branch table of section 13 is
             # 240 pt of \footnotesize tabular; when 164 pt are left, it cannot fit, and no
             # amount of prose editing changes that. So the cause is MEASURED rather than
             # assumed: take the first contiguous run of ink on the next page, and if that run
             # is taller than the slack, the slack is arithmetic and not a stranded page.
             # This does not blind the guard -- a stranded heading or an orphaned paragraph is
             # a few lines tall, far under FOOT_LIMIT, and still fails. [[a-deny-needs-a-confirmed-cause]]
             # El discriminante es la ALTURA DE CADA BLOQUE, no el hueco: un parrafo de prosa es
             # UN bloque alto (132 pt en la p.57 castellana) y una tabla son MUCHOS bloques bajos
             # que ademas se solapan verticalmente.  De modo que la racha de arriba se extiende
             # mientras los bloques sean bajos, y para en cuanto empieza la prosa.  Un titulo
             # huerfano o un parrafo abandonado no producen racha, y siguen fallando.
             tall = 0.0
             if not floated:
                 bl = sorted([b for b in nxt.get_text("blocks")
                              if b[4].strip() and b[3] < FOOTER], key=lambda b: b[1])
                 # el PRIMER bloque tambien tiene que ser bajo: sin esta comprobacion, un parrafo
                 # abandonado -- un solo bloque de 600 pt -- se absolvia solo.  Lo caza el bloque
                 # de falsacion del final.
                 if bl and (bl[0][3] - bl[0][1]) <= 70.0:
                     top, end = bl[0][1], bl[0][3]
                     for b in bl[1:]:
                         if b[1] > end + 30.0 or (b[3] - b[1]) > 70.0:
                             break
                         end = max(end, b[3])
                     tall = end - top
             if floated:
                 print("  p%-3d ink %6.1f..%6.1f   foot slack %6.1f pt, explained: p%d opens "
                       "with a float" % (i, min(rows), max(rows), slack, i + 1))
             elif tall > slack:
                 print("  p%-3d ink %6.1f..%6.1f   foot slack %6.1f pt, explained: p%d opens "
                       "with an unbreakable %.0f pt box that could not fit in it"
                       % (i, min(rows), max(rows), slack, i + 1, tall))
             else:
                 print("  p%-3d ink %6.1f..%6.1f   FOOT SLACK %6.1f pt  <-- something was "
                       "pushed to the next page" % (i, min(rows), max(rows), slack))
                 worst.append((i, slack, 0, min(rows), max(rows)))
             continue
     if not rows:
         print("  p%-3d EMPTY PAGE" % i)
         worst.append((i, h))
         continue
     ys = sorted(rows)
     gap, at = 0.0, 0
     for a, b in zip(ys, ys[1:]):
         if b - a > gap:
             gap, at = b - a, a
     worst.append((i, gap, at, ys[0], ys[-1]))
     flag = "  <-- internal band" if gap > LIMIT else ""
     print("  p%-3d ink %6.1f..%6.1f   largest internal gap %6.1f pt at y=%d%s"
           % (i, ys[0], ys[-1], gap, at, flag))
 bad = [w for w in worst if len(w) > 2 and w[1] > LIMIT]
 print(" %s: pages flagged (internal band over %.0f pt, or foot slack over %.0f pt): %d"
       % (PDF, LIMIT, FOOT_LIMIT, len(bad)))
 ALLBAD += [(PDF, b) for b in bad]
print()

# ------------------------------------------------------------------ falsacion de la absolucion
# La regla que absuelve un hueco al pie ("lo que sigue no cabia") es la unica de este fichero que
# puede CALLAR un fallo, asi que se la pone a prueba sobre casos sinteticos en cada corrida.  El
# tercero -- un parrafo abandonado, un solo bloque altisimo -- se absolvia solo hasta que este
# bloque lo caz'o.  [[a-falsification-suite-beats-a-passing-control]]
def _unbreakable(blocks, slack):
    bl = sorted(blocks)
    if not bl or (bl[0][1] - bl[0][0]) > 70.0:
        return False
    top, end = bl[0]
    for y0, y1 in bl[1:]:
        if y0 > end + 30.0 or (y1 - y0) > 70.0:
            break
        end = max(end, y1)
    return (end - top) > slack


FALSIFY = [
    ("stranded heading then prose", [(73, 90), (100, 600)], 160.0, False),
    ("orphaned display then prose", [(73, 98), (110, 620)], 160.0, False),
    ("abandoned paragraph", [(73, 700)], 160.0, False),
    ("15-row table that cannot fit", [(73, 82), (88, 100), (100, 132), (122, 165),
                                      (160, 224), (230, 274), (294, 427)], 164.0, True),
    ("short table that WOULD have fit", [(73, 82), (88, 100), (100, 120), (294, 427)],
     160.0, False),
]
bad_rule = 0
print("falsification of the foot-slack absolution:")
for label, blocks, slack, want in FALSIFY:
    got = _unbreakable(blocks, slack)
    ok = got == want
    bad_rule += 0 if ok else 1
    print("  %-34s absolves=%-5s expected=%-5s  %s" % (label, got, want, "ok" if ok else "BROKEN"))
if bad_rule:
    print("  the rule that can silence a failure is itself broken: %d case(s)" % bad_rule)
print()

# =========================================================================================
# EL DESBORDE HORIZONTAL, que hasta el 24-ago-2026 no miraba NADIE.
#
# Esta compuerta media bandas VERTICALES ---huecos en medio de una pagina, tinta por debajo del
# borde--- y no miraba ni una sola vez a la derecha.  Ninguna de las otras doce lo hacia tampoco:
# `check_formulas` lee los .log pero solo busca glifos caidos y comandos invalidos.
#
# El 24-ago-2026, al arreglar la caja de |F(0)| para que llevara su condicion de semilla, la
# linea castellana crecio 20 pt y se salio del marco: la palabra «imprimen» cruzaba el borde
# azul y el numero de ecuacion caia a una segunda linea.  Las trece compuertas pasaron en verde.
# Lo vio el PDF renderizado.  [[the-pdf-text-layer-is-not-the-page]] [[a-figure-is-an-unaudited-document]]
#
# Se mide de DOS formas, porque una sola no basta:
#
#   (a) el .log, que dice si algo es mas ancho que su caja.  Es necesario y no suficiente: un
#       Overfull de 7 pt dentro de un `keyeq` lo absorben los 9 pt de `\fboxsep` y no se ve.
#   (b) el PDF, que dice si la tinta CRUZA el marco de verdad.  Es el control que decide, y
#       puede fallar: si alguien ensancha el texto de una caja, salta.
#
# La (a) sola convertiria un defecto invisible en un fallo, y una lista blanca ancha acabaria
# tapando tambien el visible.  La (b) sola no avisaria de un desborde en prosa corriente.
#
# FALSACION, corrida el 24-ago-2026 y no supuesta.  Se devolvio la caja castellana a la version
# rota ---«para todo contenido DE BULK, y se alcanza...»---, se compilo, y las dos mitades
# saltaron:
#     (a)  el .log: 20.05 pt, fuera de la banda aceptada de 7-8 pt   -> senyalado
#     (b)  el PDF : 10.6 pt de tinta pasado el marco, pagina 19      -> senyalado
# Luego se restauro. Si alguna de las dos deja de saltar sobre ese caso, esta comprobacion ha
# vuelto a ser un adorno.  [[a-control-that-cannot-fail]] [[a-falsification-suite-beats-a-passing-control]]
OVERFULL = re.compile(r"Overfull \\hbox \((\d+(?:\.\d+)?)pt too wide\).*?at line (\d+)")

# Desbordes comprobados A MANO sobre el PDF renderizado, con su razon escrita.  Una lista blanca
# sin razones deja de ser una lista blanca.  [[a-control-that-cannot-fail]]
OVERFULL_ACEPTADOS = [
    # La caja `eq:resum`, en las dos ediciones.  El `cases` es ~7 pt mas ancho que el minipage
    # de 0.88\textwidth, pero los 9 pt de \fboxsep lo absorben: medido sobre el PDF, la tinta
    # mas a la derecha se queda a 0.0 pt del marco en las dos ediciones, no lo cruza.  Estrecho,
    # no roto.  Si algun dia crece, lo cazara la comprobacion (b) y no esta.
    (7.0, 8.0, "la caja eq:resum, absorbida por el \\fboxsep de 9 pt; verificado sobre el PDF"),
]


def _overfull_del_log(path):
    try:
        txt = open(path, encoding="utf-8", errors="replace").read()
    except OSError:
        return None
    return [(float(a), int(b)) for a, b in OVERFULL.findall(txt)]


def _tinta_fuera_del_marco(pdf):
    """Todo marco keyeq del PDF, y cuanto se sale de el la tinta que contiene."""
    doc = fitz.open(pdf)
    peor = []
    for i, page in enumerate(doc):
        marcos = [d["rect"] for d in page.get_drawings()
                  if d["rect"].width > 380 and d["rect"].height < 260]
        if not marcos:
            continue
        spans = [fitz.Rect(sp["bbox"])
                 for blk in page.get_text("dict")["blocks"]
                 for ln in blk.get("lines", []) for sp in ln["spans"]]
        for R in marcos:
            fuera = max((b.x1 - R.x1 for b in spans
                         if R.y0 - 2 < b.y0 and b.y1 < R.y1 + 2), default=0.0)
            if fuera > 0.5:
                peor.append((i + 1, round(fuera, 1)))
    return peor


print()
print("=" * 84)
print("DESBORDE HORIZONTAL -- el .log dice si algo es mas ancho que su caja, el PDF si se ve")
print("=" * 84)
bad_wide = 0
for pdf in PDFS:
    log = pdf.replace(".pdf", ".log")
    hits = _overfull_del_log(log)
    if hits is None:
        print("  %-24s *** no hay %s: compila antes ***" % (pdf, log))
        bad_wide += 1
        continue
    sin_razon = []
    for pt, line in hits:
        if not any(lo <= pt <= hi for lo, hi, _ in OVERFULL_ACEPTADOS):
            sin_razon.append((pt, line))
    print("  %-24s %d Overfull en el log, %d sin razon escrita"
          % (log, len(hits), len(sin_razon)))
    for pt, line in sin_razon:
        print("        %6.2f pt de mas, linea %d del .tex" % (pt, line))
    bad_wide += len(sin_razon)

    fuera = _tinta_fuera_del_marco(pdf)
    print("  %-24s %s" % (pdf,
                          "ninguna caja enmarcada tiene tinta fuera del marco" if not fuera
                          else "*** TINTA FUERA DEL MARCO: %s ***" % fuera))
    bad_wide += len(fuera)

if bad_wide:
    print()
    print("  Un desborde de pocos puntos puede ser inocuo y uno de veinte cruza el marco.")
    print("  La diferencia se mide sobre el PDF, no se supone.")
print()

print("both editions clean" if not (ALLBAD or bad_rule or bad_wide)
      else "FLAGGED: %s" % (ALLBAD or "rule broken" if not bad_wide else "horizontal overflow"))
raise SystemExit(1 if (ALLBAD or bad_rule or bad_wide) else 0)
