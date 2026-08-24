#!/usr/bin/env python3
"""gate_km_ghosts.py -- la ec. (68) de Komori-Maru, auditada contra SU PROPIO texto.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

POR QUE ESTE GUION.  Todo el Teorema 1 de este articulo (8D impar) se come un numero ajeno: el
-27 que sale de los coeficientes de la ec. (68) de Komori-Maru, arXiv:2503.04090.  Nosotros no lo
calculamos: lo leemos.  gauge_ghost_seed.py ya midio QUE PASA si ese numero fuera otro.  Este
guion mide algo distinto y anterior: SI HAY MOTIVO para pensar que lo es, y el motivo tiene que
salir del articulo de ellos, no de nuestra aritmetica ni de una consulta.  [[audit-the-formula-not-the-words]]

EL ARGUMENTO, EN CUATRO HECHOS DE SU TEXTO
  1. su ec. (49) define N como los grados de libertad FISICOS que corren en el bucle;
  2. su ec. (51) anade un termino de gauge-fixing -- luego trabajan en gauge covariante;
  3. en gauge covariante el determinante de Faddeev-Popov es obligatorio, y las palabras
     "ghost" y "Faddeev" no aparecen NI UNA VEZ en las 42 paginas;
  4. ponen N = 4 para A_mu y N = 1 para A_5,6 -- y 4 + 1 = 5, mientras que en seis dimensiones
     un campo gauge tiene D - 2 = 4 grados fisicos, no cinco.

El 4 es el numero de COMPONENTES de A_mu, no de grados fisicos: un vector 4D sin masa tiene 2 y
con masa tiene 3, nunca 4.  Y el 1 de A_5,6 SI esta contado despues de quitar el Goldstone -- su
propia nota al pie lo dice.  Las dos contabilidades consistentes dan cuatro:

     covariante :  4 (componentes de A_mu) - 2 (fantasmas) + 1 (A_5) + 1 (A_6)  =  4
     fisica     :  3 (vector masivo)                       + 1 (escalar ortog.) =  4

La suya mezcla una mitad de cada una.  Y Haba-Yamashita, en la formula general 5D que este mismo
articulo hereda (hep-ph/0401185), escriben el coeficiente gauge+fantasma como -(D-2) explicito.

QUE NO DEMUESTRA ESTE GUION.  Que ellos se equivoquen.  Podria haber una convencion o un
jacobiano no explicitado que lo arregle.  Lo que demuestra es que la suma de sus dos N no cumple
una condicion NECESARIA que su propio marco impone, y que la unica pieza capaz de explicar la
diferencia es la que no aparece nombrada en ningun sitio.  [[a-necessary-condition-not-a-convenient-one]]

CONTROLES
  C1  el PDF se lee y es el articulo que creemos (control de instrumento: si abrimos otro
      fichero, todo lo demas es teatro).  [[falsify-the-instrument-first]]
  C2  su ec. (49) define N como grados FISICOS -- cita literal con pagina.
  C3  su ec. (51) es un termino de gauge-fixing -- cita literal con pagina.
  C4  "ghost" y "Faddeev" aparecen cero veces.  Se buscan tambien variantes y el plural, porque
      un barrido por una sola palabra se deja la mitad.  [[a-keyword-sweep-misses-the-possessive]]
  C5  las dos asignaciones N = 4 y N = 1, literales y con pagina.
  C6  la regla de suma: 4 + 1 = 5 contra D - 2 = 4, y las dos contabilidades consistentes.
  C7  Haba-Yamashita escriben -(D-2) en hep-ph/0401185 -- cita literal.
  C8  la salida NO es una normalizacion global: un reescalado uniforme movería los tres
      coeficientes, y el cambio 4 -> 3 mueve solo el canal periodico.  Esto cierra la unica
      escapatoria de "sera una convencion".

Run:  python gate_km_ghosts.py > outputs/gate_km_ghosts.txt
"""
import pathlib
import re
import sys
from fractions import Fraction as F

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent
PAPERS = HERE.parent / "_papers"
KM = PAPERS / "SU7_GGHU_2503.04090.pdf"          # Komori-Maru, SU(7) Grand Gauge-Higgs Unification
HY = PAPERS / "hy04_hep-ph_0401185.pdf"          # Haba-Yamashita, general 5D formula
AH = PAPERS / "AHMN_2312.08608.pdf"              # Akamatsu-Hirose-Maru-Nago, 6D SU(4) en T^2/Z_2
HN = PAPERS / "hnt_hep-ph_0403106.pdf"           # Hosotani-Noda-Takenaga, 6D en el MISMO T^2/Z_2

try:
    import fitz  # PyMuPDF
except ImportError:
    raise SystemExit("FATAL: hace falta PyMuPDF (pip install pymupdf).  pdftotext NO sirve aqui: "
                     "parte las formulas y se come glifos.  [[pdftotext-lies-use-pymupdf]]")


def pages_of(path):
    """lista de (numero de pagina, texto).  1-indexado como en el PDF."""
    if not path.exists():
        raise SystemExit("FATAL: no existe %s" % path)
    with fitz.open(path) as doc:
        return [(i + 1, pg.get_text()) for i, pg in enumerate(doc)]


def find(pages, needle, flags=re.I):
    """[(pagina, linea entera)] donde aparece needle."""
    hits = []
    for n, txt in pages:
        for line in txt.splitlines():
            if re.search(needle, line, flags):
                hits.append((n, line.strip()))
    return hits


def quote(hits, limit=3):
    for n, line in hits[:limit]:
        P("         p.%-3d  %s" % (n, line[:110]))
    if len(hits) > limit:
        P("         ... y %d mas" % (len(hits) - limit))


def main():
    fails = []
    P("=" * 98)
    P("LA EC. (68) DE KOMORI-MARU, AUDITADA CONTRA SU PROPIO TEXTO")
    P("=" * 98)

    # ------------------------------------------------------------------ C1
    P("\n[C1] CONTROL DE INSTRUMENTO: ABRIMOS LOS ARTICULOS QUE CREEMOS")
    km_pages = pages_of(KM)
    hy_pages = pages_of(HY)
    km_head = km_pages[0][1][:400].replace("\n", " ")
    hy_head = hy_pages[0][1][:400].replace("\n", " ")
    ok_km = "Grand Gauge-Higgs" in km_head or "SU(7)" in km_head
    ok_hy = "general formula" in hy_head.lower() and "5D" in hy_head
    P("     KM  %-42s  %d paginas" % (KM.name, len(km_pages)))
    P("         %s" % km_head[:150])
    P("     HY  %-42s  %d paginas" % (HY.name, len(hy_pages)))
    P("         %s" % hy_head[:150])
    P("\n     ojo: hep-ph/0401185 es HABA-YAMASHITA (5D), NO Komori-Maru.  La ec. (68) es de")
    P("     arXiv:2503.04090.  Confundirlas manda al revisor al papel equivocado.")
    ok1 = ok_km and ok_hy
    P("     C1 %s" % ("PASS" if ok1 else "FAIL -- alguno de los dos no es lo que creemos"))
    if not ok1:
        fails.append("C1")
        P("\nABORTADO EN C1.")
        return 1

    # ------------------------------------------------------------------ C2
    P("\n[C2] SU EC. (49) DEFINE N COMO GRADOS DE LIBERTAD *FISICOS*")
    h = find(km_pages, r"physical degree of freedom for the fields")
    quote(h)
    ok2 = bool(h)
    P("     C2 %s -- %s" % ("PASS" if ok2 else "FAIL",
                            "N es, por su propia definicion, un conteo FISICO" if ok2
                            else "no se encuentra la definicion"))
    if not ok2:
        fails.append("C2")

    # ------------------------------------------------------------------ C3
    P("\n[C3] SU EC. (51) ES UN TERMINO DE GAUGE-FIXING -> GAUGE COVARIANTE")
    h = find(km_pages, r"gauge.fixing term|L\s*GF|Adding the following gauge")
    quote(h)
    ok3 = bool(h)
    P("     C3 %s -- %s" % ("PASS" if ok3 else "FAIL",
                            "fijan el gauge, luego el determinante de FP es obligatorio" if ok3
                            else "no se encuentra el termino de gauge-fixing"))
    if not ok3:
        fails.append("C3")

    # ------------------------------------------------------------------ C4
    P("\n[C4] CUANTAS VECES NOMBRAN AL FANTASMA")
    WORDS = ["ghost", "ghosts", "Faddeev", "Popov", "FP determinant", "BRST", "anti-?ghost"]
    total = 0
    for w in WORDS:
        n = sum(len(re.findall(w, t, re.I)) for _, t in km_pages)
        total += n
        P("       %-16s %d" % (w, n))
    ok4 = total == 0
    P("\n     C4 %s -- en %d paginas, %d menciones.  Fijan el gauge y no introducen el"
      % ("PASS" if ok4 else "FAIL", len(km_pages), total))
    P("        determinante que fijar el gauge exige.")
    if not ok4:
        fails.append("C4")

    # ------------------------------------------------------------------ C5
    P("\n[C5] LAS DOS ASIGNACIONES, LITERALES")
    h4 = find(km_pages, r"setting N = 4")
    h1 = find(km_pages, r"setting N = 1")
    hf = find(km_pages, r"physical degrees of freedom to be N = 1 not N = 2")
    quote(h4, 2)
    quote(h1, 2)
    P("     y la nota al pie que explica el 1:")
    quote(hf, 2)
    h_abs = find(km_pages, r"absorbed as the longitudinal mode")
    quote(h_abs, 2)
    ok5 = bool(h4) and bool(h1)
    P("\n     C5 %s -- A_mu lleva N = 4 ; A_5,6 lleva N = 1, y el 1 esta contado DESPUES de"
      % ("PASS" if ok5 else "FAIL"))
    P("        quitar el Goldstone, por su propia nota.  El 4 no esta contado despues de nada.")
    if not ok5:
        fails.append("C5")

    # ------------------------------------------------------------------ C6
    P("\n[C6] LA REGLA DE SUMA, QUE ES UNA CONDICION NECESARIA DE SU PROPIO MARCO")
    D = 6
    P("       en D = %d dimensiones un campo gauge tiene D - 2 = %d grados fisicos" % (D, D - 2))
    P("")
    P("       %-46s %s" % ("su conteo            :  N(A_mu) + N(A_5,6) = 4 + 1", "= 5"))
    P("       %-46s %s" % ("covariante consistente:  4 - 2 + 1 + 1", "= 4"))
    P("       %-46s %s" % ("fisico consistente    :  3 (vector masivo) + 1", "= 4"))
    ok6 = (4 + 1) != (D - 2) and (4 - 2 + 1 + 1) == (D - 2) and (3 + 1) == (D - 2)
    P("\n     C6 %s -- su suma da 5 y tiene que dar %d.  Las dos contabilidades consistentes"
      % ("PASS" if ok6 else "FAIL", D - 2))
    P("        dan %d cada una; la suya toma la mitad de una y la mitad de la otra." % (D - 2))
    P("        Un vector 4D tiene 2 grados sin masa y 3 con masa.  NUNCA cuatro.")
    if not ok6:
        fails.append("C6")

    # ------------------------------------------------------------------ C7
    P("\n[C7] HABA-YAMASHITA ESCRIBEN -(D-2) EXPLICITO EN LA FORMULA GENERAL 5D")
    h = find(hy_pages, r"D\s*[−–-]\s*2|\(D-2\)|ghost")
    quote(h, 6)
    ok7 = bool(h)
    P("     C7 %s -- %s" % ("PASS" if ok7 else "FAIL",
                            "el articulo que da la formula general SI cuenta D-2 y SI nombra fantasmas"
                            if ok7 else "no se localiza"))
    if not ok7:
        fails.append("C7")

    # ------------------------------------------------------------------ C8
    P("\n[C8] LA ESCAPATORIA 'SERA UNA CONVENCION' -- CERRADA")
    P("       Una normalizacion global multiplica los TRES coeficientes por el mismo factor.")
    P("       El cambio que discutimos mueve 4 -> 3 en el canal periodico y deja el")
    P("       antiperiodico intacto:")
    P("")
    P("         bracket con 4 dof periodicos : (2, 4, 7)")
    P("         bracket con 3 dof periodicos : (3/2, 3, 6)")
    r = [2 / 1.5, 4 / 3, 7 / 6]
    P("         cocientes                    : %.4f  %.4f  %.4f" % tuple(r))
    ok8 = not (abs(r[0] - r[2]) < 1e-9)
    P("\n     C8 %s -- los cocientes NO son iguales (%.4f contra %.4f), luego ninguna"
      % ("PASS" if ok8 else "FAIL", r[0], r[2]))
    P("        normalizacion global puede producir esta diferencia.  Si fuera convencion,")
    P("        los tres se moverian juntos.  No se mueven juntos.")
    if not ok8:
        fails.append("C8")

    # ------------------------------------------------------------------ C9
    P("\n[C9] EL ARTICULO PREVIO DEL PROPIO MARU, EN 6D Y EN EL MISMO ORBIFOLD")
    ah_pages = pages_of(AH)
    ah_head = ah_pages[0][1][:300].replace("\n", " ")
    P("     %s" % ah_head[:190])
    P("")
    h = find(ah_pages, r"number of degrees of freedom of the fields")
    quote(h, 2)
    h = find(ah_pages, r"one-loop effective potential for an SU\(4\) gauge field")
    quote(h, 2)
    gh_ah = sum(len(re.findall(w, t, re.I)) for w in ("ghost", "Faddeev", "Popov")
                for _, t in ah_pages)
    secs = [s for s in ("SU(4) gauge fields", "Fermions") if find(ah_pages, re.escape(s))]
    P("       su ec. (3.1)  :  V = ((-1)^F / 2) * N * Sum log(p^2 + M^2)")
    P("       su ec. (3.11) :  V_gauge = (4 / 2) * Sum {...}      -> N = 4 para TODO el sector")
    P("       secciones del capitulo 3 : %s   (no hay seccion de A_5,6 aparte)" % ", ".join(secs))
    P("       menciones de fantasma en ese articulo : %d" % gh_ah)
    ok9 = gh_ah == 0
    P("\n     C9 %s -- el mismo autor, en 6D y en T^2/Z_2, pone N = 4 para el sector gauge"
      % ("PASS" if ok9 else "FAIL"))
    P("        ENTERO, que es D - 2.  En KM25 el mismo sector pesa 4 + 1 = 5.")
    P("        HONESTIDAD: esto NO prueba por si solo que 4 sea lo correcto.  Cabe leerlo al reves")
    P("        -- que AHMN se dejara A_5,6 y KM25 lo anadiera.  Lo que si establece es que el")
    P("        coeficiente del grupo era 4, y que al pasar a 5 se anadio el escalar sin restar")
    P("        nunca los fantasmas.  El ancla dura sigue siendo la formula general de HY.")
    if not ok9:
        fails.append("C9")

    # ------------------------------------------------------------------ C10
    P("\n[C10] Y LA FORMULA GENERAL DE HY, CONTRASTADA CONTRA LA LITERATURA PREVIA")
    h = find(hy_pages, r"gauge sector contribution is")
    quote(h, 3)
    h = find(hy_pages, r"available for more than 5 dimensional")
    quote(h, 2)
    ok10 = bool(find(hy_pages, r"gauge sector contribution is"))
    P("\n     C10 %s -- en 5D su coeficiente gauge es -3 = -(D-2) y dicen que reproduce"
      % ("PASS" if ok10 else "FAIL"))
    P("         Refs.[7,13].  Y dicen que el metodo vale para mas de cinco dimensiones.")
    P("         En 6D eso da 4, no 5, y lo dicen ellos.")
    if not ok10:
        fails.append("C10")

    # ------------------------------------------------------------------ C11
    P("\n[C11] HOSOTANI-NODA-TAKENAGA: LA MISMA CUENTA, EN 6D Y EN EL MISMO ORBIFOLD T^2/Z_2")
    hn_pages = pages_of(HN)
    P("     %s" % hn_pages[0][1][:170].replace("\n", " "))
    P("")
    for pat in (r"Gauge fields and ghosts", r"mass matrix for ghost",
                r"4 −2 = 2 times|4 -2 = 2 times|= 2 times contributions",
                r"two extra-dimensional components"):
        quote(find(hn_pages, pat), 2)
    gh_hn = sum(len(re.findall(w, t, re.I)) for w in ("ghost",) for _, t in hn_pages)
    P("")
    P("       su §4.1 se titula 'Gauge fields and ghosts' ; menciones de 'ghost': %d" % gh_hn)
    P("       su ec. (4.6) : V_{gauge+ghost} = -(i/2) tr ln D_L D^L")
    P("       A_mu + fantasmas = 4 - 2 = 2 ; mas los DOS componentes extradimensionales = 4")
    ok11 = bool(find(hn_pages, r"mass matrix for ghost")) and gh_hn > 0
    P("\n     C11 %s -- en SEIS dimensiones y en el MISMO orbifold, la resta de fantasmas"
      % ("PASS" if ok11 else "FAIL"))
    P("         esta escrita explicitamente y da D - 2 = 4.  Ademas dicen que la matriz de masas")
    P("         de los fantasmas es LA MISMA que la de A_mu, que es justo la hipotesis de")
    P("         contorno que necesitaba el determinante de una sola orbita.")
    P("         AVISO: HNT no tiene la segunda paridad P_6 de KM, asi que fijan el TOTAL (4),")
    P("         no el reparto entre los dos sectores.  El reparto sigue siendo nuestro.")
    if not ok11:
        fails.append("C11")

    # ------------------------------------------------------------------ C12
    P("\n[C12] LA HIPOTESIS DEL TEOREMA 1, DICHA DE FORMA QUE NO SE PUEDA MALINTERPRETAR")
    P("     El teorema necesita que 8D_gauge sea IMPAR.  Con n_+ = 2 y n_- = 6 pares en la fila")
    P("     mixta, el coeficiente es 2a + 6b = 2a + 3*(2b), luego la paridad pide DOS cosas:")
    P("        (i)  2a PAR   -- o sea, a entero")
    P("        (ii) 2b IMPAR")
    P("     No basta con preguntar '.sale b = 1/2?'.  b puede ser correcto y perderse igual la")
    P("     cuantizacion, por a.  [[separate-the-group-factor-from-the-loop-coefficient]]")
    P("")
    ROWS = [(+1, 2, 1, 0), (+1, 1, 2, 0), (-1, 1, 2, 6)]      # (s, c, pares P6=+1, pares P6=-1)
    KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]

    def eightD_of(a, b):
        u = [F(0)] * 6
        for s, c, npl, nmi in ROWS:
            u[KEYS.index((s, c))] += -(npl * a + nmi * b) / 2
        A2, B2 = u[0] + 4 * u[1] + 9 * u[2], u[3] + 4 * u[4] + 9 * u[5]
        return 8 * A2 - 6 * B2

    P("     %-26s %8s %8s %10s %10s   %s"
      % ("(a, b)", "2a par", "2b impar", "8D_gauge", "8D mod 2", "HIPOTESIS TEOREMA 1"))
    verdicts = {}
    for a, b, tag in ((F(2), F(1, 2), "ec. (68) publicada"),
                      (F(3, 2), F(1, 2), "con resta de fantasmas")):
        k = eightD_of(a, b)
        c_i = (2 * a).denominator == 1 and (2 * a) % 2 == 0
        c_ii = (2 * b).denominator == 1 and (2 * b) % 2 == 1
        odd = k.denominator == 1 and k % 2 == 1
        verdicts[tag] = odd
        P("     %-26s %8s %8s %10s %10s   %s"
          % ("(%s, %s)" % (a, b), c_i, c_ii, k, k % 2 if k.denominator == 1 else "n/a",
             "SATISFECHA" if odd else "NO SATISFECHA"))
        assert odd == (c_i and c_ii), "la regla de paridad no coincide con el 8D calculado"
    P("")
    P("     control: con (2, 1/2) tiene que salir 8D_gauge = -27.  Sale %s."
      % eightD_of(F(2), F(1, 2)))
    ok12 = eightD_of(F(2), F(1, 2)) == -27 and eightD_of(F(3, 2), F(1, 2)) == -18
    P("     C12 %s" % ("PASS" if ok12 else "FAIL"))
    if not ok12:
        fails.append("C12")
    P("")
    P("     " + "-" * 88)
    for tag, odd in verdicts.items():
        P("     THEOREM_1_HYPOTHESIS (%-24s) : %s" % (tag, "SATISFIED" if odd else "NOT SATISFIED"))
    P("     " + "-" * 88)
    P("     El teorema NO esta en cuestion.  Su demostracion es correcta en los dos casos.")
    P("     Lo que esta en cuestion es si ESTE modelo satisface su premisa.")

    # ------------------------------------------------------------------ veredicto
    P("\n" + "=" * 98)
    P("VEREDICTO")
    P("=" * 98)
    P("""
  Lo que esta ESTABLECIDO, y sale entero de su texto:

    - definen N como grados FISICOS                                   (su ec. 49)
    - fijan el gauge de forma covariante                              (su ec. 51)
    - no nombran al fantasma ni una vez en 42 paginas                 (C4)
    - toman N = 4 para A_mu y N = 1 para A_5,6                        (C5)
    - 4 + 1 = 5, y en seis dimensiones tiene que dar 4                (C6)
    - el 1 esta contado tras quitar el Goldstone; el 4, no            (su nota al pie 3)
    - no puede arreglarlo una normalizacion global                    (C8)

  Lo que NO esta establecido, y hay que decirlo asi:

    que se hayan equivocado.  Cabe una convencion o un jacobiano que no escriben.  Pero la
    condicion que falla es NECESARIA en su propio marco, no comoda para nosotros, y la unica
    pieza capaz de explicar exactamente la diferencia es la que no aparece nombrada.

  Lo que esto le hace al Teorema 1:

    no toca su demostracion, que sigue siendo correcta.  Toca su HIPOTESIS.  El teorema pasa a
    enunciarse condicionado -- "dados los coeficientes gauge del potencial reducido tal y como
    se imprimen en [KM25]" -- y esa condicion tiene que estar en el resumen, no en la seccion 13.
""")
    P("=" * 98)
    if fails:
        P("FALLAN: %s" % ", ".join(fails))
        P("=" * 98)
        return 1
    P("TODOS LOS CONTROLES PASAN.")
    P("=" * 98)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
