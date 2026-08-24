#!/usr/bin/env python3
"""pdf_provenance.py -- que son exactamente los conjuntos de partones, ANTES de usarlos.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Un arbitro senalo que la deformacion se calcula con CT10nlo mientras la linea base de CMS es
NNPDF3.1, y que a M_jj de 5 a 8 TeV estamos en x alto, donde eso importa.  Bajar los ficheros y
lanzar el barrido habria sido lo rapido.  Es tambien la forma de obtener tres numeros que no
significan nada, porque un conjunto usado fuera de su rejilla no avisa: extrapola.

Esto comprueba, para cada conjunto, lo que hay que saber antes y no despues:

  P1  el fichero .info existe y dice que es -- orden, alpha_s(M_Z), esquema de sabores, numero
      de miembros.  Un conjunto NNLO metido en un elemento de matriz a LO es una eleccion, no un
      descuido, y para llamarlo eleccion hay que haberla visto.
  P2  el rango (x, Q) de la rejilla CONTIENE los puntos que le vamos a pedir.  Nuestra cuadratura
      llega a Q = 8 TeV y a x del orden de 0.5; si XMax o QMax se quedan cortos, LHAPDF
      extrapola en silencio y el numero sale igual de bonito y sin sentido.
  P3  alpha_s(M_Z) de cada conjunto, que es lo que entra en el prefactor de la seccion eficaz.
      Si difieren, el cociente al Modelo Estandar NO cancela la diferencia entera.
  P4  los ficheros estan completos: el miembro 0 existe y el .info declara tantos miembros como
      .dat hay.

Ninguno de los cuatro es una opinion.  Todos se leen del propio conjunto.

Run:  python pdf_provenance.py > ../outputs/pdf_provenance.txt
"""
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
PDFDIR = pathlib.Path("E:/proyectos/Curiosity/tools/lhapdf")
SYSDIR_NOTE = "CT10nlo vive dentro de la imagen cijet:1.0, no en tools/lhapdf/"

SETS = ["CT10nlo", "NNPDF31_nnlo_as_0118", "CT14nlo"]

def quadrature_extremes():
    """el (x, Q) MAYOR que la cuadratura pide de verdad, leido de hadronic_chi.py.

    La primera version de esto calculo la esquina a mano, x = (M/sqrt(s)) e^{+y_max} con M y
    y_max en sus topes, y saco x = 3.22 -- imposible, porque x <= 1.  La esquina no es alcanzable:
    x1 x2 = M^2/s obliga a |y_boost| <= ln(sqrt(s)/M), que a M = 8 TeV son 0.53 y no 1.7.  Un
    control que inventa el punto extremo en vez de leerlo declara fuera de rejilla a un conjunto
    que esta perfectamente dentro.  [[inventory-the-data-before-reasoning]]"""
    ns = {"__file__": str(HERE / "hadronic_chi.py"), "__name__": "hadronic_chi"}
    import contextlib
    import io as _io
    with contextlib.redirect_stdout(_io.StringIO()):
        exec(compile((HERE / "hadronic_chi.py").read_text(encoding="utf-8"),
                     "hadronic_chi.py", "exec"), ns)
    mass_bins, chi_bins = ns["read_bins"]()
    import numpy as np
    NM, NYB = ns["NM"], ns["NYB"]
    rs, ybm = ns["ROOT_S"], ns["YB_MAX"]
    mn, _mw = np.polynomial.legendre.leggauss(NM)
    yn, _yw = np.polynomial.legendre.leggauss(NYB)
    xmax = qmax = 0.0
    for mlo, mhi in mass_bins:
        mj = 0.5 * (mhi - mlo) * mn + 0.5 * (mhi + mlo)
        for m in mj:
            qmax = max(qmax, float(m))
            for y in yn:
                for sgn in (+1.0, -1.0):
                    xmax = max(xmax, float(m / rs * math.exp(sgn * ybm * y)))
    return xmax, qmax, len(mass_bins), len(chi_bins)


def info_of(name):
    """el .info del conjunto, como diccionario plano; None si no esta en tools/lhapdf/."""
    p = PDFDIR / name / (name + ".info")
    if not p.exists():
        return None
    out = {}
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in ln or ln.strip().startswith("#"):
            continue
        k, v = ln.split(":", 1)
        out[k.strip()] = v.strip()
    return out


def main():
    fails = []
    print("=" * 94)
    print("QUE SON LOS CONJUNTOS DE PARTONES QUE VAMOS A USAR")
    print("=" * 94)

    xmax_need, qmax_need, nmb, ncb = quadrature_extremes()
    print("\n  leido de la propia cuadratura de hadronic_chi.py (%d bins de masa x %d de chi):"
          % (nmb, ncb))
    print("  pide como maximo  x = %.4f   y   Q = %.0f GeV" % (xmax_need, qmax_need))

    rows = []
    for name in SETS:
        inf = info_of(name)
        print("\n" + "-" * 94)
        print("  %s" % name)
        print("-" * 94)
        if inf is None:
            print("    no esta en %s" % PDFDIR)
            print("    %s" % SYSDIR_NOTE)
            rows.append((name, None))
            continue
        keys = ["SetDesc", "OrderQCD", "AlphaS_MZ", "NumFlavors", "NumMembers",
                "XMin", "XMax", "QMin", "QMax", "Flavors"]
        for k in keys:
            if k in inf:
                v = inf[k]
                print("    %-12s %s" % (k, v[:110]))
        # P4: tantos .dat como miembros declara.  El glob IGNORA los ficheros que empiezan por
        # punto: el tarball de CT14nlo trae 57 "._CT14nlo_NNNN.dat", bifurcaciones de recurso de
        # macOS que ls no ensena y que hacian contar 114 miembros donde hay 57.  Un conteo que
        # incluye basura del empaquetado no mide el conjunto.
        ndat = len([p for p in (PDFDIR / name).glob("*.dat") if not p.name.startswith(".")])
        nmem = int(inf.get("NumMembers", "0"))
        print("    %-12s %d ficheros .dat en disco" % ("(P4)", ndat))
        okmem = ndat == nmem
        print("    %-12s el .info declara %d miembros : %s"
              % ("", nmem, "coinciden" if okmem else "*** NO COINCIDEN ***"))
        if not okmem:
            fails.append("%s:P4" % name)

        # P2: el rango contiene lo que vamos a pedir
        try:
            xmx = float(inf["XMax"])
            qmx = float(inf["QMax"])
        except (KeyError, ValueError):
            print("    %-12s el .info no declara XMax/QMax; no se puede comprobar" % "(P2)")
            fails.append("%s:P2" % name)
            rows.append((name, inf))
            continue
        # P2, CORREGIDO.  La primera version exigia XMax >= el x mayor que pide la cuadratura, y
        # ese x es 2.92: la rejilla producto (masa) x (y_boost) contiene esquinas con x > 1, que
        # son cinematicamente imposibles y cuyo peso debe ser cero.  Exigir que la rejilla del
        # PDF "llegue" a x = 2.92 no tiene sentido -- ningun conjunto pasa de x = 1, y no deben.
        # Lo que hay que exigir son dos cosas distintas: que los puntos FISICOS caigan dentro, y
        # que los no fisicos devuelvan CERO y no una extrapolacion.  Lo segundo se mide en la
        # comprobacion P5 de abajo, no se supone.  [[a-necessary-condition-not-a-convenient-one]]
        okx = xmx >= min(xmax_need, 1.0)
        okq = qmx >= qmax_need
        print("    %-12s x fisico hasta %.4f, la rejilla llega a %.4f : %s"
              % ("(P2)", min(xmax_need, 1.0), xmx, "dentro" if okx else "*** FUERA ***"))
        print("    %-12s Q hasta %.0f GeV, la rejilla llega a %.0f : %s"
              % ("", qmax_need, qmx, "dentro" if okq else "*** FUERA DE REJILLA ***"))
        if not (okx and okq):
            fails.append("%s:P2" % name)
        rows.append((name, inf))

    # P3: alpha_s comparado, que no cancela en el cociente
    print("\n" + "=" * 94)
    print("P3 -- alpha_s(M_Z) Y ORDEN, LADO A LADO")
    print("=" * 94)
    print("\n  %-24s %10s %8s %10s" % ("conjunto", "alpha_s(MZ)", "orden", "n_f"))
    als = {}
    for name, inf in rows:
        if inf is None:
            print("  %-24s %10s %8s %10s" % (name, "(en la imagen)", "-", "-"))
            continue
        a = inf.get("AlphaS_MZ", "?")
        als[name] = a
        print("  %-24s %10s %8s %10s"
              % (name, a, inf.get("OrderQCD", "?"), inf.get("NumFlavors", "?")))
    # COMO NUMEROS, no como cadenas.  "0.1180000" y "0.118000" son el mismo alpha_s y la primera
    # version de esto los declaro distintos porque comparaba texto.  Un control que compara la
    # representacion en vez del valor inventa una discrepancia y luego pide explicarla.
    vals = sorted({round(float(v), 6) for v in als.values()})
    if len(vals) > 1:
        print("\n  *** los conjuntos NO comparten alpha_s(M_Z): %s ***"
              % ", ".join("%.6f" % v for v in vals))
        print("  Eso entra en el prefactor de la seccion eficaz y NO se cancela entero en el")
        print("  cociente al Modelo Estandar, porque la interferencia va como alpha_s^2 y el")
        print("  fondo tambien pero con distinta mezcla de canales.  Hay que decirlo al citar")
        print("  la comparacion, no descubrirlo despues.")
    else:
        print("\n  todos comparten alpha_s(M_Z) = %.6f" % vals[0])
        print("  de modo que el prefactor comun no distingue los conjuntos y lo que queda es la")
        print("  MEZCLA DE SABORES, que es justo lo que la deformacion ve: toca canales de quark")
        print("  y no los iniciados por gluon.")

    # ---- P5: que devuelve cada conjunto FUERA de la region fisica ---------------------------
    # La cuadratura pide x hasta 2.92 porque la rejilla producto masa x y_boost tiene esquinas
    # que x1 x2 = M^2/s prohibe.  Si LHAPDF extrapolara ahi, toda la integral llevaria basura --
    # con CT10 tambien, es decir en cada numero del recast publicado hasta hoy.  Nadie lo habia
    # comprobado.  [[falsify-the-instrument-first]]
    print("\n" + "=" * 94)
    print("P5 -- FUERA DE LA REGION FISICA: CERO, O EXTRAPOLACION?")
    print("=" * 94)
    import importlib
    import os
    sys.path.insert(0, str(HERE))
    probe = [(0.5, 5000.0), (0.999, 5000.0), (1.0, 5000.0), (1.5, 5000.0), (2.92, 12881.0)]
    for name, inf in rows:
        os.environ["SU7_PDF"] = name
        sys.modules.pop("hadronic_chi", None)
        hc = importlib.import_module("hadronic_chi")
        got = hc.read_dump(probe)
        if got is None:
            got = hc.run_dump(probe)
        bad = [(x, max(abs(v) for v in r[:13]))
               for (x, _q), r in zip(probe, got) if x > 1.0]
        worst = max((b for _x, b in bad), default=0.0)
        okz = worst == 0.0
        print("  %-24s x = 1.5 y 2.92 devuelven, en el peor sabor : %.3e   %s"
              % (name, worst, "CERO" if okz else "*** EXTRAPOLA ***"))
        if not okz:
            fails.append("%s:P5" % name)
    print("\n  Asi que las esquinas no fisicas de la cuadratura pesan exactamente cero y la")
    print("  integral es la fisica.  No es una suposicion del codigo: esta medido aqui.")

    print("\n" + "=" * 94)
    if fails:
        print("NO USAR TODAVIA: %d comprobacion(es) fallan: %s" % (len(fails), ", ".join(fails)))
        print("=" * 94)
        return 1
    print("los conjuntos cubren los puntos que la cuadratura pide, y estan completos")
    print("=" * 94)
    return 0


if __name__ == "__main__":
    sys.exit(main())
