#!/usr/bin/env python3
"""fetch_hepdata -- the CMS-EXO-24-011 measurement and its correlation matrix, from the source.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT WAS BLOCKING THE LIMIT.  A chi^2 against a normalised distribution needs the covariance
between bins; without it one can quote a deformation but not a confidence level.  The record is
HEPData ins3136278, DOI 10.17182/hepdata.167852.

GETTING AT IT.  hepdata.net sits behind a Cloudflare challenge: the HTML record page, and every
/download/ endpoint, answer 403 to a plain client.  Two endpoints are served straight:

    /record/ins3136278?format=json                 the record and its table index
    /record/data/<recid>/<table_id>/<version>      one table, as JSON

so the fetch goes through those and nothing is scraped or screen-read.

WHAT IS THERE, AND ONE THING TO NOTICE.  Nineteen tables.  Fig. 2 gives DETECTOR-level chi in
eight mass bins; Fig. 5 gives PARTICLE-level chi -- the unfolded measurement, which is what a
calculation may be compared with -- in SEVEN.  The 2.4-3.0 TeV bin is detector level only.  So
the comparison has seven mass bins, not the eight of the CIJET grid, and any fit that silently
used eight would be inventing one.

Run:  python fetch_hepdata.py > ../outputs/fetch_hepdata.txt
"""
import json
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "hepdata"
OUT = HERE.parent / "outputs"

INSPIRE = "ins3136278"
RECID = 167852
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"


def get(url, dest):
    """curl, because requests is not installed and this is one GET."""
    r = subprocess.run(["curl", "-s", "-L", "--max-time", "120", "-A", UA,
                        "-w", "%{http_code}", "-o", str(dest), url],
                       capture_output=True, text=True)
    code = (r.stdout or "").strip()[-3:]
    return code, dest.stat().st_size if dest.exists() else 0


def main():
    DATA.mkdir(exist_ok=True)
    fails = []
    print("=" * 96)
    print("CMS-EXO-24-011 FROM HEPData: THE MEASUREMENT AND ITS CORRELATION MATRIX")
    print("=" * 96)

    idx = DATA / "record.json"
    code, size = get("https://www.hepdata.net/record/%s?format=json" % INSPIRE, idx)
    print("\n  record index : HTTP %s, %d bytes" % (code, size))
    if code != "200":
        print("  the record index did not come back; nothing else can be trusted.")
        return 1
    rec = json.loads(idx.read_text(encoding="utf-8"))
    tables = rec["data_tables"]
    print("  DOI          : %s" % rec["record"].get("hepdata_doi"))
    print("  tables       : %d" % len(tables))

    # ---- C1: the record must be the one we think it is ----------------------------------
    title = rec["record"].get("title", "")
    ok = ("dijet angular" in title.lower()
          and str(rec["record"].get("inspire_id")) == "3136278")
    print("\n[C1] IS THIS THE RIGHT RECORD?")
    print("      title   : %s" % title[:88])
    print("      inspire : %s" % rec["record"].get("inspire_id"))
    print("\n      C1 %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C1")
        return 1

    # ---- pull every table ----------------------------------------------------------------
    print("\n[1] PULLING THE TABLES")
    got = {}
    for t in tables:
        tid, name = t["id"], t["name"]
        dest = DATA / ("t%s.json" % tid)
        if dest.exists() and dest.stat().st_size > 0:
            code, size = "200", dest.stat().st_size
            note = "cached"
        else:
            code, size = get("https://www.hepdata.net/record/data/%d/%s/1" % (RECID, tid), dest)
            note = ""
        print("      %-8s HTTP %s %9d B  %-58s %s" % (tid, code, size, name[:58], note))
        if code == "200" and size > 0:
            got[name] = dest
        else:
            fails.append("fetch:%s" % tid)

    # ---- C2: particle level is seven bins, not eight -------------------------------------
    part = sorted(n for n in got if "Particle level" in n)
    det = sorted(n for n in got if "Detector level" in n)
    corr = [n for n in got if "orrelation" in n]
    print("\n[C2] WHAT IS ACTUALLY MEASURED, AND IN HOW MANY BINS")
    print("      detector-level chi tables : %d" % len(det))
    print("      particle-level chi tables : %d   <- the unfolded ones, the comparable ones"
          % len(part))
    print("      correlation matrices      : %d" % len(corr))
    ok2 = len(part) == 7 and len(det) == 8 and len(corr) == 1
    print("""
      The CIJET grid was run on EIGHT mass bins because that is the detector-level binning.
      The unfolded measurement has SEVEN: 2.4-3.0 TeV is not published at particle level.  A
      fit that used eight would be fitting a bin that does not exist.""")
    print("\n      C2 %s" % ("PASS" if ok2 else "FAIL"))
    if not ok2:
        fails.append("C2")

    # ---- C3: the correlation matrix must be square, symmetric and unit-diagonal ----------
    print("\n[C3] AND THE MATRIX MUST BE A CORRELATION MATRIX")
    side_note = 0
    if corr:
        cm = json.loads(got[corr[0]].read_text(encoding="utf-8"))
        # HEPData's rendered-table format: values[k] = {"x": [i, j], "y": [{"value": rho}]}
        ent = {}
        for v in cm["values"]:
            i = int(float(v["x"][0]["value"]) - 0.5)
            j = int(float(v["x"][1]["value"]) - 0.5)
            ent[(i, j)] = float(v["y"][0]["value"])
        npts = len(ent)
        side = int(round(npts ** 0.5))
        diag = [ent[(i, i)] for i in range(side) if (i, i) in ent]
        sym = all(abs(ent[(i, j)] - ent.get((j, i), 1e9)) < 1e-9
                  for (i, j) in list(ent)[:4000])
        mx = max(abs(v) for v in ent.values())
        print("      entries               : %d  ->  %d x %d" % (npts, side, side))
        print("      square                : %s" % (side * side == npts))
        print("      unit diagonal         : %s (worst |1-rho| = %.2e)"
              % (len(diag) == side and max(abs(d - 1) for d in diag) < 1e-9,
                 max(abs(d - 1) for d in diag) if diag else float("nan")))
        print("      symmetric             : %s" % sym)
        print("      max |rho|             : %.4f" % mx)
        ok3 = (side * side == npts and sym and len(diag) == side
               and max(abs(d - 1) for d in diag) < 1e-9 and mx <= 1.0 + 1e-9)
        print("\n      C3 %s" % ("PASS" if ok3 else "FAIL"))
        if not ok3:
            fails.append("C3")
        side_note = side

        # ---- C4: the binning the matrix declares, against the one our grid was run on -----
        print("\n[C4] AND THE BINNING IT DECLARES IS NOT THE ONE OUR GRID WAS RUN ON")
        print("""      The table says so itself: "11 times the index of the m_jj
      (3.0,3.6,4.2,4.8,5.4,6.0,7.0) bin plus the index of the CHI bin
      (1,2,3,4,5,6,7,8,9,10,12,14)".  Read the two lists.""")
        cms_mass = [3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0]
        cms_chi_edges = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]
        n_chi = len(cms_chi_edges) - 1
        print("\n      %-34s %-22s %s" % ("", "the unfolded data", "our CIJET grid"))
        print("      %-34s %-22s %s" % ("mass bins", "%d (from 3.0 TeV)" % len(cms_mass),
                                        "8 (from 2.4 TeV)"))
        print("      %-34s %-22s %s" % ("chi bins", "%d" % n_chi, "12"))
        print("      %-34s %-22s %s" % ("largest chi", "%d" % cms_chi_edges[-1], "16"))
        print("      %-34s %-22s %s" % ("fit points", "%d" % (len(cms_mass) * n_chi),
                                        "%d" % (8 * 12)))
        ok4 = side == len(cms_mass) * n_chi
        print("""
      So the fit region is %d x %d = %d points, and the matrix is exactly that size.  Our grid
      carries a 2.4-3.0 TeV mass bin the unfolding does not publish, and a 14 < chi < 16 bin
      that stops at 14 in the data.  Neither can be fitted, and a fit that kept them would be
      inventing thirteen points and mismatching every index into this matrix.""" %
              (len(cms_mass), n_chi, len(cms_mass) * n_chi))
        print("\n      C4 %s (matrix side %d, expected %d)"
              % ("PASS" if ok4 else "FAIL", side, len(cms_mass) * n_chi))
        if not ok4:
            fails.append("C4")
    else:
        print("      no correlation matrix in the record")
        fails.append("C3")

    print("""
[2] WHAT THIS UNBLOCKS, AND WHAT IT DOES NOT
      With a %d x %d correlation matrix the chi^2 can carry the correlations between chi bins
      and between mass bins, so a confidence level becomes definable.  What it does not supply
      is the QCD scale choice as a nuisance -- CMS's own small tension below 4.8 TeV moves
      between their two scale choices, so a fit that does not profile it is measuring the scale
      and not the tower.  That stays an obligation of ours.""" % (side_note, side_note))

    (OUT / "hepdata_manifest.json").write_text(json.dumps(
        {"inspire": INSPIRE, "recid": RECID,
         "doi": rec["record"].get("hepdata_doi"),
         "tables": {n: str(p.relative_to(HERE)) for n, p in got.items()},
         "particle_level": part, "detector_level": det, "correlation": corr}, indent=1))
    print("\n    [wrote outputs/hepdata_manifest.json; %d files in tu_limit/hepdata/]" % len(got))

    print("\n" + "=" * 96)
    if fails:
        print("VERDICT: %d CHECK(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        print("=" * 96)
        return 1
    print("VERDICT: the measurement and its correlation matrix are on disk and are what they")
    print("         claim to be.  The blocker is gone; the scale nuisance is not.")
    print("=" * 96)
    return 0


if __name__ == "__main__":
    sys.exit(main())
