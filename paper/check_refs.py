#!/usr/bin/env python3
"""check_refs.py - every \\ref has a \\label, every \\cite has a \\bibitem, and nothing is orphaned.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Ported from part_vi/paper/check_refs.py.  Part VI ran this over one edition; here BOTH
editions are checked, because the Spanish one is a separate file that can drift on its own
and did in Part VI (see part_vi/GATE notes on the ES longtable row count).

Run:  python check_refs.py     (from part_vii/paper/)
"""
import re

TEXS = ["su7_hierarchy.tex", "su7_hierarchy_es.tex"]

bad = []
for TEX in TEXS:
    s = open(TEX, encoding="utf-8").read()
    s = re.sub(r"(?<!\\)%[^\n]*", "", s)

    labels = set(re.findall(r"\\label\{([^}]*)\}", s))
    refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", s))
    bib = set(re.findall(r"\\bibitem\{([^}]*)\}", s))
    cited = set()
    for g in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", s):
        cited |= {x.strip() for x in g.split(",")}

    print("%s: labels %d, refs %d, bibitems %d, keys cited %d"
          % (TEX, len(labels), len(refs), len(bib), len(cited)))
    for name, s_ in (("refs with no label", refs - labels),
                     ("labels never referenced", labels - refs),
                     ("cites with no bibitem", cited - bib),
                     ("bibitems never cited", bib - cited)):
        print("  %-24s : %s" % (name, ", ".join(sorted(s_)) if s_ else "none"))
        if s_ and name != "labels never referenced":
            bad += ["%s: %s: %s" % (TEX, name, x) for x in sorted(s_)]

print()
print("both editions clean" if not bad else "FLAGGED:\n  " + "\n  ".join(bad))
raise SystemExit(1 if bad else 0)
