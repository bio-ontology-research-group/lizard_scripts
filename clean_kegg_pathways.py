#!/usr/bin/env python3
"""Rewrite kegg_parser's KEGG pathway CSV as a clean tab-separated table.

One output row per input row (one row per pathway map and KO). Fixes:
  * rows where the pathway name was split over two fields by an unquoted
    'Including: ...' description (maps 00998, 00999) are rejoined;
  * the pathway field is split into map ID, name and the KO count KEGG Mapper
    appends in parentheses;
  * transcript lists are written comma-separated without spaces;
  * output is tab-separated with LF line endings.
The KO field is located by its K##### pattern, so a split name can't shift columns.

usage: clean_kegg_pathways.py in.csv out.tsv
"""
import csv, re, sys
inp, out = sys.argv[1], sys.argv[2]
rows = list(csv.reader(open(inp, newline="", encoding="utf-8")))
hdr, data = rows[0], rows[1:]
with open(out, "w", newline="") as o:
    o.write("super_category\tcategory\tmap_id\tpathway_name\tn_KO_in_map\tKO\ttranscripts\n")
    for r in data:
        ki = [i for i, v in enumerate(r) if re.fullmatch(r"K\d{5}", v.strip())]
        if len(ki) != 1 or ki[0] < 3:
            sys.exit(f"ERROR: cannot locate KO in row {r}")
        k = ki[0]
        pw = ", ".join(x.strip() for x in r[2:k])
        m = re.fullmatch(r"(\d{5}) (.*?)\s*\((\d+)\)", pw)
        if not m:
            sys.exit(f"ERROR: unexpected pathway field {pw!r}")
        tx = ",".join(t.strip() for t in r[k + 1].split(",") if t.strip()) if len(r) > k + 1 else ""
        for v in (r[0], r[1], m.group(2), tx):
            if "\t" in v or "\n" in v:
                sys.exit("ERROR: tab or newline inside a field")
        o.write("\t".join([r[0].strip(), r[1].strip(), "map" + m.group(1), m.group(2), m.group(3), r[k].strip(), tx]) + "\n")
print(f"rows written: {len(data)}")
