#!/usr/bin/env python3
"""Build products.tsv: the GenBank product name of every final transcript.

Product name = the recommended name ('RecName: Full=') of the transcript's best SwissProt hit,
edited to GenBank conventions in three passes:
  1. take RecName Full (drop AltName/Short/Flags and the species); strip a leading
     'PUTATIVE PSEUDOGENE:'; drop {evidence} tags; 'Uncharacterized protein ...' -> 'hypothetical protein';
  2. apply the corrections table2asn made automatically, plus one for bracketed endings:
     'Probable ' / 'Novel ' -> 'putative '; 'Uncharacterized X' (X not 'protein') -> 'putative X';
     'ribonucleoprotein protein' -> 'ribonucleoprotein'; 'sulph' -> 'sulf'; trailing '[X]' -> '(X)';
  3. names containing 'open reading frame' (rejected by GenBank) -> 'hypothetical protein'.
Columns: transcript_id, product (final), swissprot_best_hit, identity_pct,
swissprot_recname (the name after pass 1). Rows sorted by transcript ID.

usage: make_products.py diamond_best_hits.tsv swissprot_titles.tsv keep_list.txt products.tsv
  diamond_best_hits.tsv: DIAMOND tabular output, one best hit per query (qseqid sseqid pident ...)
  swissprot_titles.tsv: accession and title from `blastdbcmd -outfmt "%a\\t%t"`; the separator may
                        be a tab or the literal two characters backslash-t
"""
import re, sys

def pass1(t):
    t = re.sub(r"^PUTATIVE PSEUDOGENE:\s*", "", t)
    m = re.search(r"RecName: Full=(.*?)(?:;| \[[^\[\]]*\]$|$)", t)
    n = m.group(1).strip() if m else t
    n = re.sub(r"\s*\{[^}]*\}", "", n).strip()
    if re.match(r"(?i)uncharacterized protein", n):
        return "hypothetical protein"
    return n

RULES = [(r"^Probable ", r"putative "),
         (r"^Novel ", r"putative "),
         (r"^Uncharacterized (?!protein\b)", r"putative "),
         (r"\bribonucleoprotein protein\b", r"ribonucleoprotein"),
         (r"[Ss]ulph", lambda m: "sulf" if m.group(0)[0] == "s" else "Sulf"),
         (r"\s*\[([^\[\]]+)\]$", r" (\1)")]

def main(hitsf, titlesf, keepf, out):
    hit = {}
    for r in open(hitsf):
        if r.strip():
            f = r.split("\t"); hit[f[0]] = (f[1], f[2])
    title = {}
    for l in open(titlesf):
        if l.strip():
            l = l.rstrip("\n")
            a, t = l.split("\\t", 1) if "\\t" in l else l.split("\t", 1)
            title[a] = t
    keep = [l.strip() for l in open(keepf) if l.strip()]
    rows = []
    for t in sorted(keep):
        acc, pid = hit[t]
        orig = pass1(title[acc])
        q = orig
        for pat, rep in RULES:
            q = re.sub(pat, rep, q)
        if "open reading frame" in q.lower():
            q = "hypothetical protein"
        rows.append("\t".join([t, q, acc, pid, orig]))
    with open(out, "w") as o:
        o.write("transcript_id\tproduct\tswissprot_best_hit\tidentity_pct\tswissprot_recname\n")
        o.write("\n".join(rows) + "\n")
    print(f"transcripts {len(rows)}; hypothetical protein {sum(r.split(chr(9))[1] == 'hypothetical protein' for r in rows)}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
