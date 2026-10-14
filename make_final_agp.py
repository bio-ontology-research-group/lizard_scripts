#!/usr/bin/env python3
"""Build hap1_final.agp and hap2_final.agp from the Juicebox AGP.

1. Translate each Juicebox object to its final sequence name (jbat_to_final.tsv; '-' = dropped,
   the four mitochondrial contigs) and write each pseudohaplotype's objects in order: chromosomes
   by number, then unplaced scaffolds by number. Component rows are copied unchanged.
2. Add the gaps that Hifiasm (--dual-scaf) had already written inside its output sequences: every
   N run in the final FASTAs (seqtk cutN table) that is not an AGP gap lies inside one component,
   which is split around it; the gap is 'U' if 100 bp long, otherwise 'N' with its length,
   type scaffold, linkage yes, evidence align_genus. Part numbers are renumbered.

usage: make_final_agp.py out_JBAT.FINAL.agp jbat_to_final.tsv nruns.tsv OUT_DIR
  nruns.tsv: `seqtk cutN -gn 1` output of the final FASTAs (name, 0-based start, end), both haplotypes
"""
import re, sys, collections, os

HEADER = ["##agp-version\t2.1",
          "# Scincus mitranus isolate OM-2025, pseudohaplotype {n}: final sequences (A7) described in terms of Hifiasm contigs",
          "# Derived from out_JBAT_FINAL.agp (Juicebox) with object names translated to final sequence names; removed objects dropped",
          "# Includes the 4 gaps inserted by Hifiasm 0.25.0 --dual-scaf inside its output sequences (evidence align_genus: joined by alignment to the homologous haplotype)"]

def key(n):
    return (n.split(".")[0], 0 if ".chr" in n else 1, int(re.search(r"\d+$", n).group()))

def main(jbat, mapf, nrunf, outdir):
    final = {}
    for l in open(mapf):
        a = l.rstrip("\n").split("\t")
        if a[0] != "jbat_object" and a[1] != "-":
            final[a[0]] = a[1]
    J = collections.OrderedDict()
    for l in open(jbat):
        if l.startswith("#") or not l.strip():
            continue
        r = l.rstrip("\n").split("\t"); J.setdefault(r[0], []).append(r)
    rows = collections.defaultdict(list)            # final name -> AGP rows (final name in column 1)
    for o, rs in J.items():
        if o in final:
            rows[final[o]] = [[final[o]] + r[1:] for r in rs]
    agp_gaps = {(r[0], int(r[1]), int(r[2])) for v in rows.values() for r in v if r[4] in ("N", "U")}
    for l in open(nrunf):
        if not l.strip():
            continue
        n, s0, e = l.split()[:3]; s, e = int(s0) + 1, int(e)
        if (n, s, e) in agp_gaps:
            continue
        new = []
        for r in rows[n]:
            if r[4] == "W" and int(r[1]) <= s and int(r[2]) >= e:
                os_, oe, ctg, cs, ce, st = int(r[1]), int(r[2]), r[5], int(r[6]), int(r[7]), r[8]
                if st == "+":
                    p1 = (os_, s - 1, cs, cs + (s - 1 - os_)); p2 = (e + 1, oe, cs + (e + 1 - os_), ce)
                else:
                    p1 = (os_, s - 1, cs + (oe - (s - 1)), ce); p2 = (e + 1, oe, cs, cs + (oe - (e + 1)))
                L = e - s + 1
                new += [[n, str(p1[0]), str(p1[1]), "", "W", ctg, str(p1[2]), str(p1[3]), st],
                        [n, str(s), str(e), "", "U" if L == 100 else "N", str(L), "scaffold", "yes", "align_genus"],
                        [n, str(p2[0]), str(p2[1]), "", "W", ctg, str(p2[2]), str(p2[3]), st]]
            else:
                new.append(r)
        rows[n] = new
    for hp in ("hap1", "hap2"):
        lines = [h.format(n=hp[-1]) for h in HEADER]
        for name in sorted((x for x in rows if x.startswith(hp + ".")), key=key):
            for k, r in enumerate(rows[name], 1):
                r[3] = str(k); lines.append("\t".join(r))
        with open(os.path.join(outdir, f"{hp}_final.agp"), "w") as o:
            o.write("\n".join(lines) + "\n")
        print(f"{hp}_final.agp: {sum(1 for x in rows if x.startswith(hp + '.'))} objects")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
