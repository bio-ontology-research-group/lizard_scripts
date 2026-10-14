#!/usr/bin/env python3
"""Build the final (SwissProt-filtered) GFF3 for Scincus mitranus from BRAKER3's braker.gff3.

Steps, in order:
  1. keep mRNAs listed in the keep-list, all their child features, and the parent genes
     of kept mRNAs (a gene is kept if at least one of its mRNAs is kept);
  2. recompute each kept gene's start/end as the span of its kept mRNAs;
  3. set phase to '.' on every feature that is not a CDS (GFF3 allows phase only on CDS);
  4. rename sequence IDs with a two-column map (BRAKER name -> assembly name);
  5. write '##gff-version 3' and '##sequence-region' lines, then genes sorted by
     sequence (assembly order) and start, each followed by its features in BRAKER's order.

No other attribute, coordinate or feature is changed.

usage: build_final_gff3.py braker.gff3 seqid_map.tsv lengths.tsv out.gff3 [keep_list.txt]
  lengths.tsv: name<TAB>length[<TAB>...] for the final assembly sequences (e.g. the
               fingerprint tables or a .fai index)
  keep_list.txt: transcript IDs to keep; omit it to convert the complete BRAKER3 set.
"""
import sys, re, collections

def attrs(col9):
    return dict(kv.split("=", 1) for kv in col9.strip().strip(";").split(";") if "=" in kv)

def main(gff, mapf, lenf, out, keepf=None):
    seqmap = {}
    for l in open(mapf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and not a[0].startswith(("#", "braker_seqid")):
            seqmap[a[0]] = a[1]
    lengths = {}
    for l in open(lenf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and a[1].isdigit():
            lengths[a[0]] = int(a[1])
    keep = set(l.strip() for l in open(keepf) if l.strip()) if keepf else None

    genes, gene_order = {}, []
    mrna_of_gene = collections.defaultdict(list)
    children = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            sys.exit(f"ERROR: line without 9 columns: {line[:80]}")
        a = attrs(f[8])
        if f[2] == "gene":
            genes[a["ID"]] = f; gene_order.append(a["ID"])
        elif f[2] == "mRNA":
            mrna_of_gene[a["Parent"]].append(f)
        else:
            children[a["Parent"]].append(f)

    stats = collections.Counter()
    out_genes = []
    for gid in gene_order:
        g = list(genes[gid])
        kept = [m for m in mrna_of_gene[gid] if keep is None or attrs(m[8])["ID"] in keep]
        if not kept:
            continue
        s, e = min(int(m[3]) for m in kept), max(int(m[4]) for m in kept)
        if (s, e) != (int(g[3]), int(g[4])):
            stats["gene_span_recomputed"] += 1
            g[3], g[4] = str(s), str(e)
        block = [g]
        for m in kept:
            block.append(list(m))
            block.extend(list(c) for c in children[attrs(m[8])["ID"]])
        for f in block:
            if f[2] != "CDS" and f[7] != ".":
                stats["phase_cleared_" + f[2]] += 1
                f[7] = "."
            if f[0] not in seqmap:
                sys.exit(f"ERROR: sequence {f[0]} not in the map")
            f[0] = seqmap[f[0]]
        stats["genes"] += 1; stats["mRNAs"] += len(kept)
        out_genes.append(block)

    def seq_key(n):
        hap, rest = n.split(".", 1)
        return (hap, 0 if rest.startswith("chr") else 1, int(re.search(r"\d+$", rest).group()))
    out_genes.sort(key=lambda b: (seq_key(b[0][0]), int(b[0][3]), int(b[0][4])))
    used = sorted({b[0][0] for b in out_genes}, key=seq_key)
    with open(out, "w") as o:
        o.write("##gff-version 3\n")
        for n in used:
            o.write(f"##sequence-region {n} 1 {lengths[n]}\n")
        for b in out_genes:
            for f in b:
                o.write("\t".join(f) + "\n")
    for k in sorted(stats):
        print(f"{k}\t{stats[k]}")

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
