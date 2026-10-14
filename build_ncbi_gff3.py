#!/usr/bin/env python3
"""Convert one haplotype's final GFF3 (N2a) into a GFF3 for NCBI table2asn submission.

Changes, and nothing else:
  * keeps only gene, mRNA, exon and CDS rows; drops intron, start_codon and stop_codon
    (table2asn ignores introns and does not need start/stop codon features);
  * gene:  adds locus_tag=<PREFIX>_NNNNNN, numbered in genomic order in steps of 10;
  * mRNA:  adds transcript_id=gnl|<PREFIX>|mrna.<BRAKER transcript ID>;
  * CDS:   adds protein_id=gnl|<PREFIX>|<BRAKER transcript ID> on every row of the CDS,
           and product=<name> if a product table is given (otherwise table2asn names it
           'hypothetical protein');
  * existing ID/Parent attributes, coordinates, strands and phases are unchanged.
Partialness is computed by table2asn from the CDS ends, so no partial markup is added.

usage: build_ncbi_gff3.py in.gff3 PREFIX out.gff3 out_locus_map.tsv [products.tsv]
  products.tsv: transcript_id<TAB>product name (optional; GFF3-escaped on output)
"""
import sys, collections

KEEP = {"gene", "mRNA", "exon", "CDS"}

def esc(v):
    return v.replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace("&", "%26").replace(",", "%2C")

def attrs(col9):
    return [kv.split("=", 1) for kv in col9.strip().strip(";").split(";") if "=" in kv]

def main(inp, prefix, out, mapf, prodf=None):
    prod = {}
    if prodf:
        for l in open(prodf):
            a = l.rstrip("\n").split("\t")
            if len(a) >= 2 and a[1].strip():
                prod[a[0]] = a[1].strip()
    lines = [l.rstrip("\n") for l in open(inp)]
    header = [l for l in lines if l.startswith("##")]
    feats = [l.split("\t") for l in lines if l and not l.startswith("#")]
    genes = [f for f in feats if f[2] == "gene"]           # already in genomic order (S2)
    tag = {}
    for i, g in enumerate(genes, 1):
        tag[dict(attrs(g[8]))["ID"]] = f"{prefix}_{i*10:06d}"
    mrna_gene = {}
    stats = collections.Counter()
    with open(out, "w") as o, open(mapf, "w") as m:
        o.write("\n".join(header) + "\n")
        m.write("locus_tag\tbraker_gene\tbraker_transcript\ttranscript_id\tprotein_id\tproduct\n")
        for f in feats:
            if f[2] not in KEEP:
                stats["dropped_" + f[2]] += 1
                continue
            a = attrs(f[8]); d = dict(a)
            if f[2] == "gene":
                a.append(["locus_tag", tag[d["ID"]]])
            elif f[2] == "mRNA":
                mrna_gene[d["ID"]] = d["Parent"]
                a.append(["transcript_id", f"gnl|{prefix}|mrna.{d['ID']}"])
                p = prod.get(d["ID"], "")
                m.write(f"{tag[d['Parent']]}\t{d['Parent']}\t{d['ID']}\tgnl|{prefix}|mrna.{d['ID']}\tgnl|{prefix}|{d['ID']}\t{p or 'hypothetical protein'}\n")
            elif f[2] == "CDS":
                a.append(["protein_id", f"gnl|{prefix}|{d['Parent']}"])
                if d["Parent"] in prod:
                    a.append(["product", esc(prod[d["Parent"]])])
            f[8] = ";".join(f"{k}={v}" for k, v in a) + ";"
            stats[f[2]] += 1
            o.write("\t".join(f) + "\n")
    for k in sorted(stats):
        print(f"{k}\t{stats[k]}")
    print(f"locus_tags\t{len(tag)}\tfirst {min(tag.values())}\tlast {max(tag.values())}")

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
