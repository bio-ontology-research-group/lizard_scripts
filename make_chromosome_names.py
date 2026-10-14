#!/usr/bin/env python3
"""Build chromosome_names.tsv: for each of the 32 chromosomes, its FASTA name, pseudohaplotype,
article name (Chr1-Chr18), GenBank chromosome name, Juicebox scaffold name and length.

The GenBank chromosome names are those given to table2asn in the source tables
(source_hap1.tsv, source_hap2.tsv), which follow the article's names; the article name is
'Chr' + that name. Rows are sorted by pseudohaplotype and chromosome number.

usage: make_chromosome_names.py seqid_map.tsv source_hap1.tsv source_hap2.tsv lengths.tsv out.tsv
  seqid_map.tsv: Juicebox scaffold name -> FASTA name for the chromosomes
  lengths.tsv:   name<TAB>length[...] of the final assembly (fingerprint tables or .fai)
"""
import sys

def main(mapf, src1, src2, lenf, out):
    jb = {}
    for l in open(mapf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and not a[0].startswith(("#", "braker_seqid")):
            jb[a[1]] = a[0]
    lens = {}
    for l in open(lenf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and a[1].isdigit():
            lens[a[0]] = int(a[1])
    ncbi = {}
    for f in (src1, src2):
        for l in open(f):
            a = l.rstrip("\n").split("\t")
            if a[0] != "Sequence_ID":
                ncbi[a[0]] = a[1]
    names = sorted(ncbi, key=lambda x: (x[:4], int(x.split("chr")[1])))
    with open(out, "w") as o:
        o.write("fasta_name\tpseudohaplotype\tmanuscript_chromosome\tncbi_chromosome_name\tjuicebox_scaffold\tlength_bp\n")
        for fa in names:
            o.write("\t".join([fa, fa[3], "Chr" + ncbi[fa], ncbi[fa], jb[fa], str(lens[fa])]) + "\n")
    print(f"chromosomes: {len(names)}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
