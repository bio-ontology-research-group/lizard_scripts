#!/usr/bin/env python3
"""Extract the soft-masked (lowercase) intervals of a FASTA as BED, with final sequence names.

For each sequence: renames it with a two-column map (e.g. scaffold_N -> hapN.chrM), writes every run
of lowercase bases as a BED interval (0-based start, end exclusive), and computes the MD5 of the
uppercased sequence, so the masked FASTA can be compared with the final assembly ignoring case.
Writes <prefix>.bed (sorted by final name and position) and <prefix>.check.tsv
(final name, original name, length, masked bp, MD5 of uppercase sequence), and prints the total
masked bp.

To regenerate the masked FASTA from a final FASTA:
    bedtools maskfasta -soft -fi haplotypeN.fa -bed <prefix>.bed -fo masked.fa

usage: softmask_to_bed.py seqid_map.tsv masked.fa OUT_PREFIX
"""
import sys, re, hashlib

def fasta(path):
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)

def main(mapf, fa, prefix):
    nm = {}
    for l in open(mapf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and not a[0].startswith(("#", "braker_seqid")):
            nm[a[0]] = a[1]
    rows, checks, total = [], [], 0
    low = re.compile(r"[a-z]+")
    for name, seq in fasta(fa):
        if name not in nm:
            sys.exit(f"ERROR: {name} not in the map")
        q = nm[name]; n = 0
        for m in low.finditer(seq):
            rows.append((q, m.start(), m.end())); n += m.end() - m.start()
        total += n
        checks.append((q, name, len(seq), n, hashlib.md5(seq.upper().encode()).hexdigest()))
    key = lambda q: (q.split(".")[0], int(re.search(r"\d+$", q).group()))
    rows.sort(key=lambda r: (key(r[0]), r[1]))
    with open(f"{prefix}.bed", "w") as o:
        for q, s, e in rows:
            o.write(f"{q}\t{s}\t{e}\n")
    with open(f"{prefix}.check.tsv", "w") as o:
        o.write("final_name\toriginal_name\tlength\tmasked_bp\tmd5_uppercase\n")
        for c in sorted(checks, key=lambda c: key(c[0])):
            o.write("\t".join(map(str, c)) + "\n")
    print(f"sequences: {len(checks)}; intervals: {len(rows)}; masked bp: {total}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
