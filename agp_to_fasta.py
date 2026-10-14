#!/usr/bin/env python3
"""Build sequences from an AGP file and its component FASTA (e.g. a final pseudohaplotype FASTA
from hapN_final.agp and the Hifiasm contig FASTAs).

Every AGP object is assembled from its own rows, in order: 'W' rows copy the component interval
(reverse-complemented on the '-' strand); 'N' and 'U' rows insert that many N. Each slice is taken
for its own object, so a contig split between several objects is handled correctly. Checks that
every object has the length its last row declares. Component FASTAs may be gzip-compressed.

usage: agp_to_fasta.py AGP OUT.fa COMPONENTS.fa [COMPONENTS2.fa ...] [--width 60]
"""
import sys, gzip, collections

COMP = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv", "TGCANtgcanYRSWMKVHDByrswmkvhdb")

def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, parts = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs

def main(args):
    width = 60
    if "--width" in args:
        i = args.index("--width"); width = int(args[i + 1]); del args[i:i + 2]
    if len(args) < 3:
        sys.exit(__doc__)
    agp, out, comp_files = args[0], args[1], args[2:]
    comps = {}
    for f in comp_files:
        comps.update(read_fasta(f))
    rows = collections.OrderedDict()
    for line in open(agp):
        if line.startswith("#") or not line.strip():
            continue
        r = line.rstrip("\n").split("\t")
        rows.setdefault(r[0], []).append(r)
    with open(out, "w") as o:
        for obj, rs in rows.items():
            pieces = []
            for r in rs:
                if r[4] == "W":
                    s, e = int(r[6]), int(r[7])
                    piece = comps[r[5]][s - 1:e]
                    if len(piece) != e - s + 1:
                        sys.exit(f"ERROR: {r[5]} shorter than {e}")
                    pieces.append(piece.translate(COMP)[::-1] if r[8] == "-" else piece)
                elif r[4] in ("N", "U"):
                    pieces.append("N" * int(r[5]))
            seq = "".join(pieces)
            if len(seq) != int(rs[-1][2]):
                sys.exit(f"ERROR: {obj} has length {len(seq)}, AGP says {rs[-1][2]}")
            o.write(f">{obj}\n")
            for i in range(0, len(seq), width):
                o.write(seq[i:i + width] + "\n")
    print(f"objects written: {len(rows)}")

if __name__ == "__main__":
    main(sys.argv[1:])
