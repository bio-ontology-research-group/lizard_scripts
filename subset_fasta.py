#!/usr/bin/env python3
"""Copy the FASTA records whose IDs are in a keep-list, verbatim and in source order.

Records are copied line by line exactly as they appear in the input (header and
sequence lines unchanged), so the output differs from the input only by the
records left out. The ID is the header text up to the first whitespace.
Exits with an error if any keep-list ID is missing from the input or duplicated in it.

usage: subset_fasta.py keep_list.txt in.fasta out.fasta
"""
import sys

def main(keepf, inf, outf):
    keep = set(l.strip() for l in open(keepf) if l.strip())
    seen = {}
    copying = False
    with open(inf) as fin, open(outf, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                rid = line[1:].split()[0] if line[1:].split() else ""
                copying = rid in keep
                if copying:
                    seen[rid] = seen.get(rid, 0) + 1
            if copying:
                fout.write(line)
    missing = keep - set(seen)
    dup = [k for k, v in seen.items() if v > 1]
    print(f"keep-list IDs\t{len(keep)}\nrecords written\t{sum(seen.values())}\nmissing from input\t{len(missing)}\nduplicated in input\t{len(dup)}")
    if missing or dup:
        sys.exit(f"ERROR: missing {sorted(missing)[:5]} duplicated {dup[:5]}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
