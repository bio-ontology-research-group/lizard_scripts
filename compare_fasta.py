#!/usr/bin/env python3
"""Compare two FASTA files by sequence, ignoring line wrapping, header text after the ID,
letter case and a terminal stop symbol ('*' or '.').
usage: compare_fasta.py a.fasta b.fasta
"""
import sys
def read(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l.startswith(">"):
            if n is not None: d[n] = "".join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n is not None: d[n] = "".join(b)
    return {k: v.upper().rstrip("*.") for k, v in d.items()}
a, b = read(sys.argv[1]), read(sys.argv[2])
diff = [k for k in a.keys() & b.keys() if a[k] != b[k]]
print(f"records: {len(a)} vs {len(b)} | only in first: {len(a.keys()-b.keys())} | only in second: {len(b.keys()-a.keys())} | differing sequences: {len(diff)}")
if diff: print("first differing IDs:", sorted(diff)[:10])
print("RESULT:", "IDENTICAL" if not diff and a.keys() == b.keys() else "DIFFERENT")
