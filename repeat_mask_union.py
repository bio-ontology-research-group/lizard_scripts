#!/usr/bin/env python3
"""Bases covered by each RepeatMasker run and by their union.

Reads RepeatMasker .out files (query name, begin, end of every hit, including simple repeats
and low complexity) and prints, per file and for the union of all files, the number of
distinct bases covered. If the BRAKER3 input was soft-masked with the union of the runs,
the union total equals its count of lowercase bases (1,591,780,905).

usage: repeat_mask_union.py run1.out run2.out run3.out
"""
import sys, collections

def intervals(path):
    iv = collections.defaultdict(list)
    for line in open(path):
        f = line.split()
        if len(f) >= 15 and f[0].isdigit():
            iv[f[4]].append((int(f[5]), int(f[6])))
    return iv

def covered(iv):
    tot = 0
    for seq, v in iv.items():
        v.sort(); cs = ce = None
        for s, e in v:
            if cs is None or s > ce + 1:
                if cs is not None: tot += ce - cs + 1
                cs, ce = s, e
            elif e > ce:
                ce = e
        if cs is not None: tot += ce - cs + 1
    return tot

union = collections.defaultdict(list)
for p in sys.argv[1:]:
    iv = intervals(p)
    print(f"{p}\t{covered(iv)}")
    for k, v in iv.items():
        union[k].extend(v)
print(f"UNION\t{covered(union)}")
