#!/usr/bin/env python3
"""Convert RepeatMasker .out files into per-haplotype GFF3 tracks with final sequence names.

One GFF3 feature per .out line. Feature type: 'microsatellite' for Simple_repeat,
'low_complexity_region' for Low_complexity, 'dispersed_repeat' for everything else.
Attributes: ID (unique), Name (repeat), repeat_class (class/family), rm_id (RepeatMasker element ID;
fragments of one element share it), repeat_begin / repeat_end / repeat_left (RepeatMasker's position-in-
repeat columns, verbatim), Target (consensus coordinates; written only when both are >= 1 and begin <= end,
because RepeatMasker can report positions <= 0 when an alignment runs past the consensus), divergence/
deletion/insertion (%), and overlapping_higher_score=true for lines RepeatMasker marks with '*'.
Column 6 is the Smith-Waterman score.

Sequence names are translated with a two-column map (e.g. scaffold_N -> hapN.chrM); names already
final pass through unchanged. The script stops on an unknown name or a coordinate beyond the
sequence length.

Also writes <prefix>.summary.tsv: per haplotype and in total, for each class: elements (distinct
RepeatMasker IDs), fragments (lines), bp (sum of fragment lengths) and bp_union (bases covered).

usage: repeatmasker_to_gff3.py seqid_map.tsv lengths.tsv OUT_PREFIX run1.out [run2.out ...]
  lengths.tsv: name<TAB>length[...] of the final assembly (fingerprint tables or .fai)
"""
import sys, re, collections

def esc(v):
    return re.sub(r"[;=&,%\t]", lambda m: "%{:02X}".format(ord(m.group())), v)

def so_type(cls):
    if cls.startswith("Simple_repeat"):
        return "microsatellite"
    if cls.startswith("Low_complexity"):
        return "low_complexity_region"
    return "dispersed_repeat"

def parse_out(path):
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue                       # header, blank or separator lines
            star = len(f) >= 16 and f[15] == "*"
            yield f, star

def union_bp(iv):
    iv.sort(); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None or s > ce + 1:
            if cs is not None: tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None: tot += ce - cs + 1
    return tot

def main(mapf, lenf, prefix, outs):
    nm = {}
    for l in open(mapf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and not a[0].startswith(("#", "braker_seqid")):
            nm[a[0]] = a[1]
    lens = {}
    for l in open(lenf):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2 and a[1].isdigit():
            lens[a[0]] = int(a[1])
    feats = collections.defaultdict(list)
    no_target = 0; examples = []
    for k, path in enumerate(outs, 1):
        n = 0
        for f, star in parse_out(path):
            q = nm.get(f[4], f[4])
            if q not in lens:
                sys.exit(f"ERROR: {path}: sequence {f[4]} not in the map or the final assembly")
            s, e = int(f[5]), int(f[6])
            if not (1 <= s <= e <= lens[q]):
                sys.exit(f"ERROR: {path}: {f[4]}:{s}-{e} outside 1-{lens[q]}")
            strand = "+" if f[8] == "+" else "-"
            rep, cls = f[9], f[10]
            if strand == "+":
                rb, re_, rl = f[11], f[12], f[13]
            else:
                rl, re_, rb = f[11], f[12], f[13]
            rid = f"r{k}_{f[14]}"
            n += 1
            fid = f"{rid}_{n}"
            attrs = [f"ID={fid}", f"Name={esc(rep)}", f"repeat_class={esc(cls)}", f"rm_id={rid}",
                     f"repeat_begin={esc(rb)}", f"repeat_end={esc(re_)}", f"repeat_left={esc(rl)}"]
            try:
                tb, tz = int(rb), int(re_)
            except ValueError:
                tb = tz = 0
            if 1 <= tb <= tz:
                attrs.append(f"Target={esc(rep)} {tb} {tz} {strand}")
            else:
                no_target += 1
                if len(examples) < 5: examples.append(" ".join(f))
            attrs += [f"divergence={f[1]}", f"deletion={f[2]}", f"insertion={f[3]}"]
            if star:
                attrs.append("overlapping_higher_score=true")
            feats[q].append((s, e, f[0], strand, so_type(cls), ";".join(attrs), cls, rid))
    key = lambda n: (n.split(".")[0], 0 if ".chr" in n else 1, int(re.search(r"\d+$", n).group()))
    stats = collections.defaultdict(lambda: collections.defaultdict(lambda: [set(), 0, 0, []]))
    for hap in sorted({q.split(".")[0] for q in feats}):
        seqs = sorted((q for q in feats if q.startswith(hap + ".")), key=key)
        with open(f"{prefix}.{hap}.gff3", "w") as o:
            o.write("##gff-version 3\n")
            for q in seqs:
                o.write(f"##sequence-region {q} 1 {lens[q]}\n")
            for q in seqs:
                for s, e, sc, st, ty, at, cls, rid in sorted(feats[q], key=lambda x: (x[0], x[1])):
                    o.write(f"{q}\tRepeatMasker\t{ty}\t{s}\t{e}\t{sc}\t{st}\t.\t{at}\n")
                    top = cls.split("/")[0]
                    for h in (hap, "total"):
                        r = stats[h][top]; r[0].add(rid); r[1] += 1; r[2] += e - s + 1; r[3].append((hash(q), s, e))
    with open(f"{prefix}.summary.tsv", "w") as o:
        o.write("haplotype\tclass\telements\tfragments\tbp\tbp_union\n")
        for h in sorted(stats, key=lambda x: (x == "total", x)):
            allv = []
            for c in sorted(stats[h]):
                el, fr, bp, iv = stats[h][c]
                byq = collections.defaultdict(list)
                for qh, s, e in iv: byq[qh].append((s, e))
                u = sum(union_bp(v) for v in byq.values()); allv.extend(iv)
                o.write(f"{h}\t{c}\t{len(el)}\t{fr}\t{bp}\t{u}\n")
            byq = collections.defaultdict(list)
            for qh, s, e in allv: byq[qh].append((s, e))
            o.write(f"{h}\tALL\t-\t-\t-\t{sum(union_bp(v) for v in byq.values())}\n")
    print(f"sequences with repeats: {len(feats)}; features: {sum(len(v) for v in feats.values())}; "
          f"without Target (repeat positions < 1 or reversed): {no_target}")
    for e in examples:
        print("  example .out line:", e)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
