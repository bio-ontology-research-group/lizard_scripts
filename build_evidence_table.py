#!/usr/bin/env python3
"""Build the per-transcript evidence table for all BRAKER3 transcripts of Scincus mitranus.

One row per transcript (mRNA) in the complete BRAKER3 annotation. Columns:
  transcript_id, gene_id, haplotype, sequence, start, end, strand
  in_filtered_set          1 if the transcript is in the final (SwissProt-filtered) set
  gene_in_filtered_set     1 if any transcript of its gene is in the final set
  interpro_entry           1 if InterProScan reports >=1 match integrated into an InterPro entry
  pfam_match               1 if PfamScan reports >=1 Pfam-A match (any type, gathering threshold)
  pfam_families            Pfam family names, comma-separated, sorted ('-' if none)
  pfam_te_only             1 if all of its Pfam families are transposable-element types
  transcript_evidence_model  1 if BRAKER3 built the model from transcript evidence (source 'gmst')

usage: build_evidence_table.py braker.gff3 braker.gtf interproscan.tsv pfam.txt keep_list.txt out.tsv
  braker.gff3 / braker.gtf: the renamed complete-set files (S3); interproscan.tsv may be '-' for stdin;
  pfam.txt: the cleaned PfamScan CSV (F1); keep_list.txt: final transcript IDs (D1).
"""
import sys, re, csv, collections

TE = re.compile(r"RVT|^rve|Transpos|DDE|Tnp|Retrotran|Exo_endo_phos|zf-H2C2|Integrase|RNase_H|MULE|"
                r"Dimer_Tnp|HAT|Tc1|Mariner|Helitron|Gag|gag|Peptidase_A2|HTH_Tnp|DUF4371|DUF1759|CENP-B|THAP", re.I)

def main(gff, gtf, ipr, pfam, keepf, out):
    keep = set(l.strip() for l in open(keepf) if l.strip())
    tx = {}                                   # transcript -> [gene, seq, start, end, strand]
    for l in open(gff):
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if f[2] == "mRNA":
            a = dict(kv.split("=", 1) for kv in f[8].strip(";").split(";") if "=" in kv)
            tx[a["ID"]] = [a["Parent"], f[0], f[3], f[4], f[6]]
    gmst = set()
    for l in open(gtf):
        f = l.split("\t")
        if len(f) > 8 and f[1] == "gmst":
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            gmst.add(m.group(1) if m else f[8].strip())
    interpro = set()
    fh = sys.stdin if ipr == "-" else open(ipr, encoding="utf-8", errors="replace")
    for l in fh:
        f = l.split("\t", 12)
        if len(f) > 11 and f[11] not in ("-", ""):
            interpro.add(f[0])
    fams = collections.defaultdict(set)
    for r in csv.reader(open(pfam, newline="")):
        if r and r[0] != "seq_id":
            fams[r[0]].add(r[6])
    kept_genes = {tx[t][0] for t in keep if t in tx}
    missing = keep - set(tx)
    if missing:
        sys.exit(f"ERROR: {len(missing)} keep-list IDs not in the GFF3, e.g. {sorted(missing)[:3]}")
    hap_order = lambda t: (tx[t][1].split(".")[0], int(re.search(r"\d+$", tx[t][1]).group()), int(tx[t][2]), t)
    with open(out, "w") as o:
        o.write("transcript_id\tgene_id\thaplotype\tsequence\tstart\tend\tstrand\tin_filtered_set\tgene_in_filtered_set\t"
                "interpro_entry\tpfam_match\tpfam_families\tpfam_te_only\ttranscript_evidence_model\n")
        for t in sorted(tx, key=hap_order):
            g, sq, s, e, st = tx[t]
            fs = sorted(fams.get(t, ()))
            te = int(bool(fs) and all(TE.search(x) for x in fs))
            o.write("\t".join(map(str, [t, g, sq.split(".")[0], sq, s, e, st, int(t in keep), int(g in kept_genes),
                                         int(t in interpro), int(bool(fs)), ",".join(fs) if fs else "-", te,
                                         int(t in gmst)])) + "\n")
    print(f"transcripts {len(tx)}\tin_filtered_set {len(keep)}\tinterpro {len(interpro & set(tx))}\t"
          f"pfam {sum(1 for t in tx if fams.get(t))}\tgmst {len(gmst & set(tx))}")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
