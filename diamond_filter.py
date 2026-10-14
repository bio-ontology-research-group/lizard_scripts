#!/usr/bin/env python3
"""Select BRAKER3 transcripts by the identity of their DIAMOND hits against SwissProt.

Reads DIAMOND tabular output (the 12 default BLAST-tabular columns, no header) and keeps the
hits with percent identity >= --min-identity (default 60). Writes the kept hits, in the input
format, and the sorted list of distinct query IDs. That list is the keep-list of the final
gene set: a transcript is kept if it has at least one hit at or above the threshold.
The proteins of the kept transcripts can then be extracted with subset_fasta.py:

    python3 subset_fasta.py kept_ids.txt braker.aa braker_diamond_matched.aa

Requires pandas.

usage: diamond_filter.py DIAMOND_TSV --out-tsv KEPT_HITS.tsv --out-ids KEPT_IDS.txt [--min-identity 60]
"""
import argparse
import pandas as pd

COLUMNS = ['query_id', 'subject_id', 'perc_identity', 'alignment_length', 'mismatches',
           'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diamond_tsv", help="DIAMOND tabular output (12 columns, no header)")
    ap.add_argument("--out-tsv", required=True, help="kept hits, same format as the input")
    ap.add_argument("--out-ids", required=True, help="sorted distinct query IDs of the kept hits")
    ap.add_argument("--min-identity", type=float, default=60.0, help="minimum percent identity (default 60)")
    args = ap.parse_args()

    diamond = pd.read_csv(args.diamond_tsv, sep="\t", names=COLUMNS, dtype={'query_id': str, 'subject_id': str})
    kept = diamond[diamond['perc_identity'] >= args.min_identity]

    with open(args.diamond_tsv) as fin, open(args.out_tsv, "w") as fout:
        keep_rows = set(kept.index)
        for i, line in enumerate(l for l in fin if l.strip()):
            if i in keep_rows:
                fout.write(line)
    ids = sorted(set(kept['query_id']))
    with open(args.out_ids, "w") as o:
        o.write("\n".join(ids) + "\n")
    print(f"hits: {len(diamond)}; kept hits (identity >= {args.min_identity:g}): {len(kept)}; "
          f"queries: {diamond['query_id'].nunique()}; kept queries: {len(ids)}")

if __name__ == "__main__":
    main()
