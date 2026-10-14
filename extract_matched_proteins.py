#!/usr/bin/env python3
"""Extract the proteins of the transcripts kept by diamond_filter.py.

Reads the transcript IDs from the first column of IDS (either the kept-hits table written with
--out-tsv or the ID list written with --out-ids by diamond_filter.py), and writes the BRAKER3
proteins whose ID is in that set, in the order of the input FASTA. This is the step that produced
braker_diamond_matched.aa, the protein file of the final gene set.

Requires Biopython.

usage: extract_matched_proteins.py IDS braker.aa braker_diamond_matched.aa
"""
import sys
from Bio import SeqIO

def main(ids_file, input_fasta, output_fasta):
    diamond_proteins = set()
    with open(ids_file) as fh:
        for line in fh:
            if line.strip():
                diamond_proteins.add(line.split("\t")[0].strip())
    records = [record for record in SeqIO.parse(input_fasta, "fasta")
               if record.id in diamond_proteins]
    SeqIO.write(records, output_fasta, "fasta")
    missing = len(diamond_proteins - {r.id for r in records})
    print(f"Wrote {len(records)} records to {output_fasta} ({missing} requested IDs not found in {input_fasta})")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
