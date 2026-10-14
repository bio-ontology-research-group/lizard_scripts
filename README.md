# Scripts for the *Scincus mitranus* genome assembly and annotation

Custom scripts used to produce, check and submit the assembly and annotation files of the Eastern sandfish (*Scincus mitranus*, isolate OM-2025) genome, described in the accompanying Scientific Data Data Descriptor.

- Assembly files: figshare DOI 10.6084/m9.figshare.33028703
- Annotation files: figshare DOI 10.6084/m9.figshare.31860841
- GenBank: BioProjects PRJNA1348126 (pseudohaplotype 1) and PRJNA1348127 (pseudohaplotype 2)

The Python scripts need Python 3.8 or later and only the standard library, except `diamond_filter.py` (needs pandas) and `extract_matched_proteins.py` (needs Biopython). External tools are named where a step uses them. Each script prints its usage when run without arguments.

## Scripts that produce deposited files

| Script | Output | Main inputs |
|---|---|---|
| `diamond_filter.py` | Keep-list of the final gene set (`final_transcripts_keep_list.txt`) and the kept DIAMOND hits | DIAMOND best hits of the BRAKER3 proteins against SwissProt; minimum identity 60% |
| `extract_matched_proteins.py` | `braker_diamond_matched.aa` (proteins of the final gene set; identical to `filtered_proteins.faa`) | kept IDs or kept hits from `diamond_filter.py`, BRAKER3 `braker.aa` |
| `make_final_agp.py` | `hap1_final.agp`, `hap2_final.agp` | `out_JBAT.FINAL.agp`, `data/jbat_to_final.tsv`, N-run table of the final FASTAs (`seqtk cutN -gn 1`) |
| `agp_to_fasta.py` | `haplotype1.fa`, `haplotype2.fa` (final pseudohaplotype FASTAs; compare with `fingerprint.py`) | `hap1_final.agp` / `hap2_final.agp` and the Hifiasm contig FASTAs (all in the assembly item) |
| `make_chromosome_names.py` | `chromosome_names.tsv` | `data/seqid_map.tsv`, `data/source_hap1.tsv`, `data/source_hap2.tsv`, sequence lengths |
| `build_final_gff3.py` | `filtered_annotation.gff3` (with `final_transcripts_keep_list.txt`); renamed complete `braker.gff3` (without a keep-list) | BRAKER3 `braker.gff3`, `data/seqid_map.tsv`, sequence lengths |
| `subset_fasta.py` | `filtered_proteins.faa`, `filtered_cds.fna` | `braker.aa`, `braker.codingseq`, keep-list |
| `build_evidence_table.py` | `transcript_evidence.tsv` | renamed `braker.gff3` and `braker.gtf`, `interproscan.tsv`, `pfam.txt`, keep-list |
| `clean_kegg_pathways.py` | `KEGG_reconstruct_pathways_filtered.tsv` | KEGG pathway CSV written by kegg_parser |
| `make_products.py` | `products.tsv` | DIAMOND best hits against SwissProt, SwissProt titles (`blastdbcmd`), keep-list |
| `build_ncbi_gff3.py` | GenBank submission GFF3 per pseudohaplotype; `locus_map_hap1.tsv`, `locus_map_hap2.tsv` | per-pseudohaplotype final GFF3, locus-tag prefix (ACVTY0, ACVTY1), `products.tsv` |
| `repeatmasker_to_gff3.py` | `repeats_custom.hap1.gff3`, `repeats_custom.hap2.gff3`, summary | RepeatMasker `.out` (custom library), `data/seqid_map.tsv`, sequence lengths |
| `softmask_to_bed.py` | `softmask_union.bed` (and a per-sequence check table) | soft-masked chromosomes used as BRAKER3 input, `data/seqid_map.tsv` |

## GenBank submission

`run_table2asn.sh` runs table2asn for both pseudohaplotypes with the files in `data/` (chromosome names in `source_hap*.tsv`, gap linkage evidence in `linkage_hap*.txt`). It checks every input by MD5 before running, and needs internet access for the taxonomy lookup.

## Verification scripts

| Script | Checks |
|---|---|
| `fingerprint.py` | Per-sequence length and MD5 of the uppercased sequence (compares sequences ignoring soft-masking) |
| `compare_fasta.py` | Two FASTA files by sequence, e.g. `gffread` output against the deposited protein and CDS files |
| `repeat_mask_union.py` | Bases covered by each RepeatMasker run and by their union (the gene-prediction mask) |

One-line commands used between these steps are in `COMMANDS.md`.

## Input files

Files deposited in the figshare items are not repeated here. From the assembly item: the final FASTAs, `hap1_final.agp`, `hap2_final.agp`, the Hifiasm contig FASTAs, `out_JBAT.FINAL.agp` and `chromosome_names.tsv`. From the annotation item: `final_transcripts_keep_list.txt`, the BRAKER3 outputs, `interproscan.tsv`, `pfam.txt` and the other annotation files. Tool outputs that are not deposited (DIAMOND hits, SwissProt titles, kegg_parser CSV, RepeatMasker `.out`, the soft-masked BRAKER3 input) are produced by the steps described in the article.

## Data tables

Small tables used by the scripts and not deposited elsewhere.

| File | Content |
|---|---|
| `data/jbat_to_final.tsv` | Juicebox object → final sequence name (`-`: dropped mitochondrial contig) |
| `data/seqid_map.tsv` | Juicebox/BRAKER3 scaffold name → final chromosome name |
| `data/source_hap1.tsv`, `data/source_hap2.tsv` | GenBank chromosome names given to table2asn |
| `data/linkage_hap1.txt`, `data/linkage_hap2.txt` | Gap lengths with align-genus linkage evidence (gaps introduced by Hifiasm) |