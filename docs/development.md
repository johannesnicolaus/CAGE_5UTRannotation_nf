# luscombeu/utrannotation: Development notes

## Status

The pipeline code is complete, but it has **not been run yet**. The test data are still being prepared, so the following are not done:

- `conf/test.config` and `conf/test_full.config` have `input = null`.
- The local modules have no `tests/main.nf.test` yet.
- `tests/default.nf.test` is the template's pipeline test and has no snapshot.

## Runs so far

- 2026-10-18, Osaka O10: O10 genome, O10.1.0 v1.3 GFF3 (AUGUSTUS models that already had CAGE UTRs), and cagescan `osa_MGI_220209_v4` (CAGEfightR clusters, separate SL clusters). Run directory `/flash/LuscombeU/nico/utrannotation_runs/osa_O10_v1.3`. Completed, 12 tasks.
  - 5' ends: 65% of 17,686 mRNAs got a CAGE TSS, 42% of them sl. 98.7% of sl TSSs follow an AG.
  - Proteins: all CDSs unchanged. 17,646 proteins are identical to the v1.3 peptides; the other 40 differ only in how internal stop codons are written (`*` vs `.`).
  - Operons: 2,496, with 43% of genes. Compared with the v1.3 operons, 84% of operonic genes are shared and 1,504 gene sets are identical.
  - The MultiQC length plots were not parsed (numeric first bin label); fixed afterwards.

- 2026-10-18, I69 EGAPx with the Okinawa CAGEr consensus clusters: `utrannotation_runs/oki_I69_egapx` (repackaged as cagescan output) and `oki_I69_egapx_bed` (BED input). Both completed with identical results. Also an O10 rerun with the new defaults. See rationale.md, section 10.

## Test data to add

A small, self-contained data set:

- a slice of the I25_prim genome (about 1 Mb of chr1) and the matching subset of the EGAPx `complete.genomic.gff`;
- a cagescan results directory for the slice: `tss_clusters/broad/all_samples.broad.tss_clusters.bed`, `..._tss_cluster_counts.tsv` and, for the SL case, `..._tss_cluster_sl_counts.tsv`.

## Cases the tests should cover

5' end (`ANNOTATE_UTR5`):

- [ ] plus and minus strand
- [ ] 5' UTR extended (TSS upstream of the EGAPx 5' end) and trimmed (TSS inside the EGAPx 5' UTR)
- [ ] `sl` and `nosl` TSSs; a dataset without `tss_cluster_sl_counts.tsv` (all `nosl`)
- [ ] `acceptor_AG` and the `trans_splice_acceptor_site` row for an SL TSS
- [ ] a weak exonic peak closer to the ATG ignored because of `--tss_min_rel`
- [ ] `--tss_selection dominant` versus `closest_atg`
- [ ] TSS in a 5' UTR intron ignored
- [ ] no candidate within `--tss_max_dist` (EGAPx 5' end kept)
- [ ] operon neighbour: no extension over the CDS of the upstream gene on the same strand
- [ ] 5'-partial model left alone

3' end (`ANNOTATE_UTR3`):

- [ ] model with `experiment=COORDINATES: polyA evidence` kept (and moved with `--keep_polya_evidence false`)
- [ ] PAS closest to the EGAPx end chosen over one nearer to the stop codon
- [ ] PAS inside the EGAPx 3' UTR (trim) and downstream of it (extend)
- [ ] PAS in a 3' UTR intron or straddling a splice site ignored
- [ ] no PAS within `--pas_window` (EGAPx 3' end kept)
- [ ] no extension into the CDS of the downstream gene on the same strand
- [ ] `--pas_offset`
- [ ] 3'-partial model left alone

5' end of non-coding transcripts:

- [ ] lncRNA / misc_RNA 5' end extended and trimmed; minimum remaining length
- [ ] `--tss_noncoding false`

Operons (`ANNOTATE_OPERONS`):

- [ ] Ganot et al. 2004 cyclin D3 operon (AJ628168) recovered as one operon. On OKI2018_I69 v2.3 with the old OKItoOKI consensus clusters, `bin/annotate_operons.py` gives PAR:15956186-15967604(-): POC5, IQCK, XPO7, DYNLL2, TUBD1, MBF1, CYCDC. PDCL3 is split off by a non-SL TSS of 783 reads.
- [ ] genome without SL counts: no operon step
- [ ] `--operons false` (SL labels kept, no operons) and `--trans_splicing false` (all nosl, no operons)
- [ ] several leaders (`tss_cluster_sl_leaders.tsv`): `sl_leader` attribute

Whole pipeline:

- [ ] proteins from `<id>.utr.proteins.fasta` identical to EGAPx `complete.proteins.faa`
- [ ] output GFF3 passes `gt gff3validator` and AGAT
- [ ] gzipped genome and annotation
- [ ] two genomes in one samplesheet
