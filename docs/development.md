# luscombeu/utrannotation: Development notes

## Status

The pipeline code is complete, but it has **not been run yet**. The test data are still being prepared, so the following are not done:

- `conf/test.config` and `conf/test_full.config` have `input = null`.
- The local modules have no `tests/main.nf.test` yet.
- `tests/default.nf.test` is the template's pipeline test and has no snapshot.

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
- [ ] `--skip_operons`

Whole pipeline:

- [ ] proteins from `<id>.utr.proteins.fasta` identical to EGAPx `complete.proteins.faa`
- [ ] output GFF3 passes `gt gff3validator` and AGAT
- [ ] gzipped genome and annotation
- [ ] two genomes in one samplesheet
