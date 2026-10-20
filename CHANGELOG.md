# luscombeu/utrannotation: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of luscombeu/utrannotation, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- 5' UTRs from luscombeu/cagescan TSS clusters: the strongest mature-mRNA 5' end (`--tss_selection auto`; `closest_atg` and `dominant` are alternatives) within `--tss_max_dist`, with a relative-count noise filter. TSSs are labelled `sl` / `nosl` from the spliced-leader read fraction, and SL TSSs get a `trans_splice_acceptor_site` feature and an AG-acceptor check.
- 3' UTRs from an evidence hierarchy: EGAPx poly(A)-read evidence, then the polyadenylation signal closest to the EGAPx end within `--pas_window`, then the EGAPx end.
- Explicit `five_prime_UTR` / `three_prime_UTR` features with evidence attributes; GTF and transcript / protein FASTA of the result.
- UTR statistics before and after (MultiQC table and length distributions), and AGAT statistics.
- CAGE 5' ends for lncRNAs and misc_RNAs (strongest TSS within `--tss_max_dist` of the EGAPx 5' end; `--tss_noncoding`).
- `--tss_rel_whole_utr` (default on): the `--tss_min_rel` noise filter compares with the strongest TSS in the whole EGAPx 5' UTR.
- 3' ends no longer extend into the next same-strand transcript (operon neighbours).
- Input checks: annotation sequences must be in the genome, and TSSs must lie on genome sequences.
- `--trans_splicing` and `--operons` switches; `sl_leader` attribute when cagescan searched for several leaders; coding transcripts are any transcripts with a CDS (not only EGAPx `mRNA`), so AUGUSTUS-style GFF3 works too.
- SL acceptor rule (`--sl_window_up`, `--sl_window_down`, default 5): non-SL TSSs within 5 bp of a qualifying SL acceptor are ignored for the TSS choice and operon splits (SL reads without a recognisable leader). A transcript whose chosen TSS is non-SL but that also has a qualifying SL acceptor records it (`sl_acceptor_peak`, `sl_acceptor_count`).
- Defaults for metazoans in general: `--tss_selection auto` (strongest mature-mRNA 5' end), `--pas_motifs AATAAA,ATTAAA`, `--pas_offset 15`; `nosl_tss_peak` / `nosl_tss_count` record the non-SL TSS of SL transcripts.
- TSS input from any CAGE workflow: samplesheet `tss` (cagescan directory or BED) and `tss_sl` (BED of SL acceptor clusters) replace `cagescan_dir`.
- Operon gaps between codons by default (`--operon_gap_from cds`), with `transcript` and `evidence` modes (`--operon_max_gap_tx`).
- Operons (`ANNOTATE_OPERONS`) for genomes with spliced-leader CAGE data. Same-strand genes closer than `--operon_max_gap` (500, as in Plessy et al. 2024) are grouped, and a gene with its own non-SL TSS starts a new operon. On OKI2018_I69 with the Okinawa CAGE clusters, this recovers the Ganot et al. 2004 cyclin D3 operon (XPO7/RBP16, DYNLL2, TUBD1, MBF1, CYCDC) with the POC5 leader and IQCK; with CDS gaps (the default) this needs `--operon_max_gap` of at least 404.

- Operon split rule `--operon_split_rule` (default `dominant`): a non-SL TSS splits only when it is stronger than every SL TSS between the genes (`any`: the earlier rule; `none`: distance only). Operons record `sl_upstream` per downstream gene, `downstream_sl`, `operon_evidence` (`sl`, `partial_sl`, `distance`) and `sl_source`. When SL data contain no SL TSS, AG-preceded TSSs are taken as SL acceptors (`sl_source=AG_inferred`) instead of calling no operons.

### `Fixed`

- SL acceptor rule applied before the `--tss_min_rel` filter, so a strong non-SL shadow next to an acceptor cannot filter the acceptor out; `auto` now excludes non-SL TSSs upstream of any qualifying acceptor, as documented.
- 3' ends stop before every other-gene transcript starting downstream of the stop codon, also when the 5' step moved that transcript into the current 3' UTR.
- Attributes of an earlier run are removed when the output is annotated again; transcripts whose CDS lies outside their exons are skipped with a warning.
- Non-coding transcripts of pseudogenes and 5'-partial non-coding models are left unchanged.
- TSS BEDs: a thickStart outside the cluster falls back to the 5' end; duplicate cluster names get unique IDs; all-zero scores stop the run. Gzipped TSS BEDs are rejected at samplesheet validation. Relative samplesheet paths are resolved against the samplesheet.
- Operons: SL acceptors shield non-SL TSSs only when they pass `--operon_split_min_count`; no operons are called when the TSS table has no SL TSS; operons are numbered in annotation order.

### `Dependencies`

### `Deprecated`
