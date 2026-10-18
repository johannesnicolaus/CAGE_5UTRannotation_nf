# luscombeu/utrannotation: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of luscombeu/utrannotation, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- 5' UTRs from luscombeu/cagescan TSS clusters: the TSS closest to the ATG (or the dominant one) within `--tss_max_dist`, with a relative-count noise filter. TSSs are labelled `sl` / `nosl` from the spliced-leader read fraction, and SL TSSs get a `trans_splice_acceptor_site` feature and an AG-acceptor check.
- 3' UTRs from an evidence hierarchy: EGAPx poly(A)-read evidence, then the polyadenylation signal closest to the EGAPx end within `--pas_window`, then the EGAPx end.
- Explicit `five_prime_UTR` / `three_prime_UTR` features with evidence attributes; GTF and transcript / protein FASTA of the result.
- UTR statistics before and after (MultiQC table and length distributions), and AGAT statistics.
- CAGE 5' ends for lncRNAs and misc_RNAs (strongest TSS within `--tss_max_dist` of the EGAPx 5' end; `--tss_noncoding`).
- `--tss_rel_whole_utr` (default on): the `--tss_min_rel` noise filter compares with the strongest TSS in the whole EGAPx 5' UTR.
- 3' ends no longer extend into the next same-strand transcript (operon neighbours).
- Input checks: annotation sequences must be in the genome, and TSSs must lie on genome sequences.
- `--trans_splicing` and `--operons` switches; `sl_leader` attribute when cagescan searched for several leaders; coding transcripts are any transcripts with a CDS (not only EGAPx `mRNA`), so AUGUSTUS-style GFF3 works too.
- Operons (`ANNOTATE_OPERONS`) for genomes with spliced-leader CAGE data. Same-strand genes closer than `--operon_max_gap` (500, as in Plessy et al. 2024) are grouped, and a gene with its own non-SL TSS starts a new operon. On OKI2018_I69 with the Okinawa CAGE clusters, this recovers the Ganot et al. 2004 cyclin D3 operon (XPO7/RBP16, DYNLL2, TUBD1, MBF1, CYCDC) with the POC5 leader and IQCK, for maximum gaps from 100 to 750 bp.

### `Fixed`

### `Dependencies`

### `Deprecated`
