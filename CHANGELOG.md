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

### `Fixed`

### `Dependencies`

### `Deprecated`
