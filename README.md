# luscombeu/utrannotation


[![GitHub Actions CI Status](https://github.com/luscombeu/utrannotation/actions/workflows/nf-test.yml/badge.svg)](https://github.com/luscombeu/utrannotation/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/luscombeu/utrannotation/actions/workflows/linting.yml/badge.svg)](https://github.com/luscombeu/utrannotation/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/luscombeu/utrannotation)

## Introduction

**luscombeu/utrannotation** adds UTRs to [EGAPx](https://github.com/ncbi/egapx) genome annotations. It uses CAGE TSS clusters from [luscombeu/cagescan](https://github.com/johannesnicolaus/cagescan_nextflow) for the 5' ends and polyadenylation evidence for the 3' ends. It was written for *Oikopleura dioica*, where a large share of genes are trans-spliced to a spliced leader (SL) and genes are tightly packed.

EGAPx models carry no TSS or trans-splicing information. Without long reads, their UTRs come from short-read coverage only and are often over-extended. For each coding mRNA, the pipeline does the following:

1. **5' end: moved to a CAGE TSS.** lncRNAs and misc_RNAs also get a CAGE 5' end, using the strongest TSS near their EGAPx 5' end.
   - It uses the dominant CTSS of a cagescan TSS cluster within `--tss_max_dist` (500) bases of the ATG.
   - By default it takes the TSS closest to the ATG among peaks at least `--tss_min_rel` (10%) as strong as the strongest one in the window. For trans-spliced genes, this is the SL acceptor rather than the outron TSS. `--tss_selection dominant` takes the strongest peak instead.
   - The TSS is labelled `sl` (spliced-leader acceptor) or `nosl` from the fraction of SL reads in the cluster. Datasets without an SL table are all `nosl`.
   - SL TSSs also get a `trans_splice_acceptor_site` feature and an AG-acceptor check.
2. **3' end: set from the best available evidence.**
   - (i) It keeps the EGAPx end when EGAPx found poly(A)-tail reads.
   - (ii) Otherwise it moves to the polyadenylation signal (`AATAAA`) closest to the EGAPx end, within `--pas_window` (500) bases.
   - (iii) Otherwise it keeps the EGAPx end.
3. **Operons.** For genomes with spliced-leader CAGE data, same-strand genes closer than `--operon_max_gap` (500) bases are grouped into operons. A gene with its own non-SL TSS starts a new operon, however close it is to the previous gene. Genomes without a spliced leader get UTRs only.
4. **UTR features.** Explicit `five_prime_UTR` and `three_prime_UTR` features are written, with evidence attributes. CDSs are never changed.

This pipeline is based on the CAGE 5' UTR annotation workflow in [LuscombeU_5primeanno](https://github.com/johannesnicolaus/LuscombeU_5primeanno). That workflow added SL and non-SL CAGE 5' UTRs to *O. dioica* gene models that had none (OKI2018_I69, OSKA2016 and Bar2). This pipeline brings the same approach to EGAPx annotations, which already have UTRs, as an nf-core pipeline. The `five_prime_UTR` naming and the `tss_type=sl|nosl` labels come from that workflow.

Other code it builds on:
- The 5' and 3' rules are adapted from the [annotationpolish](https://github.com/oist/LuscombeU_annotationpolish) R functions `add_five_prime_utr()` and `add_three_prime_utr()`. Those functions only extended UTRs; this pipeline can also trim the over-extended EGAPx UTRs.
- The operon rules come from annotationpolish's `calcOperons()` and `split_operons()`.

Steps:

1. Uncompress the inputs and find the TSS cluster files in the cagescan results directory.
2. Build a TSS table with SL labels (`CAGE_TSS_TABLE`).
3. Scan the genome for polyadenylation signals (`SCAN_PAS`).
4. Annotate 5' UTRs (`ANNOTATE_UTR5`), then 3' UTRs (`ANNOTATE_UTR3`), then operons (`ANNOTATE_OPERONS`).
5. Compute UTR statistics before and after (`UTR_STATS`, [`AGAT`](https://github.com/NBISweden/AGAT)).
6. Produce the GTF ([`AGAT`](https://github.com/NBISweden/AGAT)) and transcript and protein FASTA ([`gffread`](https://github.com/gpertea/gffread)).
7. Build the report ([`MultiQC`](http://multiqc.info/)).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow.

First, run [luscombeu/cagescan](https://github.com/johannesnicolaus/cagescan_nextflow) on the CAGE libraries against the **same genome** as the annotation. To get SL labels, pass `--sl_sequence` for organisms with trans-splicing.

Then prepare a samplesheet with one row per genome annotation:

`samplesheet.csv`:

```csv
id,genome,annotation,cagescan_dir
I25_prim,I25_primary_assembly_v1.0.0.fa,egapx/I25_prim/complete.genomic.gff,cagescan/I25_okinawa/results
```

`cagescan_dir` is the cagescan `--outdir`, the directory that contains `tss_clusters/`.

Now, you can run the pipeline using:

```bash
nextflow run luscombeu/utrannotation \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

See [usage](docs/usage.md) for the parameters and [output](docs/output.md) for the results and GFF3 attributes.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

luscombeu/utrannotation was originally written by Johannes Nicolaus Wibisana.

It is based on the 5' UTR annotation workflow in [LuscombeU_5primeanno](https://github.com/johannesnicolaus/LuscombeU_5primeanno), and its UTR and operon rules come from [annotationpolish](https://github.com/oist/LuscombeU_annotationpolish); both are Luscombe Unit (OIST) code.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
