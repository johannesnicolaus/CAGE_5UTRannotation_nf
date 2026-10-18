# luscombeu/utrannotation: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Annotation](#annotation) - the GFF3/GTF with CAGE 5' UTRs and refined 3' UTRs, and its sequences
- [TSS table](#tss-table) - cagescan TSS clusters labelled sl / nosl
- [Polyadenylation signals](#polyadenylation-signals) - motif positions in the genome
- [Statistics](#statistics) - UTR lengths and evidence before and after, AGAT statistics
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

`<id>` below is the `id` column of the samplesheet.

### Annotation

<details markdown="1">
<summary>Output files</summary>

- `annotation/`
  - `<id>.utr.gff3`: the annotation with new 5' and 3' ends and explicit UTR features (main result).
  - `<id>.utr.gtf`: the same in GTF (AGAT).
  - `<id>.utr.transcripts.fasta`: spliced transcript sequences (gffread `-w`).
  - `<id>.utr.proteins.fasta`: protein sequences (gffread `-y`). They must be identical to the EGAPx proteins, because CDSs are never changed.
  - `intermediate/<id>.utr5.log`, `intermediate/<id>.utr3.log`: number of mRNAs per outcome of each step.

</details>

What changes from the EGAPx GFF3:

- **5' end:** the 5' end of coding `mRNA`s and their first exon moves to the chosen TSS. Exons upstream of it are removed. `gene` rows are resized to their transcripts.
- **3' end:** the 3' end and last exon move to the polyadenylation signal, unless EGAPx had poly(A)-read evidence.
- **UTR rows:** `five_prime_UTR` and `three_prime_UTR` rows are added for every coding mRNA, one per exonic segment. Their IDs are `utr5-<transcript>-<n>` and `utr3-<transcript>-<n>`, numbered 5' to 3'.
- **SL marker:** a TSS labelled `sl` adds a 1-bp `trans_splice_acceptor_site` row (SO:0000706, source `CAGE`) at the acceptor, with Parent the mRNA.
- **TSS at the ATG:** when the TSS is the first base of the CDS (common for SL acceptors), the genomic 5' UTR has length 0. Such an mRNA has no `five_prime_UTR` row but still carries `utr5_evidence=CAGE`, `tss_type` and the other TSS attributes. Its mature 5' UTR is the spliced leader, which is not in the genome.
- **Unchanged:** CDS rows, IDs, lncRNAs, pseudogenes and all EGAPx attributes. The original exon IDs are kept, even when exons were removed.

The pipeline's evidence attributes are on the `mRNA` row only, because they describe the transcript; the `five_prime_UTR` and `three_prime_UTR` rows carry just `ID` and `Parent`:

| Attribute       | Meaning                                                                                       |
| --------------- | --------------------------------------------------------------------------------------------- |
| `utr5_evidence` | `CAGE` (5' end at a CAGE TSS) or `EGAPx` (no qualifying TSS; EGAPx 5' end kept)               |
| `tss_type`      | `sl` (spliced-leader acceptor, trans-spliced) or `nosl`                                       |
| `tss_id`        | cagescan TSS cluster (`TC…` broad, `TS…` sharp)                                               |
| `tss_peak`      | dominant CTSS, `chrom:position:strand` (1-based); the new 5' end                              |
| `tss_count`     | READ1s of the cluster, summed over samples                                                    |
| `sl_count`      | of these, READ1s carrying the spliced leader                                                  |
| `sl_fraction`   | `sl_count / tss_count`                                                                        |
| `acceptor_AG`   | `sl` only: `true` if the two genomic bases before the TSS are AG (the trans-splice acceptor) |
| `utr5_shift`    | new minus EGAPx 5' UTR length (positive: extended, negative: trimmed)                         |
| `utr3_evidence` | `polyA_reads`, `PAS` or `EGAPx` (see [usage](usage.md#3-end-annotate_utr3))                   |
| `pas_motif`     | the polyadenylation signal used                                                               |
| `pas_pos`       | its position, `chrom:start-end:strand` (1-based)                                              |
| `utr3_shift`    | new minus EGAPx 3' UTR length                                                                 |

### TSS table

<details markdown="1">
<summary>Output files</summary>

- `tss/`
  - `<id>.tss.tsv`: one row per cagescan TSS cluster. The columns are:
    - `chrom`, `peak` (dominant CTSS, 1-based), `strand`, `tss_id`;
    - `count`, `sl_count` (summed over samples);
    - `sl_fraction`, `tss_type` (`sl` / `nosl`).

</details>

### Polyadenylation signals

<details markdown="1">
<summary>Output files</summary>

- `pas/`
  - `<id>.pas.bed`: every occurrence of `--pas_motifs` in the genome, as BED6. The fields are:
    - 0-based start;
    - name = motif;
    - strand = strand of an RNA that would carry the motif.

</details>

### Statistics

<details markdown="1">
<summary>Output files</summary>

- `stats/`
  - `<id>.utr_per_transcript.tsv`: one row per coding mRNA, with 5' and 3' UTR lengths in EGAPx and in the output, plus the evidence and TSS / PAS attributes.
  - `<id>.utr_summary_mqc.tsv`: the summary table shown in MultiQC. It gives:
    - % of mRNAs with a CAGE TSS, % sl, % sl with AG;
    - extended and trimmed UTRs;
    - median and mean UTR lengths before and after;
    - the share of 3' ends set by each rule.
  - `<id>.utr5_length_mqc.tsv`, `<id>.utr3_length_mqc.tsv`: UTR length distributions before and after, shown in MultiQC.
- `stats/agat/`
  - `<id>.egapx.stats.txt`, `<id>.utr.stats.txt`: `agat_sp_statistics.pl` of the input and output annotations.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate the UTR summary table and UTR length distributions of each annotation. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
