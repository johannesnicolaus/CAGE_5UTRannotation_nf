# luscombeu/utrannotation: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline takes, for each genome:
- an [EGAPx](https://github.com/ncbi/egapx) annotation (`complete.genomic.gff`);
- the genome FASTA it was made on;
- CAGE TSS clusters on the **same genome**.

The TSS clusters can be either:
- the results of [luscombeu/cagescan](https://github.com/johannesnicolaus/cagescan_nextflow);
- BED files from any CAGE workflow, single- or paired-end, such as nf-core/cageseq, CAGEr or CAGEfightR.

The sequence names must match across all inputs; EGAPx keeps the names of the input genome.

## Samplesheet input

The samplesheet is a comma-separated file with a header row and one row per genome annotation. Several genomes can be processed in one run, for example the Okinawa, mainland and Europe assemblies.

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
id,genome,annotation,tss,tss_sl
I25_prim,/path/to/I25_primary_assembly_v1.0.0.fa,/path/to/egapx_runs/annotation_egapx/I25_prim/complete.genomic.gff,/path/to/cagescan/I25_okinawa/results,
other,/path/to/genome.fa,/path/to/complete.genomic.gff,/path/to/tss_clusters.bed,/path/to/sl_acceptor_clusters.bed
```

| Column       | Description |
| ------------ | ----------- |
| `id`         | Name of the annotation, used as the file name prefix of every output. Letters, digits, `_`, `-` and `.` only; each id once. |
| `genome`     | Genome FASTA (`.fa`, `.fasta` or `.fna`, optionally gzipped). |
| `annotation` | EGAPx GFF3 (`.gff` or `.gff3`, optionally gzipped). Other GFF3 files with gene → transcript → exon/CDS rows work too, if the CDS includes the stop codon, as in EGAPx output. |
| `tss`        | Either a cagescan `--outdir` (the directory that contains `tss_clusters/`) or a BED file of TSS clusters (see below). |
| `tss_sl`     | Optional, BED input only: BED of spliced-leader (trans-splice acceptor) clusters, built from the CAGE reads that carried the leader. Leave it empty for organisms without trans-splicing. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### TSS clusters as BED files (any CAGE workflow)

The BED files must be BED6 or longer, with these columns:
- **chrom, start, end:** the cluster coordinates.
- **name:** the cluster ID. `.` is replaced by the position.
- **score:** the count of the cluster (reads or TPM), used by the count filters.
- **strand:** `+` or `-`.

The TSS position, or peak, is taken in this order:
1. `thickStart`, when the BED has at least 8 columns. CAGEr, CAGEfightR and cagescan write the dominant CTSS there.
2. Otherwise the cluster itself, if it is 1 bp wide (a CTSS).
3. Otherwise the cluster's 5' end, with a warning.

Clusters in `tss` are non-SL (`nosl`) and clusters in `tss_sl` are SL (`sl`). Separate SL and non-SL clusters are what the CAGEr workflow of oist/LuscombeU-CAGE_libraries produces (`consensus_clusters_{no,sl}_*.bed`). SL reads can be separated from the others before clustering, for example with TagDust or cutadapt on READ1.

### What is read from a cagescan results directory

The pipeline reads cagescan results from v1.0.0dev onwards, with either kind of broad cluster set: CAGEfightR or paraclu clusters, with SL and non-SL reads clustered separately or pooled.

Files are taken from `tss_clusters/<set>/`, where `<set>` is `--cluster_set` (`broad`, the default, or `sharp`):

- `all_samples.<set>.tss_clusters.bed`: the clusters. The dominant CTSS (`thickStart`) is the TSS position.
- `all_samples.<set>.tss_cluster_counts.tsv`: READ1 counts per cluster and sample, summed over samples.
- `all_samples.<set>.tss_cluster_sl_counts.tsv`: READ1s that carried the spliced leader. cagescan writes this file only when `--sl_sequence` was given and SL reads were found. Without it, every TSS is labelled `nosl`.
  - With cagescan's separate SL clusters (`--sl_separate_clusters`, its default), `SL…` clusters hold only SL reads, so they are `sl`. `TC…` clusters hold none, so they are `nosl`. They may overlap.
- `all_samples.<set>.tss_cluster_sl_leaders.tsv` (optional, several leaders only): the main leader of each SL cluster.

A cluster is labelled `sl` when SL reads make up at least `--sl_min_fraction` (default 0.5) of its reads. The fraction and the SL count are kept as attributes, so TSSs used by both trans-spliced and non-trans-spliced RNAs can be found.

## How UTRs are assigned

Coding `mRNA`s get new 5' and 3' ends. Non-coding transcripts (lncRNAs and EGAPx misc_RNA `transcript` features) get a new 5' end only. Pseudogenes are copied unchanged.

Before anything is changed, the pipeline checks the inputs. It stops if annotation sequences are missing from the genome, or if TSSs lie on sequences or positions the genome does not have; both mean cagescan was run on another assembly.

### 5' end of coding mRNAs (`ANNOTATE_UTR5`)

1. **Candidates.** TSS peaks on the transcript's strand qualify if:
   - they give a 5' UTR of `--utr5_min_len` to `--tss_max_dist` bases, counted along the transcript (5' UTR exons, then the genome upstream of the EGAPx 5' end);
   - they are not in a 5' UTR intron;
   - the new 5' end does not overlap a CDS of another gene on the same strand.
2. **Noise filter.** Peaks below `--tss_min_count` are dropped, and so are peaks below `--tss_min_rel` × the strongest reference peak. CAGE "exon painting" (recapping or degradation) leaves weak peaks inside exons, and these would otherwise be picked as closest to the ATG.
   - By default the reference is every peak in the EGAPx 5' UTR, also beyond `--tss_max_dist`, plus the candidates. Without this, a weak peak near the ATG could trim a UTR whose strong TSS lies further away; such mRNAs now keep the EGAPx 5' end.
   - `--tss_rel_whole_utr false` compares with the candidates only.
3. **Choice.**
   - `--tss_selection auto` (default) takes the strongest mature-mRNA 5' end. SL acceptors, and non-SL TSSs that are not upstream of an acceptor, compete by count.
     - Without SL data this is the dominant TSS, the usual CAGE choice.
     - With SL data, a non-SL peak upstream of an SL acceptor cannot win. It may be the outron start, whose reads come from short-lived pre-mRNA.
   - `--tss_selection closest_atg` takes the shortest UTR. In trans-spliced genes this is the SL acceptor, and a weaker non-SL TSS downstream of an acceptor wins over it.
   - `--tss_selection dominant` takes the highest count, SL or not.
   - Ties go to the shorter UTR, then to `nosl`.
   - The strongest TSS of the other kind is recorded: `sl_acceptor_peak` for a non-SL TSS, `nosl_tss_peak` for an SL TSS.
4. **Edit.** Exons upstream of the TSS are dropped and the first exon starts at the TSS, so the EGAPx 5' UTR is **extended or trimmed**. Gene bounds follow. CDSs are never changed.
5. **Unchanged models.** Models without a qualifying TSS, and models open at the 5' end (EGAPx `partial`), keep the EGAPx 5' end (`tss_evidence=EGAPx`).
6. **SL acceptor rule.** When an SL peak passes the count filters (a "qualifying acceptor"), non-SL peaks up to `--sl_window_up` (5) bases upstream or `--sl_window_down` (5) bases downstream of it are not TSS candidates. They are taken as SL reads without a recognisable leader.
   - **Osaka data:** 104 of the 142 non-SL clusters within 5 bp of an SL cluster sit exactly 2 nt upstream, on the A of the acceptor AG. At the two sites inspected, their reads keep only the leader's final AG and come almost only from the DE1 and DE2 libraries.
   - **Why only 5 bp:** further away, non-SL TSSs can be genuine in tunicates.
     - In *Ciona*, 38 of 107 genes that make both trans-spliced and non-trans-spliced RNAs have the non-trans-spliced TSS just upstream of the acceptor in the same exon, 36 bp on average (Matsumoto et al. 2010).
     - Some *Ciona* TSS clusters even overlap acceptors (Yokomori et al. 2016).
     - In *O. dioica* testis, 89 male-specific promoters start transcription about 40-50 nt downstream of a trans-splice site, giving mRNAs without the acceptor (Danks et al. 2018).
   - **Wider windows:** `--sl_window_up 50` applies the minimum outron length of nematode experiments (Conrad et al. 1995). `--sl_window_down -1` ignores every non-SL peak inside the mature trans-spliced mRNA.
   - **What is kept when a non-SL TSS wins:** when the chosen TSS is non-SL but the transcript also has a qualifying SL acceptor, that acceptor is recorded (`sl_acceptor_peak`, `sl_acceptor_count`). The transcript is then marked as also making trans-spliced RNA.
   - **Effect on O10:** the defaults label 11 more mRNAs `sl`, and 172 `nosl` mRNAs record an SL acceptor.

   The reasoning and the evidence for every rule are collected in [rationale.md](rationale.md).
7. **TSS at the ATG.** A TSS may sit at the ATG itself (`--utr5_min_len 0`), as SL acceptors often do. The mRNA then has no `five_prime_UTR` row but still carries the TSS attributes.

### 5' end of non-coding transcripts

lncRNAs and misc_RNAs have no ATG, so the rule is different:
- **Candidates:** peaks within `--tss_max_dist` of the EGAPx 5' end. Upstream peaks count along the genome; downstream peaks count along the exons and must leave at least 100 bases of transcript. The same intron and CDS rules apply.
- **Choice:** the peak with the highest count, with ties going to the smallest move.
- **Labels:** `tss_type`, the acceptor row and the AG check work as for mRNAs.
- **Switch:** `--tss_noncoding false` leaves non-coding transcripts unchanged.

### 3' end (`ANNOTATE_UTR3`)

The first rule that applies sets the 3' end:

1. **`polyA_reads`.** EGAPx marked the model `experiment=COORDINATES: polyA evidence`: RNA-seq reads that contain the poly(A) tail mapped the cleavage site. The EGAPx end is kept. Switch this off with `--keep_polya_evidence false`.
2. **`PAS`.** A polyadenylation signal (`--pas_motifs`, default `AATAAA,ATTAAA`, the two most common metazoan signals) on the transcript's strand qualifies if it:
   - lies after the stop codon;
   - is within `--pas_window` bases of the EGAPx 3' end, measured along the transcript;
   - is not in a 3' UTR intron;
   - does not reach a CDS of another gene on the same strand;
   - does not extend the transcript into another gene's transcript on the same strand. In operons, the 3' end stops before the next gene's 5' end, which the 5' step has already set. Overlaps already present in EGAPx are left alone.

   The motif **closest to the EGAPx end** is used, so short-read coverage places the end and the motif sets the base. The new end is the motif end + `--pas_offset` (15). The signal lies within ~40 nt upstream of the cleavage site (Tian & Graber 2012). The 3' UTR is trimmed or extended.
3. **`EGAPx`.** No signal qualified, and the EGAPx 3' end is kept.

Hexamer motifs also occur by chance (about 0.6 AATAAA per kb per strand at 41% GC). `utr3_evidence` and the UTR statistics show how many ends each rule set.

## Operons (`ANNOTATE_OPERONS`)

Trans-splicing and operons have separate switches, because one does not imply the other. Eukaryotic operons depend on SL trans-splicing, which gives the downstream genes of a polycistronic transcript their cap. But trans-splicing does not imply operons: in every trans-splicing animal many SL genes are monocistronic, and operons are common in some lineages (nematodes, tunicates) and rare or unknown in others.

| Organism | Settings | Result |
| --- | --- | --- |
| With trans-splicing and operons (*Oikopleura*, *C. elegans*) | defaults | `sl` / `nosl` TSSs and operons |
| With trans-splicing, without operons | `--operons false` | `sl` / `nosl` TSSs, no operons |
| Without trans-splicing | defaults | cagescan writes no SL counts, so every TSS is `nosl` and no operons are annotated |

`--trans_splicing false` ignores the SL counts even when they exist: every TSS is `nosl` and there are no acceptor rows and no operons. Operons are never annotated without SL data.

When cagescan searched for several leaders (`--sl_sequence SL1=...,SL2=...`), each `sl` TSS gets the main leader as `sl_leader`. In *C. elegans*, SL2 marks downstream genes of operons.

1. **Genes.** Genes of `--operon_biotypes` (default `protein_coding`) are taken with their ends after UTR annotation. Genes without a `gene_biotype` count as protein-coding when they have a CDS. Genes on the other strand are ignored, so nested or interleaved genes do not interrupt an operon (as `calcOperons()` in annotationpolish).
2. **Grouping.** A gene joins the operon of the gene before it, in transcription order, when the distance from the operon's 3'-most stop codon to the gene's start codon is below `--operon_max_gap` (500) bases (`--operon_gap_from cds`). Overlapping genes join.
   - **Why codons and not transcript ends:** transcript-end distances would be the real intercistronic regions, but annotated UTR ends were too imprecise on the I69 EGAPx and O10 annotations, so they gave less SL-consistent operons. See [rationale.md](rationale.md).
   - **Alternatives:** `--operon_gap_from transcript` (transcript ends below `--operon_max_gap_tx`, 100) and `evidence` (transcript ends when both are supported by PAS or poly(A) reads and CAGE, otherwise codons).
3. **Splitting.** A non-SL (`nosl`) TSS on the same strand, between the end of the previous gene and a gene's start codon, starts a new operon at that gene, however short the intergenic distance: the gene has its own promoter (as `split_operons()`). `--operon_split_min_count` sets the minimum count of such a TSS. A non-SL TSS within the SL acceptor window (5 bp either side) of an SL TSS does not split (SL acceptor rule, above).
4. **Size.** Runs of at least two genes are reported.

In true operons the downstream genes are trans-spliced. The MultiQC table therefore reports:
- the share of downstream genes with an SL TSS between them and the previous gene;
- the same share for the genes that were split off.

## Recommended settings by organism

The defaults are meant for metazoans with or without trans-splicing. For tunicates (*Oikopleura*, *Ciona*), trans-splicing organisms without operons, nematodes, and CAGE without SL information, see [rationale.md, section 0](rationale.md#0-recommended-settings-by-organism).

## Parameters

All parameters, with defaults, are in the [parameter schema](../nextflow_schema.json) (`nextflow run luscombeu/utrannotation --help`). The main ones:

| Parameter               | Default       | Meaning                                                              |
| ----------------------- | ------------- | -------------------------------------------------------------------- |
| `--cluster_set`         | `broad`       | cagescan cluster set (`broad` or `sharp`)                            |
| `--tss_max_dist`        | `500`         | maximum 5' UTR length (non-coding: maximum 5' end move), along the transcript |
| `--tss_selection`       | `auto`        | `auto` (strongest mature-mRNA 5' end), `closest_atg` or `dominant`  |
| `--tss_min_rel`         | `0.1`         | minimum count relative to the strongest TSS (see `--tss_rel_whole_utr`) |
| `--tss_rel_whole_utr`   | `true`        | reference for `--tss_min_rel`: the whole EGAPx 5' UTR, not only the window |
| `--tss_noncoding`       | `true`        | also move the 5' end of lncRNAs and misc_RNAs                        |
| `--tss_min_count`       | `0`           | minimum pooled TSS count                                             |
| `--utr5_min_len`        | `0`           | minimum 5' UTR length (0: a TSS at the ATG is allowed)               |
| `--sl_min_fraction`     | `0.5`         | minimum SL read fraction for `tss_type=sl`                           |
| `--keep_polya_evidence` | `true`        | keep EGAPx 3' ends supported by poly(A) reads                        |
| `--pas_motifs`          | `AATAAA,ATTAAA` | comma-separated polyadenylation-signal motifs                      |
| `--pas_window`          | `500`         | maximum distance of the signal from the EGAPx 3' end                 |
| `--pas_offset`          | `15`          | bases from the motif end to the 3' end                               |
| `--trans_splicing`      | `true`        | label TSSs `sl` / `nosl` from the cagescan SL counts                 |
| `--operons`             | `true`        | annotate operons (genomes with SL data only)                         |
| `--operon_gap_from`     | `cds`         | measure operon gaps between codons (`cds`), transcript ends, or `evidence` |
| `--operon_max_gap`      | `500`         | stop-to-start codon distance below which genes join an operon        |
| `--operon_max_gap_tx`   | `100`         | transcript-end distance limit (`transcript` / `evidence` modes)      |
| `--operon_split_min_count` | `0`        | minimum count of a non-SL TSS that splits an operon                  |
| `--sl_window_up`        | `5`           | SL acceptor rule: ignore non-SL TSSs this close upstream of an SL acceptor |
| `--sl_window_down`      | `5`           | the same downstream (`-1`: the whole 5' UTR)                         |
| `--operon_biotypes`     | `protein_coding` | gene biotypes grouped into operons (`all` for every gene)        |

## References

- Conrad R, Thomas J, Spieth J, Blumenthal T (1991) Insertion of part of an intron into the 5' untranslated region of a *Caenorhabditis elegans* gene converts it into a trans-spliced gene. Mol Cell Biol 11:1921–1926. [doi:10.1128/mcb.11.4.1921-1926.1991](https://doi.org/10.1128/mcb.11.4.1921-1926.1991)
- Conrad R, Lea K, Blumenthal T (1995) SL1 trans-splicing specified by AU-rich synthetic RNA inserted at the 5' end of *Caenorhabditis elegans* pre-mRNA. RNA 1:164–170. PMID 7585246
- Saito TL et al. (2013) The transcription start site landscape of *C. elegans*. Genome Res 23:1348–1361. [doi:10.1101/gr.151571.112](https://doi.org/10.1101/gr.151571.112)
- Fejes-Toth K et al. (2009) Post-transcriptional processing generates a diversity of 5'-modified long and short RNAs. Nature 457:1028–1032. [doi:10.1038/nature07759](https://doi.org/10.1038/nature07759)
- Ganot P et al. (2004) Spliced-leader RNA trans splicing in a chordate, *Oikopleura dioica*, with a compact genome. Mol Cell Biol 24:7795–7805. [doi:10.1128/MCB.24.17.7795-7805.2004](https://doi.org/10.1128/MCB.24.17.7795-7805.2004)
- Plessy C et al. (2024) Extreme genome scrambling in marine planktonic *Oikopleura dioica* cryptic species. Genome Res 34:426–440. [doi:10.1101/gr.278295.123](https://doi.org/10.1101/gr.278295.123)
- Matsumoto J et al. (2010) High-throughput sequence analysis of *Ciona intestinalis* SL trans-spliced mRNAs: alternative expression modes and gene function correlates. Genome Res 20:636–645. [doi:10.1101/gr.100271.109](https://doi.org/10.1101/gr.100271.109)
- Yokomori R et al. (2016) Genome-wide identification and characterization of transcription start sites and promoters in the tunicate *Ciona intestinalis*. Genome Res 26:140–150. [doi:10.1101/gr.184648.114](https://doi.org/10.1101/gr.184648.114)
- Yokomori R, Kusakabe TG, Nakai K (2024) Characterization of trans-spliced chimeric RNAs: insights into the mechanism of trans-splicing. NAR Genom Bioinform 6:lqae067. [doi:10.1093/nargab/lqae067](https://doi.org/10.1093/nargab/lqae067)
- Danks GB, Navratilova P, Lenhard B, Thompson EM (2018) Distinct core promoter codes drive transcription initiation at key developmental transitions in a marine chordate. BMC Genomics 19:164. [doi:10.1186/s12864-018-4504-5](https://doi.org/10.1186/s12864-018-4504-5)
- Lasda EL, Blumenthal T (2011) Trans-splicing. WIREs RNA 2:417–434. [doi:10.1002/wrna.71](https://doi.org/10.1002/wrna.71)
- Blumenthal T (2012) Trans-splicing and operons in *C. elegans*. WormBook. [doi:10.1895/wormbook.1.5.2](https://doi.org/10.1895/wormbook.1.5.2)

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run luscombeu/utrannotation --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run luscombeu/utrannotation -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull luscombeu/utrannotation
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [luscombeu/utrannotation releases page](https://github.com/luscombeu/utrannotation/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
