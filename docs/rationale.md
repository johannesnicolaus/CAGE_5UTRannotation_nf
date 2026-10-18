# luscombeu/utrannotation: Rules, parameters and rationale

This file records every rule the pipeline applies, its parameters and defaults, why it was chosen, what evidence supports it, and what was tried and rejected. Keep it up to date whenever a rule or default changes, and add a line to the [decision log](#decision-log).

Kinds of evidence, as tagged below:

- **[lit]**: published literature. References are at the end, checked against PubMed.
- **[data]**: our own measurements, on two genomes. Parameters are checked on both.
  - **O10:** the Osaka CAGE run (cagescan `osa_MGI_220209_v4` on the O10 genome; 6 libraries: D2_4, DE1, DE2, EB1, EB2, EB3) with the O10 v1.3 annotation.
  - **I69:** the Okinawa CAGEr consensus clusters (oist/LuscombeU-CAGE_libraries, SL and non-SL clustered separately) with the **I69 EGAPx annotation** (`egapx_runs/annotation_egapx/I69_prim`).
- **[prev]**: carried over from the unit's earlier code. These are `annotationpolish` (`add_five_prime_utr()`, `add_three_prime_utr()`, `calcOperons()`, `split_operons()`) and `LuscombeU_5primeanno`.
- **[choice]**: a decision by the pipeline author (Johannes Nicolaus Wibisana), with the reason given.

## Contents

0. [Recommended settings by organism](#0-recommended-settings-by-organism)
1. [Inputs and conventions](#1-inputs-and-conventions)
2. [SL labels of TSS clusters](#2-sl-labels-of-tss-clusters)
3. [5' ends of coding mRNAs](#3-5-ends-of-coding-mrnas)
4. [SL acceptor rule](#4-sl-acceptor-rule)
5. [5' ends of non-coding transcripts](#5-5-ends-of-non-coding-transcripts)
6. [3' ends](#6-3-ends)
7. [Trans-splicing and operon switches](#7-trans-splicing-and-operon-switches)
8. [Operons](#8-operons)
9. [Input checks](#9-input-checks)
10. [Validation so far](#10-validation-so-far)
11. [Known limitations and open questions](#11-known-limitations-and-open-questions)
12. [Decision log](#decision-log)
13. [References](#references)

---

## 0. Recommended settings by organism

**The defaults are meant for metazoans in general,** with or without trans-splicing. Change parameters only for the reasons below.

| Organism group | Settings | Why |
| --- | --- | --- |
| No SL trans-splicing (most metazoans) | defaults | Without SL clusters every TSS is `nosl`, `--tss_selection auto` is the dominant TSS (the CAGE convention), and no operons are annotated. |
| Tunicates (*Oikopleura*, *Ciona*) | defaults | Settings for TSSs near SL acceptors (`--sl_window_up/down 5`) and for operons (`--operon_gap_from cds`, `--operon_max_gap 500`) were chosen and tested on *O. dioica* (I69 EGAPx, O10, Ganot operon). Tunicates have genuine non-trans-spliced TSSs tens of nt from acceptors (sections 3 and 4), so do not widen the SL windows. |
| *Ciona* specifically | defaults; `--operons` can stay on | *Ciona* has SL trans-splicing for about half of expressed genes and some polycistronic units, most of them without intercistronic sequence [lit: Satou 2006]. The pipeline has not been run on *Ciona*. |
| Trans-splicing organisms without operons (e.g. some flatworms, cnidarians) | `--operons false` | Trans-splicing does not imply operons [lit: Lasda & Blumenthal 2011]. |
| Nematodes (*C. elegans*) | defaults; see [Notes for nematodes](#notes-for-nematodes) | Benchmarked against WormBase (section 10). 5' ends and SL labels reproduce WormBase closely. Operons are only approximate, and 3' ends are weak when the input UTRs are over-extended. |
| Organisms whose CAGE has no SL information | `--trans_splicing false` | Every TSS is `nosl`; no acceptor rows and no operons. |

### Notes for nematodes

Tested on *C. elegans* (2026-10-18) with:
- SLIC-CAGE of total RNA (SRR10568191, PRJNA592867; young-adult hermaphrodites, 15.3 M read pairs), run through luscombeu/cagescan;
- the WormBase WS298 gene set in NCBI RefSeq format (GCF_000002985.6) as the reference.

Results are in section 10. Practical notes:

1. **Leaders.** Pass SL1 and the SL2 family to cagescan as several leaders (`--sl_sequence SL1=…,SL2a=…`). In SRR10568191, 52% of READ1s start with SL1 (`GGTTTAATTACCCAAGTTTGAG`) and about 6% with ten SL2 variants. The variants used are listed in `utrannotation_runs/celegans/cagescan_run/params.yaml`. Read the main leader of each TSS as `sl_leader`; SL2 marks downstream operon genes [lit: Blumenthal 2012].
2. **Annotation.** The RefSeq/WormBase GFF3 works as input, since it has the same NCBI format as EGAPx. About 10,000 RefSeq mRNAs are flagged `partial` (no known UTR end). The pipeline leaves those ends alone.
3. **5' ends.** The defaults work. SL acceptors land on WormBase's SL acceptor sites (99.1% exact). `--sl_window_up 50` (the nematode minimum outron; Conrad 1991, 1995) is an option, but was not needed in this test.
4. **3' ends.** The PAS rule is reliable only when the input 3' ends are roughly right, or absent. If they are over-extended, it tends to pick chance motifs near the wrong end (section 10). In *C. elegans*, 13% of polyadenylation sites have no detectable PAS, and trans-spliced mRNAs often have none or a variant [lit: Mangone 2010]. Two ways to handle this:
   - `--pas_motifs` could add AATGAA, the most common *C. elegans* variant [lit: Tian & Graber 2012];
   - better, use poly(A)-site data where it exists, since WormBase 3' ends come from 3P-Seq [lit: Jan 2011].
5. **Operons.** Treat `--operons` as approximate. Against WormBase operons (CEOP, WBPS18), adjacent-gene-pair precision is 0.35 and recall 0.38 with the defaults. No setting tested exceeded precision 0.42 (table in section 10). Hybrid operons (internal promoters) are split by our rule, but WormBase keeps them [lit: Blumenthal 2015].
6. **cagescan with NCBI GTFs.** cagescan's last step, `GFFCOMPARE`, failed on the RefSeq GTF ("no valid ID found for GFF record"). That step only compares cagescan's own transcripts. The TSS cluster outputs used here were complete.

## 1. Inputs and conventions

| Item | Rule | Rationale |
| --- | --- | --- |
| Annotation | GFF3 gene → transcript → exon/CDS (EGAPx). Any transcript with a CDS is coding, which also covers AUGUSTUS-style `transcript` rows. | EGAPx is the target; generalised after the O10 v1.3 test. |
| CDS | Includes the stop codon, as EGAPx writes it. The 3' UTR starts after the last CDS base. | EGAPx convention. |
| Partial models (`start_range` / `end_range`) | Kept unchanged at the open end. | The UTR of a truncated CDS is meaningless. |
| UTR rows | Always regenerated from exons and CDS. Input `five_prime_UTR`, `three_prime_UTR` and `UTR` rows are dropped. | One source of truth. |
| Evidence attributes | Kept on the transcript row only. UTR rows carry only `ID` and `Parent`. | [choice] The TSS type describes the transcript, and it must be present even when the 5' UTR has length 0. |
| CDS | Never changed. Each script asserts this. | Protein models are EGAPx's responsibility. |
| CAGE input | Samplesheet column `tss` takes either a cagescan results directory (`tss_clusters/<set>/`, `broad` by default via `--cluster_set`) or a BED of TSS clusters from any CAGE workflow; `tss_sl` takes an optional BED of SL acceptor clusters. | [choice, user] Most CAGE is single-end, and cagescan currently needs pairs (single-end mode proposed in johannesnicolaus/cagescan_nextflow#5). The BED route was checked on I69: the CAGEr BEDs used directly gave identical results to the same clusters repackaged as cagescan output. |
| BED peak | `thickStart` (BED ≥ 8 columns), otherwise the 1-bp cluster, otherwise the cluster's 5' end (with a warning). The score column is the count. | CAGEr, CAGEfightR and cagescan write the dominant CTSS in `thickStart`. |

## 2. SL labels of TSS clusters

With BED input, clusters in `tss` are `nosl` and clusters in `tss_sl` are `sl`; the threshold below does not apply.

| Parameter | Default | Rule |
| --- | --- | --- |
| `--sl_min_fraction` | 0.5 | A cluster is `sl` when SL READ1s / all READ1s (summed over samples) ≥ 0.5. Otherwise it is `nosl`. |

- **Separate clusters.** With cagescan's separate SL clusters (its default since `--sl_separate_clusters`), `SL…` clusters contain only SL reads and `TC…` clusters none, so the fractions are exactly 1 and 0. The threshold matters only for runs with pooled clusters. [data]
- **No SL table.** Every cluster is `nosl` when cagescan wrote no SL table, i.e. the organism has no trans-splicing or `--sl_sequence` was not given.
- **Several leaders.** When cagescan searched for several leaders, the main leader of each SL cluster is written as `sl_leader`. In *C. elegans*, SL2 marks downstream operon genes [lit: Blumenthal 2012].

## 3. 5' ends of coding mRNAs

| Parameter | Default | Rule and rationale |
| --- | --- | --- |
| `--tss_max_dist` | 500 | Candidate TSSs give a 5' UTR of at most 500 nt, measured along the transcript (5' UTR exons, then the genome upstream of the EGAPx 5' end). [prev] `add_five_prime_utr(dist = 500)`. [lit] Average 5' UTR lengths are 100–220 nt across taxa: invertebrates 222 nt, human 210 nt [Mignone 2002]. 500 covers typical UTRs, and longer ones keep their EGAPx end. |
| `--tss_selection` | `auto` | The strongest **mature-mRNA 5' end**: SL acceptors, and `nosl` TSSs that are not upstream of a qualifying acceptor, compete by count. Without SL data this is the dominant TSS. See below. |
| `--tss_min_rel` | 0.1 | Candidates below 10% of the strongest reference peak are dropped. CAGE "exon painting" leaves weak capped 5' ends inside exons [lit: Fejes-Toth 2009], and these would otherwise win "closest to ATG". |
| `--tss_rel_whole_utr` | true | The reference for `--tss_min_rel` is every peak in the EGAPx 5' UTR, including those beyond `--tss_max_dist`. [choice] Otherwise a weak peak near the ATG could trim a UTR whose strong TSS lies further away. Such mRNAs keep the EGAPx end. |
| `--tss_min_count` | 0 | No absolute threshold. cagescan clusters already require ≥ 5 TPM in a sample. |
| `--utr5_min_len` | 0 | A TSS at the ATG is allowed. SL acceptors can sit at or next to the ATG: in *O. dioica* the leader-to-start distance averages 22 nt, and is ≤ 5 nt for 20% of genes [lit: Ganot 2004]. With the old default of 1, such genes fell back to EGAPx and lost their `sl` label. |

**Why `auto`.** [data, lit]

- **`dominant` alone is wrong for trans-spliced genes.** A `nosl` peak upstream of an SL acceptor is often the outron start. Its reads come from capped pre-mRNA, which is removed and degraded by trans-splicing [lit: Morton & Blumenthal 2011]; Saito et al. 2013 needed nuclear enrichment to see outron TSSs at all. So it is not the 5' end of the mature mRNA. On O10, `dominant` moved 141 transcripts from an SL acceptor to an upstream `nosl` peak.
- **`closest_atg` alone is wrong for genes without trans-splicing.** It picks a minor TSS nearer the ATG over the dominant promoter.
- **Why `nosl` TSSs downstream of an acceptor still compete.** A stronger one can be the main 5' end of a non-trans-spliced mRNA. Examples are *Ciona* "infrequently trans-spliced" genes [lit: Matsumoto 2010] and *O. dioica* testis promoters [lit: Danks 2018].
- **Effect on O10.**
  - `auto`: 4,951 mRNAs `sl` (98.7% with AG).
  - `closest_atg`: 4,859. `dominant`: 4,812.
  - Under `auto`, 80 mRNAs choose a `nosl` TSS that is downstream of, and stronger than, their acceptor.
- **The other kind is always recorded.** When the chosen TSS is `nosl`, the strongest qualifying SL acceptor goes into `sl_acceptor_peak` / `sl_acceptor_count`. When it is `sl`, the strongest `nosl` TSS goes into `nosl_tss_peak` / `nosl_tss_count`; that is an outron start or a non-trans-spliced TSS.
- **History.** `closest_atg` was the first default (user preference for *O. dioica*) and remains available.

Other rules:

- **Not in 5' UTR introns.** Candidates in 5' UTR introns are excluded. A TSS inside an intron would need a different exon structure.
- **CDS guard.** The new 5' end must not overlap another gene's CDS on the same strand.
- **Ties.** Ties prefer `nosl`, then the higher count. [prev] `add_five_prime_utr()` kept non-SL peaks at the same coordinate.
- **No candidate.** If no candidate qualifies, the EGAPx 5' end is kept and `tss_evidence=EGAPx` is written.
- **Edit.** Exons upstream of the TSS are dropped, so the UTR is extended or trimmed. [choice] The earlier code only extended; EGAPx UTRs are often over-extended (I25 mean 611 nt, median 152 nt).
- **SL acceptor check.** For `sl` TSSs, `acceptor_AG` tests that the two genomic bases before the TSS are AG, and a 1-bp `trans_splice_acceptor_site` row is written (SO:0000706). [data] 98.7% of `sl` TSSs on O10 follow an AG.

## 4. SL acceptor rule

| Parameter | Default | Rule |
| --- | --- | --- |
| `--sl_window_up` | 5 | `nosl` peaks up to 5 nt upstream of a qualifying SL acceptor are ignored. `-1` means no limit. |
| `--sl_window_down` | 5 | The same for `nosl` peaks up to 5 nt downstream. `-1` means the whole 5' UTR. |

A "qualifying acceptor" is an `sl` peak that passed the count filters of section 3. The rule applies to the TSS choice (sections 3 and 5) and to operon splits (section 8).

The strongest TSS of the other kind is always recorded on the transcript (section 3): `sl_acceptor_peak` / `sl_acceptor_count`, or `nosl_tss_peak` / `nosl_tss_count`. That keeps dual trans-spliced and non-trans-spliced genes visible.

**Why the near peaks are ignored.** A `nosl` peak within a few nt of an acceptor is SL-derived. [data]

- **The −2 spike.** 142 `TC` clusters lie within 5 nt of an `SL` cluster. 104 of them peak exactly 2 nt upstream, and in 103 of those the two bases in between are the acceptor AG. The `nosl` peak is on the A of the AG.
- **The reads carry no leader.** At two such sites (PAR:7573452+, XSR:4851460−), the reads align full length with no soft clip and start with `AG` followed by the transcript. They contain no leader except its final AG; the *O. dioica* leader ends …AACAG.
- **The effect is library-specific.**
  - In DE1 and DE2 these reads outnumber the SL reads at the same sites. In EB1–3 they are absent, and D2_4 has few.
  - DE1 and DE2 also have the lowest share of recognised SL READ1s: 19.7% and 23.4%, against 34.8–43.8% in the other libraries.
  - This is a library artefact, not a promoter. The molecular cause is unknown.
- **Upstream (−3 to −5 nt): probably the same, weaker evidence.** 13 peaks. The bases in between still end in the acceptor AG, and the reads are enriched in DE1 and DE2 (60% of reads, against 23% for all `TC` clusters).
- **Downstream (+1 to +5 nt): unexplained.** 24 peaks. They are not preceded by an AG and come from the DE libraries only slightly more than usual.

**Why the window is only 5 nt.** Further from an acceptor, `nosl` TSSs can be real, especially in tunicates. [lit]

- **Upstream, in *Ciona*.**
  - Matsumoto et al. 2010, Fig. 5, pattern B: 38 of 107 genes that make both trans-spliced and non-trans-spliced RNAs have the non-trans-spliced TSS in the same exon, just upstream of the acceptor. The average distance is 36 bp. They argue for a separate promoter.
  - Yokomori et al. 2016 found TSS clusters overlapping acceptors and suggest "regions near TASs can function as alternative promoters".
- **Downstream, in *O. dioica*.** Danks et al. 2018 describe 89 male-specific promoters whose TSS lies about 40–50 nt downstream of a trans-splice site, giving mRNAs without the acceptor.
- **Why the published *Ciona* data don't settle it.** The TSS-seq studies (Yokomori 2016, 2024) excluded outrons shorter than 51 nt by design, borrowing Conrad et al. 1995. The absence of shorter cases there is a filter, not a finding.

**Outron length, for reference.** The outron is the pre-mRNA stretch from the TSS to the trans-splice acceptor, which trans-splicing removes [lit: Conrad 1991].

- **Shortest shown to work:** synthetic AU-rich outrons of 51 nt or more are trans-spliced efficiently [lit: Conrad 1995, abstract].
- **Too short:** an acceptor with only 30 nt upstream was not used [lit: Conrad 1991].
- **Natural *C. elegans* outrons:** 10th percentile 59 nt, median 369 nt [lit: Saito 2013].
- **In *Ciona*:** median 458 nt [lit: Yokomori 2024].
- **Reported by the user, not verified:** shorter synthetic outrons (e.g. 41 nt) fail. The Conrad 1995 full text was not available to check.

So a `nosl` TSS a few tens of nt upstream of an acceptor cannot start that acceptor's outron. It can still be the TSS of a separate non-trans-spliced RNA, which is why the rule does not use 50 nt.

**What the rule changes on O10** [data]:

| Version | mRNAs labelled `sl` (with AG) | Operon splits | Split-off genes with an SL TSS |
| --- | --- | --- | --- |
| No rule | 4,848 (4,786) | 1,568 | 20.3% |
| 5 nt, SL stronger (first version) | 4,861 (4,799) | 1,538 | 18.7% |
| Downstream unlimited + 50 nt upstream | 5,031 (4,967) | 1,349 | 7.3% |
| **5 nt up / 5 nt down (current window)** | **4,859 (4,796)** | **1,535** | **18.6%** |

These rows used `--tss_selection closest_atg`. With the current default `auto`, O10 has 4,951 `sl` mRNAs (4,888 with AG).

With the current rule, 172 `nosl` mRNAs record an SL acceptor (`sl_acceptor_peak`). These are mostly cases with a `nosl` peak more than 5 nt downstream of an acceptor. Under the wider rule they would have been labelled `sl`.

## 5. 5' ends of non-coding transcripts

| Parameter | Default | Rule |
| --- | --- | --- |
| `--tss_noncoding` | true | lncRNA, EGAPx misc_RNA (`transcript` without CDS) and ncRNA rows get a CAGE 5' end. [choice, user] |

- **Candidates.** Peaks within `--tss_max_dist` of the EGAPx 5' end. Upstream peaks are measured along the genome. Downstream peaks are measured along the exons and must leave at least 100 nt of transcript (`--nc-min-length` in the script).
- **Choice.** The highest count wins, with ties going to the smallest move. These transcripts have no ATG, so "closest to ATG" does not apply.
- **Other rules.** The CDS guard, `--tss_min_count` and the SL acceptor rule all apply. 3' ends are not changed.

## 6. 3' ends

The first rule that applies sets the 3' end. [choice, user, chosen over "first PAS after the stop codon"]

1. **`polyA_reads`** (`--keep_polya_evidence`, default true). EGAPx flagged the model `experiment=COORDINATES: polyA evidence`, meaning RNA-seq reads that contain the poly(A) tail map the cleavage site. This is direct evidence, so the EGAPx end is kept.
2. **`PAS`.** Short-read coverage locates the end and the motif sets the base. The candidate motifs are `--pas_motifs` (`AATAAA`) on the same strand that:
   - lie entirely after the stop codon;
   - are within `--pas_window` (500 nt) of the EGAPx end, measured along the transcript;
   - are not in a 3' UTR intron;
   - do not reach another gene's CDS on the same strand;
   - do not extend into another gene's transcript on the same strand.

   The motif **closest to the EGAPx end** is used. The new end is the motif end + `--pas_offset` (15). [prev] The window comes from `add_three_prime_utr(dist = 500)`.
3. **`EGAPx`.** Otherwise the EGAPx end is kept.

Rationale:

- **Why closest to the EGAPx end, not first after the stop codon.** Hexamer motifs occur by chance: about 0.6 AATAAA per kb per strand at 41% GC. "First after stop" would truncate real UTRs at chance motifs.
  - **Correction (2026-10-18):** an earlier version said *O. dioica* is "AT-rich". O10 is 41.4% GC, like human, so it does not have more chance AATAAA than other genomes.
- **`--pas_motifs AATAAA,ATTAAA`.** [lit: Tian & Graber 2012] In human polyA sites, AAUAAA is found at ~53–58% and AUUAAA at ~15–17%. AUUAAA has ~80% of AAUAAA's processing efficiency. In *C. elegans* only ~40% of sites have AAUAAA, and AAUGAA is the most common variant. The old code used AATAAA only.
  - On O10, adding ATTAAA raises motif sites from 136k to 216k, and PAS-set 3' ends from 64% to 71%.
  - There are no tunicate-specific data. In the *O. dioica* cyclin D3 operon, only the last cistron has AAUAAA [lit: Ganot 2004]. AATAAA alone is a conservative alternative.
- **`--pas_offset 15`.** The signal lies "within ~40 nt upstream" of the cleavage site in mammals [lit: Tian & Graber 2012], and cleavage is commonly 10–30 nt downstream of it. In *O. dioica* cyclin D, AAUAAA sits 17 nt upstream of the poly(A) tail [lit: Ganot 2004]. 0 (the old behaviour) put the end at the motif, which is always too short.
- **The next-transcript guard.** In operons the 3' end stops before the next gene's 5' end, which the 5' step has already set. Overlaps already present in EGAPx are not counted.

## 7. Trans-splicing and operon switches

| Parameter | Default | Rule |
| --- | --- | --- |
| `--trans_splicing` | true | Use the SL counts when cagescan wrote them. If false, every TSS is `nosl` and there are no acceptor rows and no operons. |
| `--operons` | true | Annotate operons. This runs only for genomes with SL data. |

Rationale [lit]:

- **Operons depend on trans-splicing.** The leader caps the downstream genes of a polycistronic pre-mRNA [Blumenthal 2004; Morton & Blumenthal 2011].
- **Trans-splicing does not imply operons.** "Operon resolution can only account for a small proportion of the total amount of SL trans-splicing" [Lasda & Blumenthal 2011]. Many trans-spliced genes are monocistronic (outron replacement).

So there are two independent switches:

- trans-splicing with operons (*Oikopleura*, *C. elegans*): use the defaults;
- trans-splicing without operons: `--operons false`;
- no trans-splicing: cagescan writes no SL counts, so everything is `nosl` and there are no operons.

## 8. Operons

| Parameter | Default | Rule |
| --- | --- | --- |
| `--operon_gap_from` | `cds` | Intergenic distance from the stop codon of the previous gene to the start codon of the next. Alternatives: `transcript` (transcript ends) and `evidence` (transcript ends when both are evidence-supported, otherwise CDS). |
| `--operon_max_gap` | 500 | CDS distance below which same-strand genes join (at most 499 bases between the codons). Genes on the other strand are ignored, and overlapping genes join. |
| `--operon_max_gap_tx` | 100 | Transcript-end distance limit for `transcript` / `evidence`. |
| `--operon_split_min_count` | 0 | Minimum count of a `nosl` TSS that splits an operon. |
| `--operon_biotypes` | `protein_coding` | Genes considered (`gene_biotype`, or "has a CDS" when the attribute is absent). |
| `--sl_window_up` / `--sl_window_down` | 5 / 5 | `nosl` TSSs this close to an SL TSS do not split (section 4). |

- **Grouping.** [lit, prev] Plessy et al. 2024 define operons as same-strand genes "separated by an intergenic distance of at most 500 bp". In their `calcOperons(window = 500)`, genes merge when start(B) − end(A) ≤ 500. Opposite-strand genes do not interrupt, and an operon needs at least 2 genes.
- **Why CDS distances, not transcript ends.** [data, choice] Plessy's 500 bp was calibrated on models without UTRs. The user asked for distances between transcript ends instead, now that UTRs are annotated: those are the true intercistronic regions, 23–30 nt in *O. dioica* [lit: Ganot 2004] and about 100 bp in *C. elegans* [lit: Blumenthal 2015]. Tested on I69 EGAPx and O10.
  - **The short-gap signal is real.** On I69 EGAPx, the share of CAGE-supported downstream genes that start at an SL acceptor is 96–97% at 0–60 nt from the previous transcript, 86% at 61–100 nt, 65% at 101–200 nt and 27% beyond 1 kb.
  - **But annotated UTR ends are too imprecise to use as gaps.** At the Ganot locus in I69 EGAPx, transcript-end gaps are 23–362 nt against true intercistronic regions of 23–30 nt. XPO7's 3' end was kept from EGAPx and is 285 nt short, and TUBD1 and CYCD3 have no CAGE TSS, so they keep EGAPx 5' ends.
  - **Transcript-end gaps of 60 nt break the Ganot operon** on I69 EGAPx and give 14–35% of genes in operons.
  - **The evidence-aware mode doesn't help.** It uses transcript ends only when the upstream 3' end comes from PAS or poly(A) reads and the downstream 5' end from CAGE. Both quality measures got worse:

    | Genome | Rule | downstream genes with SL TSS | split-off genes with SL TSS |
    | --- | --- | --- | --- |
    | I69 | CDS 500 | **51.9%** | **31.0%** |
    | I69 | evidence 100 | 47.7% | 37.2% |
    | O10 | CDS 500 | **57.9%** | **21.6%** |
    | O10 | evidence 100 | 47.8% | 31.5% |

  - **CDS windows of 200–750 nt on both genomes:** the downstream-SL share is flat (~52% I69, ~58% O10). The split-off SL share is lowest at 400–500. 500 is the smallest window that keeps the Ganot operon whole: its 403-bp CDS gap breaks it at 400, as Plessy's vignette noted. Larger windows add genes without improving either measure.
- **Splitting.** [prev, choice, user] A `nosl` TSS on the same strand between the end of the previous gene and a gene's start codon means the gene has its own promoter, so a new operon starts there, however short the gap. This is `split_operons()`. The SL acceptor rule (section 4) applies.
- **Benchmark.** [data, lit] On OKI2018_I69 v2.3 with the Okinawa CAGE consensus clusters, the Ganot et al. 2004 cyclin D3 operon (AJ628168) is recovered as one minus-strand operon at PAR:15956186–15967604: POC5 → IQCK → XPO7 (RBP16) → DYNLL2 → TUBD1 → MBF1 → CYCDC.
  - **Matches the literature:** Ganot's five genes are all there, in order, with the POC5 leader and IQCK that Plessy's analysis identified.
  - **Gaps:** 22–27 bp within the operon, against Ganot's 23–30 nt intercistronic regions.
  - **Neighbours:** PDCL3 upstream is split off by a `nosl` TSS of 783 reads.
  - **Robustness:** the result is the same for max gaps of 100–750 bp. Without the split it grows to 10 genes at 500 bp.
  - **I69 EGAPx** (CDS 500): one operon of 8 genes, POC5 → an EGAPx "uncharacterized protein" model overlapping IQCK → IQCK → XPO7 → DYNLL2 → TUBD1 → MBF1 (EDF1 homolog) → CYCD3. Its CDS gaps reproduce Plessy's (362, 65, 148, 403, 122).
- **Genome-wide.**
  - OKI: about 2,973 operons with 46.6% of genes, against 3,124 operons and about 49% in Plessy's vignette.
  - O10 (current defaults): 2,395 operons with 39.5% of genes.
  - I69 EGAPx (current defaults): 2,455 operons with 42.4% of genes.
  - Against the v1.3 operon set: 84% of operonic genes shared, and 1,504 operons with identical gene lists (measured with the earlier rule set).

## 9. Input checks

- **Annotation sequences must exist in the genome FASTA.**
- **The pipeline stops if no TSS lies on a genome sequence, or if any TSS lies beyond the end of its sequence.** Both mean cagescan was run on another assembly.
- **A missing cagescan clusters or counts file stops the run.** A missing SL counts file means "no trans-splicing data".

## 10. Validation so far

- **2026-10-18, O10 (v1.3 GFF3 + cagescan `osa_MGI_220209_v4`).** Full pipeline run: 12 tasks, 0 failed. Run directory `/flash/LuscombeU/nico/utrannotation_runs/osa_O10_v1.3`.
  - **CDSs:** all unchanged. 17,646 of 17,686 proteins are identical to the v1.3 peptides; the other 40 differ only in how internal stops are written (`*` vs `.`).
  - **Later rule changes** (SL acceptor rule and window) were measured by re-running the 5' and operon scripts on that run's outputs, not by a new full run.
- **OKI v2.3 operon benchmark:** see section 8.
- **2026-10-18, I69 EGAPx + Okinawa CAGEr consensus clusters.** Full runs completed. Run directories are `/flash/LuscombeU/nico/utrannotation_runs/oki_I69_egapx` (clusters repackaged as cagescan output) and `oki_I69_egapx_bed` (the same CAGEr BEDs given directly as `tss` / `tss_sl`). Both give identical results.
  - **5' ends:** 50.7% of 22,265 mRNAs get a CAGE TSS. 53.9% of those are `sl`, 97.3% of them with AG. 7,098 5' UTRs are trimmed and 3,999 extended: median 156 → 71 nt, mean 607 → 371 nt.
  - **3' ends:** poly(A) reads 21.5%, PAS 53.5%, EGAPx 25.0%.
  - **Non-coding:** 36.6% of 1,437 non-coding transcripts get a CAGE 5' end.
- **2026-10-18, O10 rerun with the current defaults.** Completed. 65.0% of mRNAs get a CAGE TSS, 43.0% of those `sl` (98.7% AG). PAS ends 71.3%.
- **2026-10-18, *C. elegans* accuracy benchmark** (`utrannotation_runs/celegans`; scripts in `tests/benchmark/`).
  - **Setup.**
    - Reference: WormBase WS298 models via RefSeq GCF_000002985.6, scoring only ends not flagged `partial`. That gives 18,640 known 5' ends and 20,785 known 3' ends.
    - Inputs built from the reference: "strip" (no UTRs) and "extend" (both ends pushed out by 0–500 nt).
    - CAGE: SLIC-CAGE total RNA through cagescan, giving 8,413 SL clusters (7,073 mainly SL1, 1,340 mainly SL2) and 15,805 non-SL clusters. 56.7% of read pairs carried a leader.
  - **5' ends** (output; strip / extend):

    | Class | n | exact | ≤ 5 nt | ≤ 20 nt | median error |
    | --- | --- | --- | --- | --- | --- |
    | CAGE `sl` | 7,833 / 8,039 | 60.7% / 60.7% | 75.7% / 75.9% | 89.9% / 90.4% | 0 nt |
    | CAGE `nosl` | 2,362 / 2,371 | 7.9% / 8.1% | 21.0% / 21.3% | 49.3% / 49.9% | 21 nt |
    | no CAGE TSS (input end kept) | 8,445 / 8,230 | as input | | | |

    WormBase's 5' ends of non-trans-spliced genes are themselves less certain [lit: Saito 2013].
  - **SL labels.**
    - 99.1% of `sl` TSSs (7,763 / 7,833) lie exactly on a WormBase SL acceptor site (WBPS19), and only 3% of `nosl` TSSs have one within 5 nt.
    - Our main leader agrees with WormBase's site types: of about 7,800 SL TSSs at a WormBase site, only 7 have a leader WormBase does not list there.
  - **3' ends** (PAS-set ends; strip / extend):
    - within 5 nt: 44.1% / 20.7%;
    - within 20 nt: 57.2% / 27.6%;
    - median error: 9 / 174 nt.

    The reference end lies a median 15 nt (strip) or 14 nt (extend) after the motif end, which confirms `--pas_offset 15`. With over-extended inputs, "closest motif to the input end" picks chance motifs near the wrong end. This is a known weakness (section 11).
  - **Operons vs WormBase CEOP** (WBPS18, 1,431 operons; adjacent protein-coding pairs):

    | Setting | operons | genes in operons | precision | recall |
    | --- | --- | --- | --- | --- |
    | CDS < 500 + non-SL split (default) | 1,894 | 21.3% | 0.347 | 0.383 |
    | CDS < 500, no split | 2,120 | 24.4% | 0.327 | 0.420 |
    | CDS < 300 + split | 1,073 | 11.6% | 0.394 | 0.228 |
    | transcript ends < 200 + split | 1,053 | 11.2% | 0.419 | 0.233 |
    | evidence-aware, 100 | 1,706 | 18.8% | 0.310 | 0.295 |

    58% of SL2-led genes are downstream genes of our operons, against 5% of SL1-led genes.
- **Not yet tested:**
  - the I25 and O10 EGAPx annotations (I25 lacks I25-aligned CAGE; the O10 EGAPx run was still in progress);
  - *Ciona* (ANISEED KH2012 operons are available in `/bucket/LuscombeU/nico/genomes/ciona_robusta/`, but no *Ciona* CAGE run yet);
  - nf-tests (test data pending; see [development.md](development.md)).

## 11. Known limitations and open questions

- **Tunicate dual-TSS genes.** Genes that make both trans-spliced and non-trans-spliced RNAs get one TSS per transcript model. The other form is recorded as `sl_acceptor_peak` or `nosl_tss_peak`. Making a second transcript model is a possible extension; KY2019 *Ciona* models do this [lit: Satou 2019].
- **Operon gaps ignore UTRs.** CDS distances were more reliable than annotated UTR ends on both genomes tested. With accurate 3' ends (e.g. from poly(A)-site sequencing), `--operon_gap_from transcript` with `--operon_max_gap_tx 100` (or 60 for *O. dioica*) would measure the real intercistronic regions.
- **Hybrid operons.** In *C. elegans* an operon can have an internal promoter and stay an operon [lit: Blumenthal 2012, 2015]. Our rule splits at any `nosl` TSS outside the SL window.
- **The DE1/DE2 artefact.** The cause is unknown. Reported upstream as johannesnicolaus/cagescan_nextflow#4.
- **PAS rules.** The motif and offset are general metazoan values; the offset was confirmed on *C. elegans* (median 15 nt). The rule anchors on the input 3' end. When that end is badly over-extended, it picks chance motifs near it: on *C. elegans* with ends pushed out by 0–500 nt, only 21% of PAS-set ends were within 5 nt of the truth, against 44% from inputs without UTRs. A rule that does not depend on the input end, or poly(A)-site data, would be needed to fix this. Not done yet.
- **`tss_evidence=EGAPx`** is written for kept ends even when the input annotation is not from EGAPx.
- **Duplicate isoforms.** Isoforms that become identical after the edits are not merged.

## Decision log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-10-17 | 5' rule: closest to ATG within 500 nt, with a 10% noise filter | User preference; SL acceptor biology; exon painting. |
| 2026-10-17 | 3' rule: hierarchy polyA_reads > PAS closest to the EGAPx end > EGAPx | User choice after discussing "first PAS after stop". |
| 2026-10-18 | Evidence attributes moved to the transcript only; `--utr5_min_len 0` | SL acceptors at the ATG lost their label. |
| 2026-10-18 | `--tss_rel_whole_utr`; 3' guard against the next transcript; non-coding 5' ends; input checks; operons | Review of concerns, user request. |
| 2026-10-18 | `--trans_splicing` and `--operons` as separate switches | Trans-splicing does not imply operons [Lasda & Blumenthal 2011]. |
| 2026-10-18 | SL guard v1: 5 nt, only when the SL peak is stronger | `nosl` peaks at acceptors in the Osaka data. |
| 2026-10-18 | SL guard v2: −2 always, plus 5 nt when stronger | 23 of the −2 artefact peaks were stronger than their SL peak. |
| 2026-10-18 | SL acceptor rule v3: drop `nosl` anywhere downstream, and < 50 nt upstream (`--min_outron`) | User asked for a general rule; minimum outron length [Conrad 1991, 1995; Saito 2013]. |
| 2026-10-18 | SL acceptor rule v4: 5 nt upstream and 5 nt downstream (`--sl_window_up`, `--sl_window_down`); record `sl_acceptor_peak` when a `nosl` TSS wins | Tunicates have genuine non-trans-spliced TSSs tens of nt upstream (*Ciona*, Matsumoto 2010) and downstream (*O. dioica* testis, Danks 2018) of acceptors. The *Ciona* TSS-seq studies' 51-nt minimum was a filter. The user preferred a conservative 5 nt. |
| 2026-10-18 | **Defaults for metazoans in general:** `--tss_selection auto`; `--pas_motifs AATAAA,ATTAAA`; `--pas_offset 15`; `nosl_tss_peak` recorded for `sl` TSSs | User asked for general defaults with notes for tunicates. Evidence: Tian & Graber 2012; Ganot 2004; O10 comparison of the selection modes. |
| 2026-10-18 | **Generic TSS input:** `tss` (cagescan directory or BED) and `tss_sl` replace `cagescan_dir` | Most CAGE is single-end; cagescan needs pairs (issue #5 on cagescan_nextflow). |
| 2026-10-18 | **Operon gaps:** kept CDS distances < 500 (`--operon_gap_from cds`); added `transcript` and `evidence` modes (`--operon_max_gap_tx 100`) | The user asked to use UTR-based gaps. Tested on I69 EGAPx and O10: UTR ends are too imprecise, giving less SL-consistent operons and breaking the Ganot operon. 500 is the smallest CDS window that keeps the Ganot operon whole. |
| 2026-10-18 | Correction: *O. dioica* is not AT-rich (O10 41.4% GC) | The earlier rationale for the PAS rule wrongly said it was. |
| 2026-10-18 | *C. elegans* accuracy benchmark; notes for nematodes; no default changed | 5' ends and SL labels reproduce WormBase closely. `--pas_offset 15` was confirmed. The PAS rule is weak on over-extended inputs, and operons are approximate in nematodes. |

## References

All checked against PubMed unless noted.

- Mangone M et al. (2010) The landscape of *C. elegans* 3'UTRs. Science 329:432–435. [doi:10.1126/science.1191244](https://doi.org/10.1126/science.1191244)
- Jan CH, Friedman RC, Ruby JG, Bartel DP (2011) Formation, regulation and evolution of *Caenorhabditis elegans* 3'UTRs. Nature 469:97–101. [doi:10.1038/nature09616](https://doi.org/10.1038/nature09616)
- Tian B, Graber JH (2012) Signals for pre-mRNA cleavage and polyadenylation. WIREs RNA 3:385–396. [doi:10.1002/wrna.116](https://doi.org/10.1002/wrna.116)
- Mignone F, Gissi C, Liuni S, Pesole G (2002) Untranslated regions of mRNAs. Genome Biol 3:reviews0004. [doi:10.1186/gb-2002-3-3-reviews0004](https://doi.org/10.1186/gb-2002-3-3-reviews0004)
- Blumenthal T (2004) Operons in eukaryotes. Brief Funct Genomic Proteomic 3:199–211. [doi:10.1093/bfgp/3.3.199](https://doi.org/10.1093/bfgp/3.3.199)
- Blumenthal T (2012) Trans-splicing and operons in *C. elegans*. WormBook. [doi:10.1895/wormbook.1.5.2](https://doi.org/10.1895/wormbook.1.5.2)
- Blumenthal T, Davis P, Garrido-Lecca A (2015) Operon and non-operon gene clusters in the *C. elegans* genome. WormBook. [doi:10.1895/wormbook.1.175.1](https://doi.org/10.1895/wormbook.1.175.1)
- Conrad R, Thomas J, Spieth J, Blumenthal T (1991) Insertion of part of an intron into the 5' untranslated region of a *Caenorhabditis elegans* gene converts it into a trans-spliced gene. Mol Cell Biol 11:1921–1926. [doi:10.1128/mcb.11.4.1921-1926.1991](https://doi.org/10.1128/mcb.11.4.1921-1926.1991)
- Conrad R, Lea K, Blumenthal T (1995) SL1 trans-splicing specified by AU-rich synthetic RNA inserted at the 5' end of *Caenorhabditis elegans* pre-mRNA. RNA 1:164–170. PMID 7585246 (abstract only read)
- Danks GB, Navratilova P, Lenhard B, Thompson EM (2018) Distinct core promoter codes drive transcription initiation at key developmental transitions in a marine chordate. BMC Genomics 19:164. [doi:10.1186/s12864-018-4504-5](https://doi.org/10.1186/s12864-018-4504-5)
- Fejes-Toth K et al. (2009) Post-transcriptional processing generates a diversity of 5'-modified long and short RNAs. Nature 457:1028–1032. [doi:10.1038/nature07759](https://doi.org/10.1038/nature07759)
- Ganot P et al. (2004) Spliced-leader RNA trans splicing in a chordate, *Oikopleura dioica*, with a compact genome. Mol Cell Biol 24:7795–7805. [doi:10.1128/MCB.24.17.7795-7805.2004](https://doi.org/10.1128/MCB.24.17.7795-7805.2004)
- Lasda EL, Blumenthal T (2011) Trans-splicing. WIREs RNA 2:417–434. [doi:10.1002/wrna.71](https://doi.org/10.1002/wrna.71)
- Matsumoto J et al. (2010) High-throughput sequence analysis of *Ciona intestinalis* SL trans-spliced mRNAs: alternative expression modes and gene function correlates. Genome Res 20:636–645. [doi:10.1101/gr.100271.109](https://doi.org/10.1101/gr.100271.109)
- Morton JJ, Blumenthal T (2011) RNA processing in *C. elegans*. Methods Cell Biol 106:187–217. [doi:10.1016/B978-0-12-544172-8.00007-4](https://doi.org/10.1016/B978-0-12-544172-8.00007-4)
- Plessy C et al. (2024) Extreme genome scrambling in marine planktonic *Oikopleura dioica* cryptic species. Genome Res 34:426–440. [doi:10.1101/gr.278295.123](https://doi.org/10.1101/gr.278295.123)
- Saito TL et al. (2013) The transcription start site landscape of *C. elegans*. Genome Res 23:1348–1361. [doi:10.1101/gr.151571.112](https://doi.org/10.1101/gr.151571.112)
- Satou Y et al. (2006) Genomic overview of mRNA 5'-leader trans-splicing in the ascidian *Ciona intestinalis*. Nucleic Acids Res 34:3378–3388. [doi:10.1093/nar/gkl418](https://doi.org/10.1093/nar/gkl418)
- Satou Y et al. (2019) A nearly complete genome of *Ciona intestinalis* type A (*C. robusta*) reveals the contribution of inversion to chromosomal evolution in the genus *Ciona*. Genome Biol Evol 11:3144–3157. [doi:10.1093/gbe/evz228](https://doi.org/10.1093/gbe/evz228)
- Yokomori R et al. (2016) Genome-wide identification and characterization of transcription start sites and promoters in the tunicate *Ciona intestinalis*. Genome Res 26:140–150. [doi:10.1101/gr.184648.114](https://doi.org/10.1101/gr.184648.114)
- Yokomori R, Kusakabe TG, Nakai K (2024) Characterization of trans-spliced chimeric RNAs: insights into the mechanism of trans-splicing. NAR Genom Bioinform 6:lqae067. [doi:10.1093/nargab/lqae067](https://doi.org/10.1093/nargab/lqae067)
