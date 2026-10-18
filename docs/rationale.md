# luscombeu/utrannotation: Rules, parameters and rationale

This file records every rule the pipeline applies, its parameters and defaults, why it was chosen, what evidence supports it, and what was tried and rejected. Keep it up to date whenever a rule or default changes, and add a line to the [decision log](#decision-log).

Kinds of evidence, as tagged below:

- **[lit]**: published literature. References are at the end, checked against PubMed.
- **[data]**: our own measurements. Unless stated otherwise they come from the Osaka CAGE run (cagescan `osa_MGI_220209_v4` on the O10 genome; 6 libraries: D2_4, DE1, DE2, EB1, EB2, EB3) and the O10 v1.3 annotation.
- **[prev]**: carried over from the unit's earlier code. These are `annotationpolish` (`add_five_prime_utr()`, `add_three_prime_utr()`, `calcOperons()`, `split_operons()`) and `LuscombeU_5primeanno`.
- **[choice]**: a decision by the pipeline author (Johannes Nicolaus Wibisana), with the reason given.

## Contents

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

## 1. Inputs and conventions

| Item | Rule | Rationale |
| --- | --- | --- |
| Annotation | GFF3 gene → transcript → exon/CDS (EGAPx). Any transcript with a CDS is coding, which also covers AUGUSTUS-style `transcript` rows. | EGAPx is the target; generalised after the O10 v1.3 test. |
| CDS | Includes the stop codon, as EGAPx writes it. The 3' UTR starts after the last CDS base. | EGAPx convention. |
| Partial models (`start_range` / `end_range`) | Kept unchanged at the open end. | The UTR of a truncated CDS is meaningless. |
| UTR rows | Always regenerated from exons and CDS. Input `five_prime_UTR`, `three_prime_UTR` and `UTR` rows are dropped. | One source of truth. |
| Evidence attributes | Kept on the transcript row only. UTR rows carry only `ID` and `Parent`. | [choice] The TSS type describes the transcript, and it must be present even when the 5' UTR has length 0. |
| CDS | Never changed. Each script asserts this. | Protein models are EGAPx's responsibility. |
| CAGE input | cagescan `tss_clusters/<set>/`. The `broad` set is the default (`--cluster_set`). | `sharp` (RECLU) needs replicates. |

## 2. SL labels of TSS clusters

| Parameter | Default | Rule |
| --- | --- | --- |
| `--sl_min_fraction` | 0.5 | A cluster is `sl` when SL READ1s / all READ1s (summed over samples) ≥ 0.5. Otherwise it is `nosl`. |

- **Separate clusters.** With cagescan's separate SL clusters (its default since `--sl_separate_clusters`), `SL…` clusters contain only SL reads and `TC…` clusters none, so the fractions are exactly 1 and 0. The threshold matters only for runs with pooled clusters. [data]
- **No SL table.** Every cluster is `nosl` when cagescan wrote no SL table, i.e. the organism has no trans-splicing or `--sl_sequence` was not given.
- **Several leaders.** When cagescan searched for several leaders, the main leader of each SL cluster is written as `sl_leader`. In *C. elegans*, SL2 marks downstream operon genes [lit: Blumenthal 2012].

## 3. 5' ends of coding mRNAs

| Parameter | Default | Rule and rationale |
| --- | --- | --- |
| `--tss_max_dist` | 500 | Candidate TSSs give a 5' UTR of at most 500 nt, measured along the transcript (5' UTR exons, then the genome upstream of the EGAPx 5' end). [prev] `add_five_prime_utr(dist = 500)`. |
| `--tss_selection` | `closest_atg` | Choose the candidate giving the shortest 5' UTR. [choice, user] In trans-spliced genes the mature 5' end is the SL acceptor, which lies downstream of the outron TSS, so "closest to the ATG" picks the acceptor. `dominant` (highest count) is the alternative. |
| `--tss_min_rel` | 0.1 | Candidates below 10% of the strongest reference peak are dropped. CAGE "exon painting" leaves weak capped 5' ends inside exons [lit: Fejes-Toth 2009], and these would otherwise win "closest to ATG". |
| `--tss_rel_whole_utr` | true | The reference for `--tss_min_rel` is every peak in the EGAPx 5' UTR, including those beyond `--tss_max_dist`. [choice] Otherwise a weak peak near the ATG could trim a UTR whose strong TSS lies further away. Such mRNAs keep the EGAPx end. |
| `--tss_min_count` | 0 | No absolute threshold. cagescan clusters already require ≥ 5 TPM in a sample. |
| `--utr5_min_len` | 0 | A TSS at the ATG is allowed. SL acceptors can sit at or next to the ATG: in *O. dioica* the leader-to-start distance averages 22 nt, and is ≤ 5 nt for 20% of genes [lit: Ganot 2004]. With the old default of 1, such genes fell back to EGAPx and lost their `sl` label. |

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

When the chosen TSS is still `nosl` but the transcript has a qualifying SL acceptor, the transcript records that acceptor as `sl_acceptor_peak` and `sl_acceptor_count`. The transcript also makes trans-spliced RNA, and this information is not lost.

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
| **5 nt up / 5 nt down (current)** | **4,859 (4,796)** | **1,535** | **18.6%** |

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

   The motif **closest to the EGAPx end** is used. The new end is the motif end + `--pas_offset` (0). [prev] The window and motif come from `add_three_prime_utr(dist = 500, motif = "AATAAA")`.
3. **`EGAPx`.** Otherwise the EGAPx end is kept.

Rationale:

- **Why closest to the EGAPx end, not first after the stop codon.** *O. dioica* is AT-rich, so AATAAA occurs by chance. "First after stop" would truncate real UTRs at chance motifs.
- **`--pas_offset 0`.** Cleavage usually happens 10–30 nt downstream of AATAAA. 0 keeps the old behaviour; the offset is exposed for tuning.
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
| `--operon_max_gap` | 500 | Same-strand genes join when the intergenic distance is below 500 nt (at most 499 bases between them). Genes on the other strand are ignored, and overlapping genes join. |
| `--operon_split_min_count` | 0 | Minimum count of a `nosl` TSS that splits an operon. |
| `--operon_biotypes` | `protein_coding` | Genes considered (`gene_biotype`, or "has a CDS" when the attribute is absent). |
| `--sl_window_up` / `--sl_window_down` | 5 / 5 | `nosl` TSSs this close to an SL TSS do not split (section 4). |

- **Grouping.** [lit, prev] Plessy et al. 2024 define operons as same-strand genes "separated by an intergenic distance of at most 500 bp". In their `calcOperons(window = 500)`, genes merge when start(B) − end(A) ≤ 500. Opposite-strand genes do not interrupt, and an operon needs at least 2 genes.
- **Gene ends.** Gene ends are those after UTR annotation. Plessy used CDS-only models, so gaps here are smaller and the same window is more permissive.
- **Splitting.** [prev, choice, user] A `nosl` TSS on the same strand between the end of the previous gene and a gene's start codon means the gene has its own promoter, so a new operon starts there, however short the gap. This is `split_operons()`. The SL acceptor rule (section 4) applies.
- **Benchmark.** [data, lit] On OKI2018_I69 v2.3 with the Okinawa CAGE consensus clusters, the Ganot et al. 2004 cyclin D3 operon (AJ628168) is recovered as one minus-strand operon at PAR:15956186–15967604: POC5 → IQCK → XPO7 (RBP16) → DYNLL2 → TUBD1 → MBF1 → CYCDC.
  - **Matches the literature:** Ganot's five genes are all there, in order, with the POC5 leader and IQCK that Plessy's analysis identified.
  - **Gaps:** 22–27 bp within the operon, against Ganot's 23–30 nt intercistronic regions.
  - **Neighbours:** PDCL3 upstream is split off by a `nosl` TSS of 783 reads.
  - **Robustness:** the result is the same for max gaps of 100–750 bp. Without the split it grows to 10 genes at 500 bp.
- **Genome-wide.**
  - OKI: about 2,973 operons with 46.6% of genes, against 3,124 operons and about 49% in Plessy's vignette.
  - O10 (current rules): 2,500 operons with 43.1% of genes.
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
- **Not yet tested:**
  - the I25 and O10 EGAPx annotations (the O10 EGAPx run was still in progress);
  - non-coding 5' ends (v1.3 has no lncRNAs);
  - nf-tests (test data pending; see [development.md](development.md)).

## 11. Known limitations and open questions

- **Tunicate dual-TSS genes.** Genes that make both trans-spliced and non-trans-spliced RNAs get one TSS per transcript model. The other form is recorded only as `sl_acceptor_peak` / `sl_acceptor_count`, and only when the chosen TSS is `nosl`. Making a second transcript model is a possible extension; KY2019 *Ciona* models do this [lit: Satou 2019].
- **Hybrid operons.** In *C. elegans* an operon can have an internal promoter and stay an operon [lit: Blumenthal 2012, 2015]. Our rule splits at any `nosl` TSS outside the SL window.
- **The DE1/DE2 artefact.** The cause is unknown. Reported upstream as johannesnicolaus/cagescan_nextflow#4.
- **PAS rules.** `--pas_offset 0` puts the end at the motif, not at the cleavage site. Chance AATAAA motifs are common in the AT-rich genome.
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
| 2026-10-18 | **SL acceptor rule v4: 5 nt upstream and 5 nt downstream (`--sl_window_up`, `--sl_window_down`); record `sl_acceptor_peak` when a `nosl` TSS wins** | Tunicates have genuine non-trans-spliced TSSs tens of nt upstream (*Ciona*, Matsumoto 2010) and downstream (*O. dioica* testis, Danks 2018) of acceptors. The *Ciona* TSS-seq studies' 51-nt minimum was a filter. The user preferred a conservative 5 nt. |

## References

All checked against PubMed unless noted.

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
