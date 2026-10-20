#!/usr/bin/env python3
"""Annotate operons: runs of closely spaced same-strand genes, split at non-SL TSSs.

  annotate_operons.py --gff PREFIX.utr.gff3 --tss PREFIX.tss.tsv --prefix PREFIX [options]

1. Genes of --biotypes (gene_biotype attribute; default protein_coding; genes without
   the attribute count as protein_coding when they have a CDS) are sorted along
   each strand of each sequence. Genes on the other strand are ignored, so nested or
   interleaved genes do not interrupt an operon (as calcOperons() in annotationpolish).
2. A gene joins the operon of the gene before it (in transcription order) when the
   intergenic distance, from the 3'-most stop codon of the operon so far to the gene's
   5'-most start codon, is below --max-gap-cds (500) bases (--gap-from cds, default).
   Distances between transcript ends would be the true intercistronic regions (23-30 nt
   in O. dioica), but annotated UTR ends are not precise enough: on the I69 EGAPx and O10
   annotations they gave less SL-consistent operons and broke the Ganot operon (see
   docs/rationale.md). --gap-from transcript (< --max-gap-tx, 100) and --gap-from
   evidence (transcript ends when both are supported by PAS/poly(A) reads and CAGE,
   otherwise CDS) are available. Non-coding genes use their transcript ends. Overlapping
   genes join.
3. An operon is split before a downstream gene when a non-SL (nosl) TSS with at least
   --split-min-count reads (no relative filter) lies on its strand between the 3' edge of
   the operon so far (the stop codon with --gap-from cds; the transcript 3' end with the
   other modes) and the gene's start codon (transcript 5' end for non-coding genes): the
   gene has its own promoter, however short the intergenic distance (as split_operons()).
   SL acceptor rule (as in annotate_utr5.py): a nosl TSS within --sl-window-up bases
   upstream or --sl-window-down bases downstream of an sl TSS with at least
   --split-min-count reads (default 5) does not split.
4. Runs of at least --min-genes (>= 2) genes are reported.

Outputs:
  PREFIX.operons.gff3       one 'operon' row (SO:0000178) per operon; attributes: genes
                            (transcription order), leader, n_genes, intergenic (bp)
  PREFIX.operons.bed        BED6; name = operon ID, score = number of genes
  PREFIX.operon_genes.tsv   one row per gene in a candidate run: operon (or NA after a
                            split left it alone), position, distance to the previous gene,
                            TSS type of the gene, whether an SL TSS lies between the previous
                            gene and its start codon, and the nosl TSS that split the run
  PREFIX.operon_summary_mqc.tsv  MultiQC table
"""

import argparse
import collections
import statistics
import sys

from gff3_utils import PointIndex, Annotation, escape, read_fasta


def read_tss_rows(path):
    items = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in handle:
            f = line.rstrip("\n").split("\t")
            items.append((f[col["chrom"]], f[col["strand"]], int(f[col["peak"]]),
                          {"tss_id": f[col["tss_id"]], "count": float(f[col["count"]]),
                           "tss_type": f[col["tss_type"]], "seqid": f[col["chrom"]], "strand": f[col["strand"]]}))
    return items


class GeneInfo:
    __slots__ = ("id", "name", "seqid", "strand", "start", "end", "atg", "tss_type", "sl_fraction",
                 "span_start", "span_end", "cds_start", "cds_end", "end5_evidence", "end3_evidence")

    def __init__(self, gene, gap_from="cds"):
        f = gene.feature
        self.id = f.id
        self.name = f.get("gene") or f.get("Name") or f.id
        self.seqid, self.strand, self.start, self.end = f.seqid, f.strand, f.start, f.end
        coding = [t for t in gene.transcripts if t.is_coding]
        if coding:
            # 5'-most start codon of the gene
            self.atg = min(t.cds_bounds()[0] for t in coding) if f.strand == "+" else max(t.cds_bounds()[1] for t in coding)
        else:
            self.atg = f.start if f.strand == "+" else f.end
        # CDS span (start codon to stop codon); the whole gene for non-coding genes
        if coding:
            self.cds_start = min(t.cds_bounds()[0] for t in coding)
            self.cds_end = max(t.cds_bounds()[1] for t in coding)
        else:
            self.cds_start, self.cds_end = f.start, f.end
        # span used for operon gaps with --gap-from cds / transcript
        if gap_from == "cds":
            self.span_start, self.span_end = self.cds_start, self.cds_end
        else:
            self.span_start, self.span_end = f.start, f.end
        # are the gene's 5' and 3' ends supported by evidence (set by this pipeline)?
        self.end5_evidence = any(t.feature.get("tss_evidence") == "CAGE" for t in gene.transcripts)
        self.end3_evidence = any(t.feature.get("utr3_evidence") in ("PAS", "polyA_reads") for t in gene.transcripts)
        types = [t.feature.get("tss_type") for t in gene.transcripts if t.feature.get("tss_type")]
        self.tss_type = "sl" if "sl" in types else ("nosl" if types else "NA")
        fractions = [t.feature.get("sl_fraction") for t in gene.transcripts if t.feature.get("sl_fraction")]
        self.sl_fraction = fractions[0] if fractions else "NA"

    @property
    def five_prime(self):
        return self.start if self.strand == "+" else self.end

    @property
    def three_prime(self):
        """3' end of the span used for operon gaps (stop codon with --gap-from cds)."""
        return self.span_end if self.strand == "+" else self.span_start


def pair_gap(prev_tx_edge, prev_cds_edge, prev_end3_ok, g, strand, gap_from, max_gap_tx, max_gap_cds):
    """(gap, joins) between the run so far and gene g.

    gap_from 'cds': stop-to-start codon distance < max_gap_cds.
    gap_from 'transcript': distance between transcript ends < max_gap_tx.
    gap_from 'evidence': the transcript-end distance (< max_gap_tx) when both ends are
    supported (upstream 3' end by PAS or poly(A) reads, downstream 5' end by CAGE),
    otherwise the CDS distance (< max_gap_cds).
    """
    if strand == "+":
        tx_gap, cds_gap = g.start - prev_tx_edge - 1, g.cds_start - prev_cds_edge - 1
    else:
        tx_gap, cds_gap = prev_tx_edge - g.end - 1, prev_cds_edge - g.cds_end - 1
    if gap_from == "transcript":
        return tx_gap, tx_gap < max_gap_tx
    if gap_from == "cds":
        return cds_gap, cds_gap < max_gap_cds
    if prev_end3_ok and g.end5_evidence:
        return tx_gap, tx_gap < max_gap_tx
    return cds_gap, cds_gap < max_gap_cds


def candidate_runs(genes, max_gap_tx, gap_from="cds", max_gap_cds=500):
    """Runs of same-strand genes closer than the gap limits, in transcription order.

    Yields lists of (gene, distance to the previous gene or None).
    """
    by_strand = collections.defaultdict(list)
    for g in genes:
        by_strand[(g.seqid, g.strand)].append(g)
    for (seqid, strand), items in sorted(by_strand.items()):
        items.sort(key=lambda g: (g.cds_start, g.cds_end))
        if strand == "-":
            items.sort(key=lambda g: (-g.cds_end, -g.cds_start))
        three = (lambda g: (g.end, g.cds_end)) if strand == "+" else (lambda g: (g.start, g.cds_start))
        further = max if strand == "+" else min
        run, tx_edge, cds_edge, end3_ok = [], None, None, False
        for g in items:
            if run:
                gap, joins = pair_gap(tx_edge, cds_edge, end3_ok, g, strand, gap_from, max_gap_tx, max_gap_cds)
                if joins:
                    run.append((g, gap))
                    tx3, cds3 = three(g)
                    if further(tx_edge, tx3) == tx3:
                        end3_ok = g.end3_evidence
                    tx_edge, cds_edge = further(tx_edge, tx3), further(cds_edge, cds3)
                    continue
                yield run
            run = [(g, None)]
            tx_edge, cds_edge = three(g)
            end3_ok = g.end3_evidence
        if run:
            yield run


def inside_acceptor_zone(pos, strand, sl_positions, window_up, window_down):
    """True if pos is within window_up bases upstream or window_down bases downstream of an
    sl TSS (-1: no limit on that side)."""
    for a in sl_positions:
        d = (a - pos) if strand == "+" else (pos - a)  # > 0: pos is upstream of the acceptor
        if d > 0 and (window_up < 0 or d <= window_up):
            return True
        if d <= 0 and (window_down < 0 or -d <= window_down):
            return True
    return False


def acceptor_dinucleotide(genome, seqid, strand, pos):
    """True if the two genomic bases before a TSS (in transcription order) are AG."""
    seq = genome.get(seqid)
    if seq is None:
        return False
    if strand == "+":
        return pos >= 3 and seq[pos - 3:pos - 1].upper() == "AG"
    return pos + 2 <= len(seq) and seq[pos:pos + 2].upper() == "CT"


def split_run(run, tss_index, split_min_count, window_up, window_down, rule="dominant", is_acceptor=None):
    """Split a run before genes that have their own promoter: a nosl TSS upstream of their start codon.

    The window runs from the 3' edge of the run so far (the stop codon with --gap-from cds,
    the transcript 3' end otherwise) to the gene's start codon. A nosl TSS is a candidate if
    its count is at least split_min_count and it is not in the acceptor zone of an sl TSS
    with count >= split_min_count (as the SL acceptor rule of ANNOTATE_UTR5).
      rule 'dominant': split only if the strongest candidate is stronger than every sl TSS
                       in the window (a gene mainly trans-spliced but with a minor promoter,
                       as in hybrid operons, stays in the operon);
      rule 'any':      split at any candidate;
      rule 'none':     never split.
    is_acceptor(row, pos): which TSSs count as SL acceptors; by default the sl label. Without
    SL reads, AG-preceded TSSs are taken as undetected acceptors.

    Returns (list of operons, {gene id: splitting TSS row}, {gene id: sl TSS upstream}).
    """
    is_acceptor = is_acceptor or (lambda row, pos: row["tss_type"] == "sl")
    pieces, splits, sl_upstream = [[run[0]]], {}, {}
    edge = run[0][0].three_prime
    for g, gap in run[1:]:
        if g.strand == "+":
            lo, hi = edge + 1, g.atg
        else:
            lo, hi = g.atg, edge - 1
        # acceptors just outside the window can still shadow a nosl TSS inside it
        pad = hi - lo + 1 if min(window_up, window_down) < 0 else max(window_up, window_down, 0)
        nearby = [(p, row, is_acceptor(row, p))
                  for p, row in tss_index.window(g.seqid, g.strand, lo - pad, hi + pad) if row["count"] > 0]
        inside = [(p, row, acc) for p, row, acc in nearby if lo <= p <= hi]
        acceptors = [p for p, row, acc in nearby if acc and row["count"] >= split_min_count]
        nosl = [row for p, row, acc in inside
                if not acc and row["count"] >= split_min_count
                and not inside_acceptor_zone(p, g.strand, acceptors, window_up, window_down)]
        sl_inside = [row["count"] for _, row, acc in inside if acc]
        sl_upstream[g.id] = bool(sl_inside)
        best = max(nosl, key=lambda r: r["count"]) if nosl else None
        if rule == "none" or (rule == "dominant" and best and best["count"] <= max(sl_inside, default=0)):
            best = None
        if best:
            splits[g.id] = best
            pieces.append([(g, gap)])
        else:
            pieces[-1].append((g, gap))
        edge = max(edge, g.span_end) if g.strand == "+" else min(edge, g.span_start)
    return pieces, splits, sl_upstream


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tss", required=True, help="table from cage_tss_table.py")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--max-gap-cds", type=int, default=500,
                        help="stop-to-start codon distances below this join genes (--gap-from cds, evidence)")
    parser.add_argument("--max-gap-tx", type=int, default=100,
                        help="distances between transcript ends below this join genes (--gap-from transcript, evidence)")
    parser.add_argument("--split-min-count", type=float, default=0, help="minimum pooled count of a splitting nosl TSS")
    parser.add_argument("--split-rule", choices=["dominant", "any", "none"], default="dominant",
                        help="dominant: split at a nosl TSS stronger than every sl TSS between the genes (default); "
                             "any: at any nosl TSS; none: distance only")
    parser.add_argument("--fasta", help="genome; without SL reads in the TSS table, AG-preceded TSSs are taken as "
                                        "undetected SL acceptors (without it, operons are not split)")
    parser.add_argument("--min-genes", type=int, default=2)
    parser.add_argument("--gap-from", choices=["cds", "transcript", "evidence"], default="cds",
                        help="cds: stop to start codon (default); transcript: transcript ends; evidence: "
                             "transcript ends when both are supported by evidence, otherwise cds")
    parser.add_argument("--sl-window-up", type=int, default=5,
                        help="a nosl TSS up to this many bases upstream of an sl TSS does not split (-1: no limit)")
    parser.add_argument("--sl-window-down", type=int, default=5,
                        help="a nosl TSS up to this many bases downstream of an sl TSS does not split (-1: no limit)")
    parser.add_argument("--biotypes", default="protein_coding",
                        help="comma-separated gene_biotype values of the genes considered ('all' for every gene)")
    args = parser.parse_args()
    if args.min_genes < 2:
        parser.error("--min-genes must be at least 2")

    ann = Annotation.read(args.gff)
    tss_items = read_tss_rows(args.tss)
    tss_index = PointIndex(tss_items)
    biotypes = None if args.biotypes == "all" else set(args.biotypes.split(","))
    def biotype(gene):
        value = gene.feature.get("gene_biotype")
        if value:
            return value
        return "protein_coding" if any(t.is_coding for t in gene.transcripts) else "other"

    genes = [GeneInfo(g, args.gap_from) for g in ann.genes
             if g.transcripts and g.feature.type == "gene" and (biotypes is None or biotype(g) in biotypes)]
    sl_source, is_acceptor, split_rule = "SL_reads", None, args.split_rule
    if not any(row["tss_type"] == "sl" and row["count"] > 0 for _, _, _, row in tss_items):
        # SL data were given but no SL read was recognised (e.g. a wrong leader sequence): the
        # trans-splice acceptors then appear as nosl TSSs. Acceptors follow an AG (96-99% of sl
        # TSSs on O10, I69 and C. elegans, against 2-13% of nosl TSSs), so AG-preceded TSSs are
        # taken as undetected acceptors; without the genome, operons are not split.
        if args.fasta:
            genome = dict(read_fasta(args.fasta))
            sl_source = "AG_inferred"
            is_acceptor = lambda row, pos, _g=genome: acceptor_dinucleotide(_g, row["seqid"], row["strand"], pos)
            print("WARNING: no sl TSS with reads in the TSS table (was the spliced leader detected?); "
                  "operons use AG-preceded TSSs as SL acceptors (sl_source=AG_inferred)", file=sys.stderr)
        else:
            sl_source, split_rule = "none", "none"
            print("WARNING: no sl TSS with reads in the TSS table (was the spliced leader detected?); "
                  "operons are called from distances only, without splits (sl_source=none)", file=sys.stderr)

    operons, rows = [], []
    for run in candidate_runs(genes, args.max_gap_tx, args.gap_from, args.max_gap_cds):
        if len(run) < 2:
            continue
        pieces, splits, sl_upstream = split_run(run, tss_index, args.split_min_count,
                                                args.sl_window_up, args.sl_window_down, split_rule, is_acceptor)
        for piece in pieces:
            operon_index = len(operons) if len(piece) >= args.min_genes else None
            if operon_index is not None:
                operons.append([(g, gap, sl_upstream.get(g.id)) for g, gap in piece])
            for pos, (g, gap) in enumerate(piece, 1):
                rows.append([operon_index, pos, g, gap, splits.get(g.id), sl_upstream.get(g.id)])

    # number operons along the genome
    seq_rank = {}
    for g in genes:  # sequences in the order of the annotation
        seq_rank.setdefault(g.seqid, len(seq_rank))
    order = sorted(range(len(operons)), key=lambda i: (seq_rank[operons[i][0][0].seqid], min(g.start for g, _, _ in operons[i])))
    width = max(5, len(str(len(operons))))
    names = [None] * len(operons)
    for rank, i in enumerate(order, 1):
        names[i] = f"operon_{rank:0{width}d}"

    evidence_counts = collections.Counter()
    with open(f"{args.prefix}.operons.gff3", "w") as gff, open(f"{args.prefix}.operons.bed", "w") as bed:
        gff.write("##gff-version 3\n")
        for name, piece in sorted(zip(names, operons)):
            first = piece[0][0]
            start = min(g.start for g, _, _ in piece)
            end = max(g.end for g, _, _ in piece)
            n_sl = sum(1 for _, _, up in piece[1:] if up)
            n_down = len(piece) - 1
            evidence = "distance" if n_sl == 0 else ("sl" if n_sl == n_down else "partial_sl")
            attrs = [  # list values are comma-separated GFF3 multi-values
                ("ID", [name]),
                ("n_genes", [len(piece)]),
                ("genes", [g.id for g, _, _ in piece]),
                ("leader", [first.id]),
                ("intergenic", [gap for _, gap, _ in piece[1:]]),
                ("tss_types", [g.tss_type for g, _, _ in piece]),
                ("sl_upstream", [str(bool(up)).lower() for _, _, up in piece[1:]]),
                ("downstream_sl", [f"{n_sl}/{n_down}"]),
                ("operon_evidence", [evidence]),
                ("sl_source", [sl_source]),
            ]
            gff.write("\t".join([first.seqid, "utrannotation", "operon", str(start), str(end), ".", first.strand, ".",
                                 ";".join(f"{k}={','.join(escape(v) for v in vs)}" for k, vs in attrs)]) + "\n")
            bed.write(f"{first.seqid}\t{start - 1}\t{end}\t{name}\t{len(piece)}\t{first.strand}\n")
            evidence_counts[evidence] += 1

    with open(f"{args.prefix}.operon_genes.tsv", "w") as out:
        out.write("operon_id\tposition\tgene_id\tgene_name\tseqid\tstrand\tstart\tend\t"
                  "distance_to_previous\ttss_type\tsl_fraction\tsl_tss_upstream\tsplit_by_tss\tsplit_tss_count\n")
        for operon_index, pos, g, gap, split, sl_up in rows:
            operon = names[operon_index] if operon_index is not None else "NA"
            out.write("\t".join(str(x) for x in [
                operon, pos if operon_index is not None else "NA", g.id, g.name, g.seqid, g.strand, g.start, g.end,
                "NA" if gap is None else gap, g.tss_type, g.sl_fraction, "NA" if sl_up is None else str(sl_up).lower(),
                split["tss_id"] if split else "NA", f"{split['count']:g}" if split else "NA",
            ]) + "\n")

    in_operons = [(pos, g) for operon_index, pos, g, _, _, _ in rows if operon_index is not None]
    downstream_sl = [sl_up for operon_index, pos, _, _, _, sl_up in rows if operon_index is not None and pos > 1]
    split_genes_sl = [sl_up for _, _, _, _, split, sl_up in rows if split]
    sizes = [len(p) for p in operons]
    pct = lambda n, d: f"{100 * n / d:.1f}" if d else "NA"
    with open(f"{args.prefix}.operon_summary_mqc.tsv", "w") as out:
        out.write("# id: 'operon_summary'\n")
        out.write("# section_name: 'Operons'\n")
        out.write('# description: "Runs of same-strand genes closer than the maximum intergenic distance, split '
                  'before genes with a non-SL TSS. In true operons the downstream genes are trans-spliced (sl)."\n')
        out.write("# plot_type: 'table'\n")
        out.write("# pconfig:\n#     id: 'operon_summary_table'\n#     title: 'Operons'\n")
        out.write("Sample\tGenes considered\tOperons\tGenes in operons (%)\tMedian genes per operon\t"
                  "Largest operon\tOperons with SL at every downstream gene\tOperons with SL at some\t"
                  "Operons by distance only\tSplits at nosl TSS\tDownstream genes with SL TSS (%)\t"
                  "Split genes with SL TSS (%)\tSL source\n")
        out.write("\t".join(str(x) for x in [
            args.prefix, len(genes), len(operons), pct(len(in_operons), len(genes)),
            f"{statistics.median(sizes):g}" if sizes else "NA", max(sizes) if sizes else "NA",
            evidence_counts["sl"], evidence_counts["partial_sl"], evidence_counts["distance"],
            sum(1 for r in rows if r[4]),
            pct(sum(downstream_sl), len(downstream_sl)),
            pct(sum(split_genes_sl), len(split_genes_sl)),
            sl_source,
        ]) + "\n")

    print(f"{len(genes)} genes, {len(operons)} operons with {len(in_operons)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
