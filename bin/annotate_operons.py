#!/usr/bin/env python3
"""Annotate operons: runs of closely spaced same-strand genes, split at non-SL TSSs.

  annotate_operons.py --gff PREFIX.utr.gff3 --tss PREFIX.tss.tsv --prefix PREFIX [options]

1. Genes of --biotypes (gene_biotype attribute; default protein_coding; genes without
   the attribute count as protein_coding when they have a CDS) are sorted along
   each strand of each sequence. Genes on the other strand are ignored, so nested or
   interleaved genes do not interrupt an operon (as calcOperons() in annotationpolish).
2. A gene joins the operon of the gene before it (in transcription order) when the
   intergenic distance, from the 3'-most end of the operon so far to the 5' end of the
   gene, is below --max-gap bases. Overlapping genes have a negative distance and join.
   Gene ends are those of the input annotation, i.e. after the CAGE 5' and PAS 3' ends.
3. An operon is split before a downstream gene when a non-SL (nosl) TSS with at least
   --split-min-count reads lies on its strand between the end of the previous gene and
   the gene's start codon (its first transcript base for non-coding genes): the gene
   has its own promoter, however short the intergenic distance (as split_operons()).
   A nosl TSS within --sl-shadow bases of a stronger sl TSS is taken as SL reads whose
   leader was not recognised, and does not split.
4. Runs of at least --min-genes genes are reported.

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

from gff3_utils import PointIndex, Annotation, escape


def read_tss(path):
    items = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in handle:
            f = line.rstrip("\n").split("\t")
            items.append((f[col["chrom"]], f[col["strand"]], int(f[col["peak"]]),
                          {"tss_id": f[col["tss_id"]], "count": float(f[col["count"]]),
                           "tss_type": f[col["tss_type"]]}))
    return PointIndex(items)


class GeneInfo:
    __slots__ = ("id", "name", "seqid", "strand", "start", "end", "atg", "tss_type", "sl_fraction")

    def __init__(self, gene):
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
        types = [t.feature.get("tss_type") for t in gene.transcripts if t.feature.get("tss_type")]
        self.tss_type = "sl" if "sl" in types else ("nosl" if types else "NA")
        fractions = [t.feature.get("sl_fraction") for t in gene.transcripts if t.feature.get("sl_fraction")]
        self.sl_fraction = fractions[0] if fractions else "NA"

    @property
    def five_prime(self):
        return self.start if self.strand == "+" else self.end

    @property
    def three_prime(self):
        return self.end if self.strand == "+" else self.start


def candidate_runs(genes, max_gap):
    """Runs of same-strand genes with intergenic distances below max_gap, in transcription order.

    Yields lists of (gene, distance to the previous gene or None).
    """
    by_strand = collections.defaultdict(list)
    for g in genes:
        by_strand[(g.seqid, g.strand)].append(g)
    for (seqid, strand), items in sorted(by_strand.items()):
        items.sort(key=lambda g: (g.start, g.end))
        if strand == "-":
            items.sort(key=lambda g: (-g.end, -g.start))
        run, edge = [], None  # edge: 3'-most coordinate of the run so far
        for g in items:
            if run:
                gap = (g.start - edge - 1) if strand == "+" else (edge - g.end - 1)
                if gap < max_gap:
                    run.append((g, gap))
                    edge = max(edge, g.end) if strand == "+" else min(edge, g.start)
                    continue
                yield run
            run, edge = [(g, None)], g.three_prime
        if run:
            yield run


def shadowed(pos, row, upstream_rows, shadow):
    """True if a stronger sl TSS lies within `shadow` bases of this nosl TSS."""
    return any(r["tss_type"] == "sl" and abs(p - pos) <= shadow and r["count"] > row["count"]
               for p, r in upstream_rows)


def split_run(run, tss_index, split_min_count, shadow):
    """Split a run before genes that have a nosl TSS upstream of their start codon.

    Returns (list of operons, {gene id: splitting TSS row}, {gene id: sl TSS upstream}).
    """
    pieces, splits, sl_upstream = [[run[0]]], {}, {}
    edge = run[0][0].three_prime
    for g, gap in run[1:]:
        if g.strand == "+":
            lo, hi = edge + 1, g.atg
        else:
            lo, hi = g.atg, edge - 1
        upstream = [(p, row) for p, row in tss_index.window(g.seqid, g.strand, lo - shadow, hi + shadow)
                    if row["count"] > 0]
        inside = [(p, row) for p, row in upstream if lo <= p <= hi]
        nosl = [row for p, row in inside
                if row["tss_type"] == "nosl" and row["count"] >= split_min_count
                and not shadowed(p, row, upstream, shadow)]
        sl_upstream[g.id] = any(row["tss_type"] == "sl" for _, row in inside)
        if nosl:
            splits[g.id] = max(nosl, key=lambda r: r["count"])
            pieces.append([(g, gap)])
        else:
            pieces[-1].append((g, gap))
        edge = max(edge, g.end) if g.strand == "+" else min(edge, g.start)
    return pieces, splits, sl_upstream


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tss", required=True, help="table from cage_tss_table.py")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--max-gap", type=int, default=500, help="intergenic distances below this join genes")
    parser.add_argument("--split-min-count", type=float, default=0, help="minimum pooled count of a splitting nosl TSS")
    parser.add_argument("--min-genes", type=int, default=2)
    parser.add_argument("--sl-shadow", type=int, default=5,
                        help="a nosl TSS this close to a stronger sl TSS does not split (0 to disable)")
    parser.add_argument("--biotypes", default="protein_coding",
                        help="comma-separated gene_biotype values of the genes considered ('all' for every gene)")
    args = parser.parse_args()

    ann = Annotation.read(args.gff)
    tss_index = read_tss(args.tss)
    biotypes = None if args.biotypes == "all" else set(args.biotypes.split(","))
    def biotype(gene):
        value = gene.feature.get("gene_biotype")
        if value:
            return value
        return "protein_coding" if any(t.is_coding for t in gene.transcripts) else "other"

    genes = [GeneInfo(g) for g in ann.genes
             if g.transcripts and g.feature.type == "gene" and (biotypes is None or biotype(g) in biotypes)]

    operons, rows = [], []
    for run in candidate_runs(genes, args.max_gap):
        if len(run) < 2:
            continue
        pieces, splits, sl_upstream = split_run(run, tss_index, args.split_min_count, args.sl_shadow)
        for piece in pieces:
            operon_index = len(operons) if len(piece) >= args.min_genes else None
            if operon_index is not None:
                operons.append(piece)
            for pos, (g, gap) in enumerate(piece, 1):
                rows.append([operon_index, pos, g, gap, splits.get(g.id), sl_upstream.get(g.id)])

    # number operons along the genome
    order = sorted(range(len(operons)), key=lambda i: (operons[i][0][0].seqid, min(g.start for g, _ in operons[i])))
    width = max(5, len(str(len(operons))))
    names = [None] * len(operons)
    for rank, i in enumerate(order, 1):
        names[i] = f"operon_{rank:0{width}d}"

    with open(f"{args.prefix}.operons.gff3", "w") as gff, open(f"{args.prefix}.operons.bed", "w") as bed:
        gff.write("##gff-version 3\n")
        for name, piece in sorted(zip(names, operons)):
            first = piece[0][0]
            start = min(g.start for g, _ in piece)
            end = max(g.end for g, _ in piece)
            attrs = [  # list values are comma-separated GFF3 multi-values
                ("ID", [name]),
                ("n_genes", [len(piece)]),
                ("genes", [g.id for g, _ in piece]),
                ("leader", [first.id]),
                ("intergenic", [gap for _, gap in piece[1:]]),
                ("tss_types", [g.tss_type for g, _ in piece]),
            ]
            gff.write("\t".join([first.seqid, "utrannotation", "operon", str(start), str(end), ".", first.strand, ".",
                                 ";".join(f"{k}={','.join(escape(v) for v in vs)}" for k, vs in attrs)]) + "\n")
            bed.write(f"{first.seqid}\t{start - 1}\t{end}\t{name}\t{len(piece)}\t{first.strand}\n")

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
                  "Largest operon\tSplits at nosl TSS\tDownstream genes with SL TSS (%)\tSplit genes with SL TSS (%)\n")
        out.write("\t".join(str(x) for x in [
            args.prefix, len(genes), len(operons), pct(len(in_operons), len(genes)),
            f"{statistics.median(sizes):g}" if sizes else "NA", max(sizes) if sizes else "NA",
            sum(1 for r in rows if r[4]),
            pct(sum(downstream_sl), len(downstream_sl)),
            pct(sum(split_genes_sl), len(split_genes_sl)),
        ]) + "\n")

    print(f"{len(genes)} genes, {len(operons)} operons with {len(in_operons)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
