#!/usr/bin/env python3
"""Move the 5' end of each coding mRNA to a CAGE TSS and annotate its 5' UTR.

  annotate_utr5.py --gff annotation.gff3 --tss PREFIX.tss.tsv --fasta genome.fa \\
                   --output PREFIX.utr5.gff3 [options]

For every mRNA with a CDS (models open at the 5' end are left alone):

1. Candidates are TSS peaks on the same strand that
   - give a 5' UTR of --min-utr-len to --max-dist bases, measured along the transcript
     (5' UTR exons, then genomic sequence upstream of the annotated 5' end);
   - are not in an intron of the 5' UTR;
   - do not make the transcript overlap a CDS of another gene on the same strand.
2. Peaks weaker than --min-rel x the strongest candidate, or below --min-count, are
   dropped (CAGE "exon painting" leaves weak peaks inside exons).
3. One peak is chosen:
   closest_atg  the shortest 5' UTR; in trans-spliced genes this is the SL acceptor,
                which lies downstream of the outron TSS (default);
   dominant     the highest count.
   Ties prefer nosl over sl (as add_five_prime_utr() did), then the higher count.
4. Exons upstream of the TSS are dropped and the first exon starts at the TSS, so the
   UTR is extended or trimmed. CDS rows are never changed.

Attributes added to the mRNA, and copied onto its five_prime_UTR rows:
  utr5_evidence  CAGE, or EGAPx when no peak qualified (EGAPx 5' end kept)
  tss_type       sl (spliced-leader acceptor) or nosl
  tss_id, tss_peak (chrom:pos:strand), tss_count, sl_count, sl_fraction
  utr5_shift     new minus old 5' UTR length
  acceptor_AG    sl only: true if the two genomic bases before the TSS are AG
An sl TSS also gets a 1-bp trans_splice_acceptor_site row (SO:0000706).
"""

import argparse
import collections
import sys

from gff3_utils import (TSA, Annotation, Feature, PointIndex, containing_interval, exonic_length_between,
                        feature_core_id, read_fasta)


def read_tss(path):
    items = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for name in ("chrom", "peak", "strand", "tss_id", "count", "sl_count", "sl_fraction", "tss_type"):
            if name not in col:
                sys.exit(f"ERROR: column {name} missing from {path}")
        for line in handle:
            f = line.rstrip("\n").split("\t")
            row = {
                "tss_id": f[col["tss_id"]],
                "count": float(f[col["count"]]),
                "sl_count": float(f[col["sl_count"]]),
                "sl_fraction": f[col["sl_fraction"]],
                "tss_type": f[col["tss_type"]],
            }
            items.append((f[col["chrom"]], f[col["strand"]], int(f[col["peak"]]), row))
    return PointIndex(items)


def candidates(t, gene_id, tss_index, cds_index, max_dist, min_utr_len):
    """(pos, tss_row, utr_len) for every positional candidate of transcript t."""
    exons = t.exon_intervals()
    cds_start, cds_end = t.cds_bounds()
    key = (t.seqid, t.strand)
    out = []
    if t.strand == "+":
        first = exons[0][0]
        exonic_utr = exonic_length_between(exons, first, cds_start - 1)
        hits = tss_index.window(t.seqid, "+", first - max_dist, cds_start - min_utr_len)
        for pos, row in hits:
            if pos < first:
                utr_len = exonic_utr + (first - pos)
                guard = (pos, first - 1)
            else:
                if containing_interval(exons, pos) is None:
                    continue  # in a 5' UTR intron
                utr_len = exonic_length_between(exons, pos, cds_start - 1)
                guard = (pos, pos)
            if min_utr_len <= utr_len <= max_dist and not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, utr_len))
    else:
        last = exons[-1][1]
        exonic_utr = exonic_length_between(exons, cds_end + 1, last)
        hits = tss_index.window(t.seqid, "-", cds_end + min_utr_len, last + max_dist)
        for pos, row in hits:
            if pos > last:
                utr_len = exonic_utr + (pos - last)
                guard = (last + 1, pos)
            else:
                if containing_interval(exons, pos) is None:
                    continue
                utr_len = exonic_length_between(exons, cds_end + 1, pos)
                guard = (pos, pos)
            if min_utr_len <= utr_len <= max_dist and not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, utr_len))
    return out


def choose(cands, selection, min_rel, min_count):
    if not cands:
        return None
    best_count = max(row["count"] for _, row, _ in cands)
    threshold = max(min_count, min_rel * best_count)
    kept = [c for c in cands if c[1]["count"] >= threshold and c[1]["count"] > 0]
    if not kept:
        return None
    nosl_first = lambda c: 0 if c[1]["tss_type"] == "nosl" else 1
    if selection == "dominant":
        return min(kept, key=lambda c: (-c[1]["count"], c[2], nosl_first(c)))
    return min(kept, key=lambda c: (c[2], nosl_first(c), -c[1]["count"]))


def set_five_prime_end(t, pos):
    if t.strand == "+":
        kept = [e for e in t.exons if e.end >= pos]
        kept[0].start = pos
    else:
        kept = [e for e in t.exons if e.start <= pos]
        kept[-1].end = pos
    t.exons = kept
    t.sync_bounds()


def acceptor_ag(genome, seqid, strand, pos):
    seq = genome.get(seqid)
    if seq is None:
        return None
    if strand == "+":
        return pos >= 3 and seq[pos - 3:pos - 1].upper() == "AG"
    return seq[pos:pos + 2].upper() == "CT"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tss", required=True, help="table from cage_tss_table.py")
    parser.add_argument("--fasta", required=True, help="genome, for the AG check of SL acceptors")
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-dist", type=int, default=500, help="maximum 5' UTR length (bases along the transcript)")
    parser.add_argument("--selection", choices=["closest_atg", "dominant"], default="closest_atg")
    parser.add_argument("--min-rel", type=float, default=0.1, help="minimum count relative to the strongest candidate")
    parser.add_argument("--min-count", type=float, default=0, help="minimum pooled count of a TSS")
    parser.add_argument("--min-utr-len", type=int, default=1, help="minimum 5' UTR length")
    args = parser.parse_args()

    ann = Annotation.read(args.gff)
    tss_index = read_tss(args.tss)
    cds_index = ann.cds_index()
    genome = None

    stats = collections.Counter()
    for gene, t in ann.coding_transcripts():
        stats["mRNA"] += 1
        f = t.feature
        # drop acceptor rows from an earlier run, so the script can be re-applied
        t.others = [o for o in t.others if o.type != TSA]
        if t.five_prime_partial():
            stats["5' partial, kept"] += 1
            f.set("utr5_evidence", "EGAPx")
            continue
        cds_before = [(c.start, c.end) for c in t.cds]
        old_len = t.utr_lengths()[0]
        pick = choose(
            candidates(t, gene.feature.id, tss_index, cds_index, args.max_dist, args.min_utr_len),
            args.selection, args.min_rel, args.min_count,
        )
        if pick is None:
            stats["no CAGE TSS, EGAPx kept"] += 1
            f.set("utr5_evidence", "EGAPx")
            continue

        pos, row, utr_len = pick
        set_five_prime_end(t, pos)
        assert [(c.start, c.end) for c in t.cds] == cds_before, f"CDS changed in {t.id}"
        assert t.utr_lengths()[0] == utr_len, f"5' UTR length mismatch in {t.id}"
        shift = utr_len - old_len
        stats["CAGE TSS"] += 1
        stats["extended" if shift > 0 else "trimmed" if shift < 0 else "unchanged"] += 1

        f.set("utr5_evidence", "CAGE")
        f.set("tss_type", row["tss_type"])
        f.set("tss_id", row["tss_id"])
        f.set("tss_peak", f"{t.seqid}:{pos}:{t.strand}")
        f.set("tss_count", f"{row['count']:g}")
        f.set("sl_count", f"{row['sl_count']:g}")
        f.set("sl_fraction", row["sl_fraction"])
        f.set("utr5_shift", shift)
        if row["tss_type"] == "sl":
            stats["sl"] += 1
            if genome is None:
                genome = dict(read_fasta(args.fasta))
            ag = acceptor_ag(genome, t.seqid, t.strand, pos)
            if ag is not None:
                f.set("acceptor_AG", "true" if ag else "false")
                stats["sl with AG"] += bool(ag)
            site = Feature(t.seqid, "CAGE", TSA, pos, pos, ".", t.strand, ".", collections.OrderedDict())
            site.set("ID", f"tsa-{feature_core_id(t)}")
            site.set("Parent", t.id)
            site.set("tss_id", row["tss_id"])
            site.set("sl_fraction", row["sl_fraction"])
            t.others.append(site)

    ann.write(args.output)
    for key, value in stats.items():
        print(f"{key}\t{value}", file=sys.stderr)


if __name__ == "__main__":
    main()
