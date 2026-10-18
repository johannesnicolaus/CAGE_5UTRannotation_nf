#!/usr/bin/env python3
"""Move the 5' end of each transcript to a CAGE TSS (5' UTRs of coding mRNAs).

  annotate_utr5.py --gff annotation.gff3 --tss PREFIX.tss.tsv --fasta genome.fa \\
                   --output PREFIX.utr5.gff3 [options]

Coding mRNAs (models open at the 5' end are left alone):

1. Candidates are TSS peaks on the same strand that
   - give a 5' UTR of --min-utr-len to --max-dist bases, measured along the transcript
     (5' UTR exons, then genomic sequence upstream of the annotated 5' end);
   - are not in an intron of the 5' UTR;
   - do not make the transcript overlap a CDS of another gene on the same strand.
2. Peaks below --min-count, or weaker than --min-rel x the strongest reference peak, are
   dropped (CAGE "exon painting" leaves weak peaks inside exons). The reference peaks
   are the candidates plus every peak in the EGAPx 5' UTR, also beyond --max-dist, so a
   weak peak near the ATG cannot trim a UTR whose strong TSS lies further away
   (--rel-window-only compares with the candidates only).
3. One peak is chosen:
   closest_atg  the shortest 5' UTR; in trans-spliced genes this is the SL acceptor,
                which lies downstream of the outron TSS (default);
   dominant     the highest count.
   Ties prefer nosl over sl (as add_five_prime_utr() did), then the higher count.
   A nosl peak is not a candidate (an "SL shadow") when it lies exactly 2 bases upstream
   of an sl peak, on the A of the acceptor AG, or within --sl-shadow bases of a stronger
   sl peak. Such peaks come from reads that start with the leader's last two bases (AG)
   and no other leader sequence; they align to the acceptor AG and look like a TSS 2
   bases upstream. So an SL acceptor with such a shadow stays sl.
4. Exons upstream of the TSS are dropped and the first exon starts at the TSS, so the
   UTR is extended or trimmed. CDS rows are never changed.

Non-coding transcripts (lncRNA; 'transcript', i.e. EGAPx misc_RNA; ncRNA) have no ATG:
candidates lie within --max-dist of the EGAPx 5' end, upstream (genome) or downstream
(exons, keeping at least --nc-min-length bases of transcript), and the peak with the
highest count is used (ties: the smallest move). --min-count applies.

Attributes added to the transcript (five_prime_UTR rows carry only ID and Parent):
  tss_evidence   CAGE, or EGAPx when no peak qualified (EGAPx 5' end kept)
  tss_type       sl (spliced-leader acceptor) or nosl
  tss_id, tss_peak (chrom:pos:strand), tss_count, sl_count, sl_fraction
  tss_shift      bases added (positive) or removed (negative) at the 5' end
  acceptor_AG    sl only: true if the two genomic bases before the TSS are AG
  sl_leader      sl only, when several leaders were searched for: the main leader
An sl TSS also gets a 1-bp trans_splice_acceptor_site row (SO:0000706). When the TSS of
an mRNA is at the ATG (5' UTR of length 0, common for SL acceptors) there is no
five_prime_UTR row; the attributes and the acceptor row still record it.

The script stops if annotation sequences are missing from the genome, or if TSSs fall
outside the genome's sequences (cagescan run on another assembly).
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
                "sl_leader": f[col["sl_leader"]] if "sl_leader" in col else "NA",
            }
            items.append((f[col["chrom"]], f[col["strand"]], int(f[col["peak"]]), row))
    return items


def check_sequences(ann, genome, tss_items, tss_path):
    """Stop on sequence-name or coordinate mismatches between annotation, genome and TSSs."""
    missing = sorted(ann.seqids() - set(genome))
    if missing:
        sys.exit(f"ERROR: annotation sequences not in the genome FASTA: {', '.join(missing[:10])}")
    if not tss_items:
        print(f"WARNING: no TSS in {tss_path}; every 5' end is kept", file=sys.stderr)
        return
    off_genome = [i for i in tss_items if i[0] not in genome]
    beyond_end = [i for i in tss_items if i[0] in genome and i[2] > len(genome[i[0]])]
    if beyond_end:
        chrom, _, pos, _ = beyond_end[0]
        sys.exit(f"ERROR: {len(beyond_end)} TSSs lie beyond the end of their sequence (e.g. {chrom}:{pos}); "
                 "was cagescan run on another assembly?")
    if len(off_genome) == len(tss_items):
        sys.exit("ERROR: no TSS is on a genome sequence (TSS sequences e.g. "
                 f"{', '.join(sorted({i[0] for i in tss_items})[:5])}); was cagescan run on another assembly?")
    if off_genome:
        print(f"WARNING: {len(off_genome)} of {len(tss_items)} TSSs are on sequences absent from the genome",
              file=sys.stderr)


def coding_candidates(t, gene_id, tss_index, cds_index, max_dist, min_utr_len):
    """(pos, tss_row, utr_len, in_window) for peaks in the EGAPx 5' UTR or within max_dist upstream.

    in_window: the peak gives a 5' UTR of min_utr_len..max_dist bases (a selection candidate);
    the others are exonic peaks further than max_dist from the ATG, used only as reference.
    """
    exons = t.exon_intervals()
    cds_start, cds_end = t.cds_bounds()
    key = (t.seqid, t.strand)
    out = []
    if t.strand == "+":
        first = exons[0][0]
        exonic_utr = exonic_length_between(exons, first, cds_start - 1)
        for pos, row in tss_index.window(t.seqid, "+", first - max_dist, cds_start - min_utr_len):
            if pos < first:
                utr_len = exonic_utr + (first - pos)
                guard = (pos, first - 1)
            else:
                if containing_interval(exons, pos) is None:
                    continue  # in a 5' UTR intron
                utr_len = exonic_length_between(exons, pos, cds_start - 1)
                guard = (pos, pos)
            in_window = min_utr_len <= utr_len <= max_dist
            if (in_window or pos >= first) and not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, utr_len, in_window))
    else:
        last = exons[-1][1]
        exonic_utr = exonic_length_between(exons, cds_end + 1, last)
        for pos, row in tss_index.window(t.seqid, "-", cds_end + min_utr_len, last + max_dist):
            if pos > last:
                utr_len = exonic_utr + (pos - last)
                guard = (last + 1, pos)
            else:
                if containing_interval(exons, pos) is None:
                    continue
                utr_len = exonic_length_between(exons, cds_end + 1, pos)
                guard = (pos, pos)
            in_window = min_utr_len <= utr_len <= max_dist
            if (in_window or pos <= last) and not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, utr_len, in_window))
    return out


def is_sl_shadow(pos, row, strand, nearby, shadow):
    """True for a nosl peak on the A of an SL acceptor (2 bases upstream of an sl peak), or
    within `shadow` bases of a stronger sl peak. `nearby`: (pos, row) within `shadow` bases."""
    if row["tss_type"] != "nosl" or shadow <= 0:
        return False
    acceptor_a = pos + 2 if strand == "+" else pos - 2
    return any(r["tss_type"] == "sl" and (p == acceptor_a or r["count"] > row["count"]) for p, r in nearby)


def drop_sl_shadows(peaks, tss_index, seqid, strand, shadow):
    """Remove nosl peaks that are SL shadows (leftover SL reads at an acceptor)."""
    if shadow <= 0:
        return peaks
    return [p for p in peaks
            if not is_sl_shadow(p[0], p[1], strand,
                                tss_index.window(seqid, strand, p[0] - max(shadow, 2), p[0] + max(shadow, 2)), shadow)]


def choose_coding(peaks, selection, min_rel, min_count, rel_window_only):
    """(pos, row, utr_len) of the chosen peak, or None."""
    cands = [p for p in peaks if p[3]]
    if not cands:
        return None
    reference = cands if rel_window_only else peaks
    threshold = max(min_count, min_rel * max(p[1]["count"] for p in reference))
    kept = [c for c in cands if c[1]["count"] >= threshold and c[1]["count"] > 0]
    if not kept:
        return None
    nosl_first = lambda c: 0 if c[1]["tss_type"] == "nosl" else 1
    if selection == "dominant":
        best = min(kept, key=lambda c: (-c[1]["count"], c[2], nosl_first(c)))
    else:
        best = min(kept, key=lambda c: (c[2], nosl_first(c), -c[1]["count"]))
    return best[:3]


def noncoding_candidates(t, gene_id, tss_index, cds_index, max_dist, min_length):
    """(pos, tss_row, shift) for peaks within max_dist of the EGAPx 5' end of a non-coding transcript."""
    exons = t.exon_intervals()
    key = (t.seqid, t.strand)
    total = exonic_length_between(exons, exons[0][0], exons[-1][1])
    out = []
    if t.strand == "+":
        first = exons[0][0]
        for pos, row in tss_index.window(t.seqid, "+", first - max_dist, exons[-1][1]):
            if pos < first:
                shift, guard = first - pos, (pos, first - 1)
            else:
                if containing_interval(exons, pos) is None:
                    continue
                removed = exonic_length_between(exons, first, pos - 1)
                if removed > max_dist or total - removed < min_length:
                    continue
                shift, guard = -removed, (pos, pos)
            if not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, shift))
    else:
        last = exons[-1][1]
        for pos, row in tss_index.window(t.seqid, "-", exons[0][0], last + max_dist):
            if pos > last:
                shift, guard = pos - last, (last + 1, pos)
            else:
                if containing_interval(exons, pos) is None:
                    continue
                removed = exonic_length_between(exons, pos + 1, last)
                if removed > max_dist or total - removed < min_length:
                    continue
                shift, guard = -removed, (pos, pos)
            if not cds_index.overlaps(key, *guard, exclude_owner=gene_id):
                out.append((pos, row, shift))
    return out


def choose_noncoding(cands, min_count):
    kept = [c for c in cands if c[1]["count"] >= min_count and c[1]["count"] > 0]
    if not kept:
        return None
    return min(kept, key=lambda c: (-c[1]["count"], abs(c[2]), 0 if c[1]["tss_type"] == "nosl" else 1))


def set_five_prime_end(t, pos):
    if t.strand == "+":
        kept = [e for e in t.exons if e.end >= pos]
        kept[0].start = pos
    else:
        kept = [e for e in t.exons if e.start <= pos]
        kept[-1].end = pos
    t.exons = kept
    t.sync_bounds()


def transcript_length(t):
    return sum(e.end - e.start + 1 for e in t.exons)


def acceptor_ag(genome, seqid, strand, pos):
    seq = genome.get(seqid)
    if seq is None:
        return None
    if strand == "+":
        return pos >= 3 and seq[pos - 3:pos - 1].upper() == "AG"
    return seq[pos:pos + 2].upper() == "CT"


def record_tss(t, pos, row, shift, genome, stats):
    """Set the TSS attributes of transcript t; add the acceptor row for an sl TSS."""
    f = t.feature
    f.set("tss_evidence", "CAGE")
    f.set("tss_type", row["tss_type"])
    f.set("tss_id", row["tss_id"])
    f.set("tss_peak", f"{t.seqid}:{pos}:{t.strand}")
    f.set("tss_count", f"{row['count']:g}")
    f.set("sl_count", f"{row['sl_count']:g}")
    f.set("sl_fraction", row["sl_fraction"])
    f.set("tss_shift", shift)
    stats["extended" if shift > 0 else "trimmed" if shift < 0 else "unchanged"] += 1
    if row["tss_type"] != "sl":
        return
    stats["sl"] += 1
    if row["sl_leader"] != "NA":
        f.set("sl_leader", row["sl_leader"])
    ag = acceptor_ag(genome, t.seqid, t.strand, pos)
    if ag is not None:
        f.set("acceptor_AG", "true" if ag else "false")
        stats["sl with AG"] += bool(ag)
    site = Feature(t.seqid, "CAGE", TSA, pos, pos, ".", t.strand, ".", collections.OrderedDict())
    site.set("ID", f"tsa-{feature_core_id(t)}")
    site.set("Parent", t.id)
    site.set("tss_id", row["tss_id"])
    site.set("sl_fraction", row["sl_fraction"])
    if row["sl_leader"] != "NA":
        site.set("sl_leader", row["sl_leader"])
    t.others.append(site)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tss", required=True, help="table from cage_tss_table.py")
    parser.add_argument("--fasta", required=True, help="genome, for sequence checks and the AG check of SL acceptors")
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-dist", type=int, default=500,
                        help="maximum 5' UTR length (coding) or 5' end move (non-coding), along the transcript")
    parser.add_argument("--selection", choices=["closest_atg", "dominant"], default="closest_atg")
    parser.add_argument("--min-rel", type=float, default=0.1, help="minimum count relative to the strongest reference peak")
    parser.add_argument("--rel-window-only", action="store_true",
                        help="compare --min-rel with the strongest peak within --max-dist only, not the whole EGAPx 5' UTR")
    parser.add_argument("--min-count", type=float, default=0, help="minimum pooled count of a TSS")
    parser.add_argument("--min-utr-len", type=int, default=0, help="minimum 5' UTR length; 0 accepts a TSS at the ATG")
    parser.add_argument("--nc-min-length", type=int, default=100,
                        help="minimum length left of a non-coding transcript after trimming its 5' end")
    parser.add_argument("--coding-only", action="store_true", help="leave non-coding transcripts unchanged")
    parser.add_argument("--sl-shadow", type=int, default=5,
                        help="nosl peaks this close to a stronger sl peak are ignored (0 to disable)")
    args = parser.parse_args()

    ann = Annotation.read(args.gff)
    genome = dict(read_fasta(args.fasta))
    tss_items = read_tss(args.tss)
    check_sequences(ann, genome, tss_items, args.tss)
    tss_index = PointIndex(tss_items)
    cds_index = ann.cds_index()

    stats = {"mRNA": collections.Counter(), "non-coding": collections.Counter()}
    for gene, t in ann.coding_transcripts():
        st = stats["mRNA"]
        st["total"] += 1
        # drop acceptor rows from an earlier run, so the script can be re-applied
        t.others = [o for o in t.others if o.type != TSA]
        if t.five_prime_partial():
            st["5' partial, kept"] += 1
            t.feature.set("tss_evidence", "EGAPx")
            continue
        cds_before = [(c.start, c.end) for c in t.cds]
        old_len = t.utr_lengths()[0]
        peaks = coding_candidates(t, gene.feature.id, tss_index, cds_index, args.max_dist, args.min_utr_len)
        pick = choose_coding(drop_sl_shadows(peaks, tss_index, t.seqid, t.strand, args.sl_shadow),
                             args.selection, args.min_rel, args.min_count, args.rel_window_only)
        if pick is None:
            st["no CAGE TSS, EGAPx kept"] += 1
            t.feature.set("tss_evidence", "EGAPx")
            continue
        pos, row, utr_len = pick
        set_five_prime_end(t, pos)
        assert [(c.start, c.end) for c in t.cds] == cds_before, f"CDS changed in {t.id}"
        assert t.utr_lengths()[0] == utr_len, f"5' UTR length mismatch in {t.id}"
        st["CAGE TSS"] += 1
        record_tss(t, pos, row, utr_len - old_len, genome, st)

    if not args.coding_only:
        for gene, t in ann.noncoding_transcripts():
            st = stats["non-coding"]
            st["total"] += 1
            t.others = [o for o in t.others if o.type != TSA]
            cands = noncoding_candidates(t, gene.feature.id, tss_index, cds_index, args.max_dist, args.nc_min_length)
            pick = choose_noncoding(drop_sl_shadows(cands, tss_index, t.seqid, t.strand, args.sl_shadow),
                                    args.min_count)
            if pick is None:
                st["no CAGE TSS, EGAPx kept"] += 1
                t.feature.set("tss_evidence", "EGAPx")
                continue
            pos, row, shift = pick
            old_len = transcript_length(t)
            set_five_prime_end(t, pos)
            assert transcript_length(t) - old_len == shift, f"5' end shift mismatch in {t.id}"
            st["CAGE TSS"] += 1
            record_tss(t, pos, row, shift, genome, st)

    ann.write(args.output)
    for kind, counter in stats.items():
        for key, value in counter.items():
            print(f"{kind}\t{key}\t{value}", file=sys.stderr)


if __name__ == "__main__":
    main()
