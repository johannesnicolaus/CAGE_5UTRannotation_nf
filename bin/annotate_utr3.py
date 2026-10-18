#!/usr/bin/env python3
"""Refine the 3' end of each coding mRNA with polyadenylation evidence.

  annotate_utr3.py --gff PREFIX.utr5.gff3 --pas PREFIX.pas.bed --output PREFIX.utr.gff3 [options]

For every mRNA with a CDS (models open at the 3' end are left alone), the first rule
that applies sets the 3' end:

1. polyA_reads  EGAPx marked the model 'experiment=COORDINATES: polyA evidence': reads
                with the poly(A) tail mapped the cleavage site. The EGAPx end is kept
                (disable with --no-keep-polya-evidence).
2. PAS          a polyadenylation signal (scan_pas.py) on the same strand
                - entirely downstream of the stop codon,
                - within --window bases of the EGAPx 3' end, measured along the transcript
                  (3' UTR exons inside the model, genomic sequence beyond it),
                - not in a 3' UTR intron, not reaching a CDS of another gene on the same
                  strand, and not extending the transcript into another gene's transcript
                  on the same strand (in operons the 3' end stops before the next gene's
                  5' end, which ANNOTATE_UTR5 has already set; overlaps already present
                  in EGAPx are not counted).
                The motif closest to the EGAPx end is used, and the new 3' end is the
                motif end + --offset bases (cleavage lies ~10-30 nt downstream of AATAAA).
                The 3' UTR is trimmed or extended.
3. EGAPx        no PAS qualified; the EGAPx 3' end is kept.

Attributes added to the mRNA (the three_prime_UTR rows carry only ID and Parent):
  utr3_evidence (polyA_reads, PAS, EGAPx), pas_motif, pas_pos (chrom:start-end:strand of
  the motif), utr3_shift (new minus old 3' UTR length).
"""

import argparse
import collections
import sys

from gff3_utils import Annotation, PointIndex, containing_interval, exonic_length_between

POLYA_EVIDENCE = "polyA evidence"


def read_pas(path):
    items = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            start, end = int(f[1]) + 1, int(f[2])
            items.append((f[0], f[5], start, (start, end, f[3])))
    return PointIndex(items)


def pas_candidates(t, gene_id, pas_index, cds_index, tx_index, window, offset):
    """(new_end, signed shift along the transcript, motif) for every qualifying PAS."""
    exons = t.exon_intervals()
    cds_start, cds_end = t.cds_bounds()
    key = (t.seqid, t.strand)
    out = []
    if t.strand == "+":
        end0 = exons[-1][1]
        limit = end0 + window
        blockers = [s for s, _, _ in cds_index.overlaps(key, end0 + 1, limit, exclude_owner=gene_id)]
        blockers += [s for s, _, _ in tx_index.overlaps(key, end0 + 1, limit, exclude_owner=gene_id) if s > end0]
        if blockers:
            limit = min(blockers) - 1
        for ms, (_, me, motif) in pas_index.window(t.seqid, "+", cds_end + 1, limit):
            new_end = me + offset
            i = containing_interval(exons, ms)
            if i is None:
                if ms <= end0:
                    continue  # in a 3' UTR intron
                segment_end = limit
            else:
                segment_end = limit if i == len(exons) - 1 else exons[i][1]
            if me > segment_end or new_end > segment_end:
                continue
            shift = new_end - end0 if new_end >= end0 else -exonic_length_between(exons, new_end + 1, end0)
            if abs(shift) <= window:
                out.append((new_end, shift, (ms, me, motif)))
    else:
        end0 = exons[0][0]
        limit = end0 - window
        blockers = [e for _, e, _ in cds_index.overlaps(key, limit, end0 - 1, exclude_owner=gene_id)]
        blockers += [e for _, e, _ in tx_index.overlaps(key, limit, end0 - 1, exclude_owner=gene_id) if e < end0]
        if blockers:
            limit = max(blockers) + 1
        for ms, (_, me, motif) in pas_index.window(t.seqid, "-", limit, cds_start - 1):
            if me >= cds_start:
                continue
            new_end = ms - offset
            i = containing_interval(exons, me)
            if i is None:
                if me >= end0:
                    continue
                segment_start = limit
            else:
                segment_start = limit if i == 0 else exons[i][0]
            if ms < segment_start or new_end < max(segment_start, 1):
                continue
            shift = end0 - new_end if new_end <= end0 else -exonic_length_between(exons, end0, new_end - 1)
            if abs(shift) <= window:
                out.append((new_end, shift, (ms, me, motif)))
    return out


def set_three_prime_end(t, pos):
    if t.strand == "+":
        kept = [e for e in t.exons if e.start <= pos]
        kept[-1].end = pos
    else:
        kept = [e for e in t.exons if e.end >= pos]
        kept[0].start = pos
    t.exons = kept
    t.sync_bounds()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--pas", required=True, help="BED6 from scan_pas.py")
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=500, help="maximum distance of the PAS from the EGAPx 3' end")
    parser.add_argument("--offset", type=int, default=0, help="bases from the motif end to the 3' end")
    parser.add_argument("--no-keep-polya-evidence", dest="keep_polya", action="store_false",
                        help="also move 3' ends that EGAPx supports with poly(A) reads")
    args = parser.parse_args()

    ann = Annotation.read(args.gff)
    pas_index = read_pas(args.pas)
    cds_index = ann.cds_index()
    tx_index = ann.transcript_index()

    stats = collections.Counter()
    for gene, t in ann.coding_transcripts():
        stats["mRNA"] += 1
        f = t.feature
        if t.three_prime_partial():
            stats["3' partial, kept"] += 1
            f.set("utr3_evidence", "EGAPx")
            continue
        if args.keep_polya and POLYA_EVIDENCE in (f.get("experiment") or ""):
            stats["polyA_reads"] += 1
            f.set("utr3_evidence", "polyA_reads")
            continue
        cands = pas_candidates(t, gene.feature.id, pas_index, cds_index, tx_index, args.window, args.offset)
        if not cands:
            stats["EGAPx"] += 1
            f.set("utr3_evidence", "EGAPx")
            continue

        new_end, shift, (ms, me, motif) = min(cands, key=lambda c: (abs(c[1]), c[1]))
        cds_before = [(c.start, c.end) for c in t.cds]
        old_len = t.utr_lengths()[1]
        set_three_prime_end(t, new_end)
        assert [(c.start, c.end) for c in t.cds] == cds_before, f"CDS changed in {t.id}"
        assert t.utr_lengths()[1] - old_len == shift, f"3' UTR length mismatch in {t.id}"
        stats["PAS"] += 1
        stats["PAS extended" if shift > 0 else "PAS trimmed" if shift < 0 else "PAS unchanged"] += 1
        f.set("utr3_evidence", "PAS")
        f.set("pas_motif", motif)
        f.set("pas_pos", f"{t.seqid}:{ms}-{me}:{t.strand}")
        f.set("utr3_shift", shift)

    ann.write(args.output)
    for key, value in stats.items():
        print(f"{key}\t{value}", file=sys.stderr)


if __name__ == "__main__":
    main()
