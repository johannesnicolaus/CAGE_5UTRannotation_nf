#!/usr/bin/env python3
"""Score the 5' and 3' ends of a pipeline output against a reference with trusted UTRs.

  compare_utrs.py --reference reference.gff3 --input degraded.gff3 --output pipeline.utr.gff3 \\
                  [--pas pas.bed] --prefix PREFIX

Coding mRNAs are matched by transcript ID. Only reference ends that are known are scored:
ends of models open at that end (NCBI/EGAPx start_range / end_range, e.g. RefSeq C. elegans
models without UTR evidence) are skipped. For each end, the signed error is the distance
(nt, along the genome, positive = the end lies downstream in transcription direction) between
the end and the reference end. Writes:
  PREFIX.per_transcript.tsv  errors of the input and of the output, with the evidence class
  PREFIX.summary.tsv         share of ends within 0 / 5 / 20 / 50 nt, median absolute error,
                             per evidence class, for the input (baseline) and the output
  PREFIX.pas_offset.tsv      (with --pas) for PAS-set 3' ends: distance from the motif end to
                             the reference 3' end, i.e. the offset that would have been exact
"""

import argparse
import collections
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from gff3_utils import Annotation  # noqa: E402


def ends(path):
    """transcript ID -> (strand, 5' end, 3' end, attributes, 5' known, 3' known) for coding mRNAs."""
    out = {}
    for _g, t in Annotation.read(path).coding_transcripts():
        f = t.feature
        five, three = (f.start, f.end) if t.strand == "+" else (f.end, f.start)
        out[t.id] = (t.strand, five, three, f, not t.five_prime_partial(), not t.three_prime_partial())
    return out


def signed(strand, pos, ref):
    return (pos - ref) if strand == "+" else (ref - pos)


def summarise(errors):
    a = [abs(e) for e in errors]
    n = len(a)
    if not n:
        return [0, "NA", "NA", "NA", "NA", "NA"]
    return [n] + [f"{100 * sum(x <= k for x in a) / n:.1f}" for k in (0, 5, 20, 50)] + [f"{statistics.median(a):g}"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--input", required=True, help="the degraded annotation given to the pipeline")
    parser.add_argument("--output", required=True, help="the pipeline's <id>.utr.gff3")
    parser.add_argument("--pas", help="the pipeline's <id>.pas.bed, for the offset analysis")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    ref, inp, out = ends(args.reference), ends(args.input), ends(args.output)
    rows, groups = [], collections.defaultdict(list)
    for tid, (strand, r5, r3, _, known5, known3) in ref.items():
        if tid not in out or tid not in inp:
            continue
        _, i5, i3, _, _, _ = inp[tid]
        _, o5, o3, f, _, _ = out[tid]
        ev5 = f.get("tss_evidence") or "NA"
        tt = f.get("tss_type") or "NA"
        ev3 = f.get("utr3_evidence") or "NA"
        e = [signed(strand, i5, r5) if known5 else "NA", signed(strand, o5, r5) if known5 else "NA",
             signed(strand, i3, r3) if known3 else "NA", signed(strand, o3, r3) if known3 else "NA"]
        rows.append([tid, strand, ev5, tt, ev3, f.get("pas_pos") or "NA"] + e)
        if known5:
            for cls in ("all", f"{ev5}/{tt}"):
                groups[("5'", "input", cls)].append(e[0])
                groups[("5'", "output", cls)].append(e[1])
        if known3:
            for cls in ("all", ev3):
                groups[("3'", "input", cls)].append(e[2])
                groups[("3'", "output", cls)].append(e[3])

    with open(f"{args.prefix}.per_transcript.tsv", "w") as fh:
        fh.write("transcript_id\tstrand\ttss_evidence\ttss_type\tutr3_evidence\tpas_pos\t"
                 "err5_input\terr5_output\terr3_input\terr3_output\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(f"{args.prefix}.summary.tsv", "w") as fh:
        fh.write("end\tannotation\tclass\tn\twithin_0\twithin_5\twithin_20\twithin_50\tmedian_abs_error\n")
        for key in sorted(groups):
            fh.write("\t".join(list(key) + [str(x) for x in summarise(groups[key])]) + "\n")

    if args.pas:
        # PAS-set ends: where would the exact end have been relative to the motif end?
        offs = []
        for r in rows:
            if r[4] != "PAS" or r[5] == "NA" or r[9] == "NA":
                continue
            tid, strand = r[0], r[1]
            span = r[5].split(":")[1]
            ms, me = map(int, span.split("-"))
            motif_end = me if strand == "+" else ms
            offs.append(signed(strand, ref[tid][2], motif_end))
        c = collections.Counter(offs)
        with open(f"{args.prefix}.pas_offset.tsv", "w") as fh:
            fh.write("offset_motif_end_to_reference_end\tn\n")
            for k in sorted(c):
                fh.write(f"{k}\t{c[k]}\n")
        if offs:
            near = [o for o in offs if -10 <= o <= 60]
            print(f"PAS-set 3' ends: {len(offs)}; reference end 0-40 nt after the motif: "
                  f"{100 * sum(0 <= o <= 40 for o in offs) / len(offs):.1f}%; median offset among -10..60: "
                  f"{statistics.median(near) if near else 'NA'}", file=sys.stderr)
    print(f"{len(rows)} coding mRNAs compared", file=sys.stderr)


if __name__ == "__main__":
    main()
