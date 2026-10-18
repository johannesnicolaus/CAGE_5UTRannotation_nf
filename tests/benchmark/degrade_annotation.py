#!/usr/bin/env python3
"""Make an EGAPx-like test input from a reference annotation with trusted UTRs.

  degrade_annotation.py --gff reference.gff3 --mode strip|extend [--max-extend 500]
                        [--seed 1] --output degraded.gff3

Only coding mRNAs (transcripts with a CDS) are changed; CDS rows never change.
  strip   every model starts at its start codon and ends at its stop codon (no UTRs):
          tests whether the pipeline adds the right UTRs.
  extend  the 5' and 3' ends move outward by a random 0..--max-extend bases each, like the
          over-extended UTRs of EGAPx models: tests whether the pipeline trims them back.
Models open at an end (start_range / end_range) are left alone at that end. Attributes
that carry the reference's evidence (e.g. 'experiment' poly(A) evidence) are removed, so
the pipeline cannot keep reference ends through them.
"""

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from gff3_utils import Annotation  # noqa: E402


def set_ends(t, new5, new3):
    """Move the transcript's 5' and 3' ends to genomic positions new5 / new3."""
    exons = sorted(t.exons, key=lambda e: e.start)
    if t.strand == "+":
        lo, hi = new5, new3
    else:
        lo, hi = new3, new5
    # exons overlapping the new span; the outer ones are cut or stretched to it. New ends are
    # either exonic (strip: start/stop codons) or outside the model (extend).
    kept = [e for e in exons if e.end >= lo and e.start <= hi]
    kept[0].start = lo
    kept[-1].end = hi
    t.exons = kept
    t.sync_bounds()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--mode", choices=["strip", "extend"], required=True)
    parser.add_argument("--max-extend", type=int, default=500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    ann = Annotation.read(args.gff)
    n = 0
    for _gene, t in ann.coding_transcripts():
        t.feature.attrs.pop("experiment", None)
        for e in t.exons:
            e.attrs.pop("experiment", None)
        cds_start, cds_end = t.cds_bounds()
        cur5 = t.feature.start if t.strand == "+" else t.feature.end
        cur3 = t.feature.end if t.strand == "+" else t.feature.start
        atg = cds_start if t.strand == "+" else cds_end
        stop = cds_end if t.strand == "+" else cds_start
        if args.mode == "strip":
            new5 = cur5 if t.five_prime_partial() else atg
            new3 = cur3 if t.three_prime_partial() else stop
        else:
            d5 = 0 if t.five_prime_partial() else rng.randint(0, args.max_extend)
            d3 = 0 if t.three_prime_partial() else rng.randint(0, args.max_extend)
            new5 = max(1, cur5 - d5) if t.strand == "+" else cur5 + d5
            new3 = cur3 + d3 if t.strand == "+" else max(1, cur3 - d3)
        set_ends(t, new5, new3)
        n += 1
    ann.write(args.output)
    print(f"{n} coding mRNAs degraded ({args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
