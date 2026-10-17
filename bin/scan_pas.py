#!/usr/bin/env python3
"""Find polyadenylation-signal motifs (default AATAAA) on both strands of a genome.

  scan_pas.py --fasta genome.fa --motifs AATAAA[,ATTAAA] --output PREFIX.pas.bed

Output: BED6 (0-based start), name = motif, score = 0, strand = strand of the RNA that
would carry the motif. Overlapping occurrences are all reported. Case is ignored, so
soft-masked genomes are scanned completely. Port of scan_motif() in annotationpolish.
"""

import argparse
import re

from gff3_utils import read_fasta

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--motifs", default="AATAAA", help="comma-separated DNA motifs (RNA sense)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    motifs = [m.strip().upper().replace("U", "T") for m in args.motifs.split(",") if m.strip()]
    for m in motifs:
        if not re.fullmatch(r"[ACGT]+", m):
            raise SystemExit(f"ERROR: motif {m!r} must contain only A, C, G, T")
    # lookahead patterns report overlapping matches
    patterns = [(m, re.compile(f"(?={m})"), re.compile(f"(?={revcomp(m)})")) for m in motifs]

    with open(args.output, "w") as out:
        for name, seq in read_fasta(args.fasta):
            seq = seq.upper()
            hits = []
            for motif, fwd, rev in patterns:
                k = len(motif)
                hits += [(m.start(), m.start() + k, motif, "+") for m in fwd.finditer(seq)]
                hits += [(m.start(), m.start() + k, motif, "-") for m in rev.finditer(seq)]
            hits.sort()
            for start, end, motif, strand in hits:
                out.write(f"{name}\t{start}\t{end}\t{motif}\t0\t{strand}\n")


if __name__ == "__main__":
    main()
