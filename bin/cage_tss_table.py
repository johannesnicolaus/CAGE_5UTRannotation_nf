#!/usr/bin/env python3
"""Turn cagescan TSS clusters into one TSS table with a spliced-leader label per cluster.

  cage_tss_table.py --clusters all_samples.<set>.tss_clusters.bed \\
                    --counts   all_samples.<set>.tss_cluster_counts.tsv \\
                    [--sl-counts all_samples.<set>.tss_cluster_sl_counts.tsv] \\
                    [--sl-leaders all_samples.<set>.tss_cluster_sl_leaders.tsv] \\
                    --sl-min-fraction 0.5 --output PREFIX.tss.tsv

Inputs (luscombeu/cagescan, tss_clusters/<set>/):
  tss_clusters.bed       BED9; name = tss_id, thickStart = dominant CTSS (0-based)
  tss_cluster_counts.tsv tss_id, location, strand, peak, <one column per sample>
  tss_cluster_sl_counts  same layout, READ1s carrying the spliced leader; cagescan writes
                         it only when the dataset has SL reads. With cagescan's separate
                         SL clusters (default), 'SL...' clusters hold only SL reads and
                         'TC...' clusters none, so their fractions are 1 and 0.
  tss_cluster_sl_leaders main leader of each cluster with SL reads (several leaders only)

Output columns:
  chrom, peak (1-based), strand, tss_id, count, sl_count, sl_fraction, tss_type, sl_leader
count and sl_count are summed over samples; tss_type is 'sl' when
sl_count / count >= --sl-min-fraction, otherwise 'nosl'. Without an SL table every
TSS is 'nosl'. sl_leader is the main leader of an sl TSS when several leaders were
searched for, otherwise NA.

This is the only script that knows the cagescan output format.
"""

import argparse
import sys


def read_matrix(path):
    """tss_id -> sum of the sample columns."""
    totals = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:4] != ["tss_id", "location", "strand", "peak"]:
            sys.exit(f"ERROR: unexpected header in {path}: {header[:4]}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            totals[fields[0]] = sum(float(x) for x in fields[4:] if x not in ("", "NA"))
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--sl-counts", default=None)
    parser.add_argument("--sl-leaders", default=None)
    parser.add_argument("--sl-min-fraction", type=float, default=0.5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = read_matrix(args.counts)
    sl_counts = read_matrix(args.sl_counts) if args.sl_counts else {}
    leaders = {}
    if args.sl_leaders:
        with open(args.sl_leaders) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            if "main_leader" not in header:
                sys.exit(f"ERROR: no main_leader column in {args.sl_leaders}")
            i = header.index("main_leader")
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                leaders[fields[0]] = fields[i]

    n = n_sl = n_missing = 0
    with open(args.clusters) as bed, open(args.output, "w") as out:
        out.write("chrom\tpeak\tstrand\ttss_id\tcount\tsl_count\tsl_fraction\ttss_type\tsl_leader\n")
        for line in bed:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                sys.exit(f"ERROR: {args.clusters} is not BED9 (thickStart needed): {line!r}")
            chrom, tss_id, strand = fields[0], fields[3], fields[5]
            peak = int(fields[6]) + 1
            if tss_id not in counts:
                n_missing += 1
            count = counts.get(tss_id, 0)
            sl = sl_counts.get(tss_id, 0)
            fraction = sl / count if count > 0 else 0.0
            tss_type = "sl" if sl > 0 and fraction >= args.sl_min_fraction else "nosl"
            n += 1
            n_sl += tss_type == "sl"
            leader = leaders.get(tss_id, "NA") if tss_type == "sl" else "NA"
            out.write(f"{chrom}\t{peak}\t{strand}\t{tss_id}\t{count:g}\t{sl:g}\t{fraction:.3f}\t{tss_type}\t{leader}\n")

    if n_missing:
        print(f"WARNING: {n_missing} clusters have no row in {args.counts}; their count is 0", file=sys.stderr)
    source = args.sl_counts if args.sl_counts else "no SL table: all nosl"
    print(f"{n} TSS clusters, {n_sl} labelled sl ({source})", file=sys.stderr)


if __name__ == "__main__":
    main()
