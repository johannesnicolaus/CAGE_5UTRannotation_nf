#!/usr/bin/env python3
"""Turn CAGE TSS clusters into one TSS table with a spliced-leader label per cluster.

Two input modes; this is the only script that knows the input formats.

cagescan mode (luscombeu/cagescan results, tss_clusters/<set>/):
  cage_tss_table.py --clusters all_samples.<set>.tss_clusters.bed \\
                    --counts   all_samples.<set>.tss_cluster_counts.tsv \\
                    [--sl-counts  all_samples.<set>.tss_cluster_sl_counts.tsv] \\
                    [--sl-leaders all_samples.<set>.tss_cluster_sl_leaders.tsv] \\
                    --output PREFIX.tss.tsv
  tss_clusters.bed       BED9; name = tss_id, thickStart = dominant CTSS (0-based)
  tss_cluster_counts.tsv tss_id, location, strand, peak, <one column per sample>
  tss_cluster_sl_counts  same layout, READ1s carrying the spliced leader; cagescan writes
                         it only when the dataset has SL reads. With cagescan's separate
                         SL clusters (default), 'SL...' clusters hold only SL reads and
                         'TC...' clusters none, so their fractions are 1 and 0.
  tss_cluster_sl_leaders main leader of each cluster with SL reads (several leaders only)
  A cluster is 'sl' when sl_count / count >= --sl-min-fraction.

BED mode (TSS clusters from any CAGE workflow, single- or paired-end, e.g.
nf-core/cageseq, CAGEr, CAGEfightR):
  cage_tss_table.py --clusters tss_clusters.bed [--sl-bed sl_acceptor_clusters.bed] \\
                    --output PREFIX.tss.tsv
  BED6 or longer: chrom, start, end, name, score, strand[, thickStart, thickEnd, ...].
  - peak: thickStart when the BED has >= 8 columns and thickStart < thickEnd (the
    dominant CTSS, as written by CAGEr/CAGEfightR/cagescan); otherwise the cluster must
    be 1 bp wide (a CTSS), or its 5' end is used with a warning;
  - count: the score column (reads or TPM); clusters with score 0 or '.' get count 0;
  - name: the cluster ID ('.' -> generated from the position).
  Clusters from --clusters are 'nosl'; clusters from --sl-bed (SL trans-splice acceptor
  clusters, from SL reads) are 'sl'. Without --sl-bed every TSS is 'nosl'.

Output columns:
  chrom, peak (1-based), strand, tss_id, count, sl_count, sl_fraction, tss_type, sl_leader
count and sl_count are summed over samples. sl_leader is the main leader of an sl TSS
when several leaders were searched for, otherwise NA.
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


def read_bed(path):
    """Yield (chrom, peak 1-based, strand, name, score, fields) from a BED6+ file."""
    warned = 0
    with open(path) as bed:
        for line in bed:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[5] not in ("+", "-"):
                sys.exit(f"ERROR: {path} needs BED6 or longer with a strand (+/-): {line!r}")
            chrom, start, end, strand = f[0], int(f[1]), int(f[2]), f[5]
            if len(f) >= 8 and f[6].isdigit() and f[7].isdigit() and int(f[6]) < int(f[7]):
                peak = int(f[6]) + 1
            else:
                if end - start > 1:
                    warned += 1
                peak = start + 1 if strand == "+" else end
            name = f[3] if f[3] not in ("", ".") else f"{chrom}:{peak}:{strand}"
            score = float(f[4]) if f[4] not in ("", ".") else 0.0
            yield chrom, peak, strand, name, score, f
    if warned:
        print(f"WARNING: {warned} clusters in {path} are wider than 1 bp and have no thickStart; "
              "their 5' end was used as the TSS", file=sys.stderr)


HEADER = "chrom\tpeak\tstrand\ttss_id\tcount\tsl_count\tsl_fraction\ttss_type\tsl_leader\n"


def cagescan_mode(args, out):
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
    for chrom, peak, strand, tss_id, _score, f in read_bed(args.clusters):
        if len(f) < 8:
            sys.exit(f"ERROR: {args.clusters} is not BED9 (thickStart needed)")
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
    return n, n_sl, source


def bed_mode(args, out):
    n = n_sl = 0
    seen = set()
    for path, tss_type in ((args.clusters, "nosl"), (args.sl_bed, "sl")):
        if not path:
            continue
        for chrom, peak, strand, tss_id, score, _f in read_bed(path):
            if tss_id in seen:
                tss_id = f"{tss_id}_{tss_type}"
            seen.add(tss_id)
            sl = score if tss_type == "sl" else 0.0
            fraction = 1.0 if tss_type == "sl" else 0.0
            n += 1
            n_sl += tss_type == "sl"
            out.write(f"{chrom}\t{peak}\t{strand}\t{tss_id}\t{score:g}\t{sl:g}\t{fraction:.3f}\t{tss_type}\tNA\n")
    source = args.sl_bed if args.sl_bed else "no SL clusters: all nosl"
    return n, n_sl, source


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clusters", required=True, help="TSS clusters (BED); non-SL clusters in BED mode")
    parser.add_argument("--counts", default=None, help="cagescan count matrix (selects cagescan mode)")
    parser.add_argument("--sl-counts", default=None)
    parser.add_argument("--sl-leaders", default=None)
    parser.add_argument("--sl-bed", default=None, help="BED mode: SL trans-splice acceptor clusters")
    parser.add_argument("--sl-min-fraction", type=float, default=0.5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.counts and args.sl_bed:
        sys.exit("ERROR: --sl-bed is for BED mode; cagescan results carry SL counts instead")

    with open(args.output, "w") as out:
        out.write(HEADER)
        n, n_sl, source = cagescan_mode(args, out) if args.counts else bed_mode(args, out)
    print(f"{n} TSS clusters, {n_sl} labelled sl ({source})", file=sys.stderr)


if __name__ == "__main__":
    main()
