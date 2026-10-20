#!/usr/bin/env python3
"""Build the test dataset (-profile test) from public C. elegans data.

  make_testdata.py --genome GCF_000002985.6_WBcel235_genomic.fna.gz \\
                   --annotation ce_extend.gff3 \\
                   --cagescan <cagescan results dir> --outdir tests/testdata

Sources (all public):
  genome      WBcel235, NCBI RefSeq GCF_000002985.6
  annotation  WormBase WS298 gene set in RefSeq format (GCF_000002985.6), with UTRs pushed
              outward by 0-500 nt (tests/benchmark/degrade_annotation.py --mode extend --seed 1),
              to mimic over-extended EGAPx models
  CAGE        SLIC-CAGE of total RNA, young-adult hermaphrodites (SRR10568191, PRJNA592867),
              processed with luscombeu/cagescan (main@80098bc) with SL1 and SL2 leaders

The test region is chromosome III (NC_003281.10) 4,800,001-5,200,000: 144 genes, 18 WormBase
operons, 124 WormBase SL2 acceptor sites. Coordinates are shifted so that the region starts
at 1, on a sequence named 'III_4.8M'. Only genes lying entirely in the region are kept.

Outputs (in --outdir):
  genome.fa.gz, annotation.gff3.gz
  cagescan/tss_clusters/broad/all_samples.broad.*   (cagescan input mode)
  tss_clusters.bed, tss_sl_clusters.bed              (BED input mode)
  samplesheet.csv                                    (one sample per mode; relative paths)
"""

import argparse
import gzip
import os

CHROM, START, END, NAME = "NC_003281.10", 4_800_001, 5_200_000, "III_4.8M"
OFF = START - 1


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def inside(a, b):
    return START <= a and b <= END


def write_genome(src, dst):
    seq, keep = [], False
    with opener(src) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == CHROM
            elif keep:
                seq.append(line.strip())
    s = "".join(seq)[START - 1:END]
    with gzip.open(dst, "wt") as out:
        out.write(f">{NAME} {CHROM}:{START}-{END} (WBcel235)\n")
        for i in range(0, len(s), 80):
            out.write(s[i:i + 80] + "\n")


def write_annotation(src, dst):
    rows, parents = [], {}
    with opener(src) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != CHROM:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append((f, attrs))
            if "ID" in attrs:
                parents.setdefault(attrs["ID"], attrs.get("Parent"))

    def top(i):
        while parents.get(i):
            i = parents[i]
        return i

    spans = {}
    for f, attrs in rows:
        if "Parent" not in attrs and "ID" in attrs:
            a, b = spans.get(attrs["ID"], (10**12, 0))
            spans[attrs["ID"]] = (min(a, int(f[3])), max(b, int(f[4])))
    keep = {g for g, (a, b) in spans.items() if inside(a, b)}
    n = 0
    with gzip.open(dst, "wt") as out:
        out.write("##gff-version 3\n")
        out.write(f"##sequence-region {NAME} 1 {END - OFF}\n")
        for f, attrs in rows:
            gid = top(attrs.get("ID") or attrs.get("Parent", "").split(",")[0])
            if gid not in keep:
                continue
            f = list(f)
            f[0] = NAME
            f[3], f[4] = str(int(f[3]) - OFF), str(int(f[4]) - OFF)
            f[8] = ";".join(shift_range(kv) for kv in f[8].split(";"))
            out.write("\t".join(f) + "\n")
            n += 1
    return len(keep), n


def shift_range(kv):
    if kv.startswith(("start_range=", "end_range=")):
        key, val = kv.split("=", 1)
        parts = [p if p == "." else str(int(p) - OFF) for p in val.split(",")]
        return f"{key}={','.join(parts)}"
    return kv


def write_cage(cagescan, outdir):
    src = os.path.join(cagescan, "tss_clusters", "broad")
    dst = os.path.join(outdir, "cagescan", "tss_clusters", "broad")
    os.makedirs(dst, exist_ok=True)
    keep = set()
    bed_rows = []
    with open(os.path.join(src, "all_samples.broad.tss_clusters.bed")) as fh, \
            open(os.path.join(dst, "all_samples.broad.tss_clusters.bed"), "w") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != CHROM or not inside(int(f[1]) + 1, int(f[2])):
                continue
            keep.add(f[3])
            f[0] = NAME
            for i in (1, 2, 6, 7):
                f[i] = str(int(f[i]) - OFF)
            out.write("\t".join(f) + "\n")
            bed_rows.append(f)
    totals = {}
    for name in ("tss_cluster_counts", "tss_cluster_sl_counts", "tss_cluster_sl_leaders"):
        path = os.path.join(src, f"all_samples.broad.{name}.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as fh, open(os.path.join(dst, f"all_samples.broad.{name}.tsv"), "w") as out:
            out.write(fh.readline())
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] not in keep:
                    continue
                chrom, rng = f[1].rsplit(":", 1)
                a, b = rng.split("-")
                f[1] = f"{NAME}:{int(a) - OFF}-{int(b) - OFF}"
                f[3] = str(int(f[3]) - OFF)
                out.write("\t".join(f) + "\n")
                if name == "tss_cluster_counts":
                    totals[f[0]] = sum(float(x) for x in f[4:])
    # BED input mode: non-SL clusters and SL acceptor clusters, score = pooled READ1 count
    with open(os.path.join(outdir, "tss_clusters.bed"), "w") as nosl, \
            open(os.path.join(outdir, "tss_sl_clusters.bed"), "w") as sl:
        for f in bed_rows:
            row = f[:]
            row[4] = f"{totals.get(f[3], 0):g}"
            (sl if f[3].startswith("SL") else nosl).write("\t".join(row) + "\n")
    return len(keep)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--cagescan", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    write_genome(a.genome, os.path.join(a.outdir, "genome.fa.gz"))
    genes, rows = write_annotation(a.annotation, os.path.join(a.outdir, "annotation.gff3.gz"))
    clusters = write_cage(a.cagescan, a.outdir)
    with open(os.path.join(a.outdir, "samplesheet.csv"), "w") as out:
        out.write("id,genome,annotation,tss,tss_sl\n")
        out.write("ce_cagescan,genome.fa.gz,annotation.gff3.gz,cagescan,\n")
        out.write("ce_bed,genome.fa.gz,annotation.gff3.gz,tss_clusters.bed,tss_sl_clusters.bed\n")
    print(f"{genes} genes ({rows} rows), {clusters} TSS clusters")


if __name__ == "__main__":
    main()
