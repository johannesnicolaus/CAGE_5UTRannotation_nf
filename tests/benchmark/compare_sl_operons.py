#!/usr/bin/env python3
"""Score SL TSS labels and operons against reference SL acceptor sites and operons.

  compare_sl_operons.py --utr pipeline.utr.gff3 --operons pipeline.operons.gff3 \\
      --ref-sl SL_acceptor_sites.gff3 --ref-operons operons.gff3 \\
      [--chrom-map assembly_report.txt] [--gene-xref WormBase] --prefix PREFIX

Written for C. elegans (WormBase SL1_acceptor_site / SL2_acceptor_site features and CEOP
operons with 'genes=WBGene...' lists) on an NCBI-style annotation whose genes carry
'Dbxref=...,WormBase:WBGene...'. --chrom-map translates the reference sequence names
(I, II, ...) to the annotation's (NC_...), using an NCBI assembly report.

Reports (stdout and PREFIX.sl_operons.tsv):
  SL labels  for transcripts whose CAGE TSS is sl: distance to the nearest reference acceptor
             site on the same strand (offset distribution), share within 0 / 2 / 5 nt, and
             agreement of sl_leader (SL1 vs SL2*) with the type of that site; for nosl CAGE
             TSSs: share with a reference acceptor site within 5 nt (possibly missed SL).
  Operons    adjacent co-operonic gene pairs (protein-coding genes known to both): shared,
             only in the pipeline, only in the reference; precision and recall. And the share
             of genes with an SL2-labelled TSS that are downstream genes of a pipeline operon,
             against SL1-labelled genes.
"""

import argparse
import bisect
import collections
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from gff3_utils import Annotation  # noqa: E402


def chrom_map(path):
    m = {}
    if not path:
        return m
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 7:
            m[f[0]] = f[6]
    return m


def read_sites(path, cmap):
    """(seqid, strand) -> sorted list of (pos, kinds). pos: the first transcript base is taken
    as the end of the 2-bp feature on +, its start on - (reported offsets show the convention).
    A position can carry both an SL1 and an SL2 site."""
    sites = collections.defaultdict(lambda: collections.defaultdict(set))
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or "acceptor_site" not in f[2]:
            continue
        seqid = cmap.get(f[0], f[0])
        pos = int(f[4]) if f[6] == "+" else int(f[3])
        sites[(seqid, f[6])][pos].add("SL2" if f[2].startswith("SL2") else "SL1")
    return {k: sorted(v.items()) for k, v in sites.items()}


def nearest(sites, seqid, strand, pos):
    lst = sites.get((seqid, strand))
    if not lst:
        return None
    keys = [p for p, _ in lst]
    i = bisect.bisect_left(keys, pos)
    best = None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(lst):
            d = pos - lst[j][0] if strand == "+" else lst[j][0] - pos  # >0: TSS downstream of the site
            if best is None or abs(d) < abs(best[0]):
                best = (d, lst[j][1])
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--utr", required=True)
    parser.add_argument("--operons", required=True)
    parser.add_argument("--ref-sl", required=True)
    parser.add_argument("--ref-operons", required=True)
    parser.add_argument("--chrom-map")
    parser.add_argument("--gene-xref", default="WormBase", help="Dbxref prefix that links genes to the reference")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    cmap = chrom_map(args.chrom_map)
    sites = read_sites(args.ref_sl, cmap)
    ann = Annotation.read(args.utr)

    # ---- SL labels ----
    off = collections.Counter()
    agree = collections.Counter()
    nosl_near = n_nosl = 0
    gene_leader = {}
    gene_xref = {}
    for gene, t in ann.transcripts():
        dbx = gene.feature.get("Dbxref") or ""
        m = re.search(args.gene_xref + r":([^,;]+)", dbx)
        if m:
            gene_xref[gene.feature.id] = m.group(1)
        f = t.feature
        if f.get("tss_evidence") != "CAGE" or not t.is_coding:
            continue
        seqid, pos, strand = f.get("tss_peak").split(":")
        hit = nearest(sites, seqid, strand, int(pos))
        if f.get("tss_type") == "sl":
            if hit is None:
                off["none"] += 1
                continue
            d, kind = hit
            off[d if abs(d) <= 10 else ("<-10" if d < 0 else ">10")] += 1
            leader = f.get("sl_leader") or "NA"
            if abs(d) <= 2 and leader != "NA":
                ref_kind = "SL1+SL2" if len(kind) == 2 else next(iter(kind))
                agree[(leader[:3], ref_kind)] += 1
            if leader != "NA":
                gene_leader.setdefault(gene.feature.id, set()).add(leader[:3])
        else:
            n_nosl += 1
            if hit is not None and abs(hit[0]) <= 5:
                nosl_near += 1
    n_sl = sum(off.values())
    within = lambda k: sum(v for d, v in off.items() if isinstance(d, int) and abs(d) <= k)

    # ---- operons ----
    def pairs_from(gene_lists):
        out = set()
        for genes in gene_lists:
            for a, b in zip(genes, genes[1:]):
                out.add(frozenset((a, b)))
        return out

    ours_lists, ours_pos = [], {}
    for line in open(args.operons):
        if line.startswith("#") or not line.strip():
            continue
        g = re.search(r";genes=([^;\n]+)", line).group(1).split(",")
        ours_lists.append([gene_xref.get(x, x) for x in g])
        for i, x in enumerate(g):
            ours_pos[x] = i + 1
    ref_lists = []
    for line in open(args.ref_operons):
        if line.startswith("#") or not line.strip():
            continue
        m = re.search(r"genes=([^;\n]+)", line)
        if m:
            ref_lists.append(m.group(1).split(","))
    known = {gene_xref[g.feature.id] for g in ann.genes
             if g.feature.id in gene_xref and any(t.is_coding for t in g.transcripts)}
    ours = {p for p in pairs_from(ours_lists) if p <= known}
    ref = {p for p in pairs_from(ref_lists) if p <= known}
    shared = ours & ref
    prec = len(shared) / len(ours) if ours else float("nan")
    rec = len(shared) / len(ref) if ref else float("nan")

    downstream = collections.Counter()
    for gid, leaders in gene_leader.items():
        kind = "SL2" if "SL2" in leaders else "SL1"
        pos = ours_pos.get(gid)
        downstream[(kind, "downstream" if pos and pos > 1 else ("leader" if pos == 1 else "not in operon"))] += 1

    lines = [
        ("sl_tss_with_reference_site_within_0", within(0), n_sl),
        ("sl_tss_with_reference_site_within_2", within(2), n_sl),
        ("sl_tss_with_reference_site_within_5", within(5), n_sl),
        ("nosl_tss_with_reference_site_within_5", nosl_near, n_nosl),
        ("operon_pairs_shared", len(shared), ""),
        ("operon_pairs_pipeline_only", len(ours - ref), ""),
        ("operon_pairs_reference_only", len(ref - ours), ""),
        ("operon_pair_precision", f"{prec:.3f}", ""),
        ("operon_pair_recall", f"{rec:.3f}", ""),
    ]
    with open(f"{args.prefix}.sl_operons.tsv", "w") as fh:
        fh.write("metric\tvalue\tof\n")
        for k, v, n in lines:
            fh.write(f"{k}\t{v}\t{n}\n")
        for k, v in sorted(agree.items()):
            fh.write(f"leader_vs_site_{k[0]}_{k[1]}\t{v}\t\n")
        for k, v in sorted(downstream.items()):
            fh.write(f"genes_{k[0]}_{k[1].replace(' ', '_')}\t{v}\t\n")
        for k in sorted((d for d in off if isinstance(d, int))):
            fh.write(f"sl_offset_{k}\t{off[k]}\t\n")
        for k in ("<-10", ">10", "none"):
            fh.write(f"sl_offset_{k}\t{off.get(k, 0)}\t\n")
    for k, v, n in lines:
        print(f"{k}\t{v}\t{n}")
    print("leader (ours) vs reference site type:", dict(agree))
    print("SL2/SL1 genes by operon position:", dict(downstream))


if __name__ == "__main__":
    main()
