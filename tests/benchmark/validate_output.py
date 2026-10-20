#!/usr/bin/env python3
"""Check structural invariants of a pipeline output GFF3 against its input annotation.

  validate_output.py --input annotation.gff3 --output <id>.utr.gff3

Checks, for every transcript of the output:
  - exons are sorted, do not overlap, have length >= 1, and the transcript spans them exactly;
  - every CDS row lies inside an exon, and the CDS rows are identical to the input's;
  - five_prime_UTR / three_prime_UTR rows equal the exonic parts outside the CDS;
  - trans_splice_acceptor_site rows sit on the transcript's 5' end and the TSS is tss_peak;
  - the gene spans all its transcripts;
  - IDs are unique (CDS rows may share an ID within one transcript) and every Parent exists;
  - transcripts and genes of the input are all present.
Prints a count per check and examples of violations; exits 1 if any check fails.
"""

import argparse
import collections
import gzip
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from gff3_utils import UTR3, UTR5, TSA, Annotation, Feature  # noqa: E402


def rows(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            yield Feature.from_line(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    fails = collections.Counter()
    examples = collections.defaultdict(list)

    def fail(check, what):
        fails[check] += 1
        if len(examples[check]) < 3:
            examples[check].append(what)

    # IDs and parents (raw rows)
    ids = collections.Counter()
    seen_multi = {}
    cds_owner = {}
    all_rows = list(rows(args.output))
    for f in all_rows:
        if f.id is None:
            continue
        if f.type == "CDS":
            owner = cds_owner.setdefault(f.id, f.parent)
            if owner != f.parent:
                fail("CDS ID shared across transcripts", f.id)
        else:
            # GFF3 multi-row (discontinuous) features share an ID: same type and Parent
            key = (f.type, f.parent, f.seqid)
            if f.id in seen_multi and seen_multi[f.id] == key:
                continue
            seen_multi[f.id] = key
            ids[f.id] += 1
    for i, n in ids.items():
        if n > 1:
            fail("duplicate ID", i)
    known = set(ids) | set(cds_owner)
    for f in all_rows:
        if f.parent and any(p not in known for p in f.parent.split(",")):
            fail("Parent not found", f"{f.type} {f.id} -> {f.parent}")

    inp = Annotation.read(args.input)
    out = Annotation.read(args.output)
    in_tx = {t.id: t for _, t in inp.transcripts()}
    out_tx = {t.id: t for _, t in out.transcripts()}
    for tid in in_tx:
        if tid not in out_tx:
            fail("input transcript missing", tid)
    in_genes = {g.feature.id for g in inp.genes}
    out_genes = {g.feature.id for g in out.genes}
    for gid in in_genes - out_genes:
        fail("input gene missing", gid)

    n_tx = 0
    for gene, t in out.transcripts():
        n_tx += 1
        ex = sorted(t.exons, key=lambda e: e.start)
        if any(e.end < e.start for e in ex):
            fail("exon length < 1", t.id)
        if any(a.end >= b.start for a, b in zip(ex, ex[1:])):
            fail("overlapping exons", t.id)
        if ex and (t.feature.start != ex[0].start or t.feature.end != ex[-1].end):
            fail("transcript does not span its exons", f"{t.id} {t.feature.start}-{t.feature.end} vs {ex[0].start}-{ex[-1].end}")
        for c in t.cds:
            if not any(e.start <= c.start and c.end <= e.end for e in ex):
                fail("CDS outside exons", f"{t.id} {c.start}-{c.end}")
        if t.id in in_tx:
            a = sorted((c.start, c.end, c.phase) for c in in_tx[t.id].cds)
            b = sorted((c.start, c.end, c.phase) for c in t.cds)
            if a != b:
                fail("CDS changed", t.id)
        if not (gene.feature.start <= t.feature.start and t.feature.end <= gene.feature.end):
            fail("gene does not span transcript", f"{gene.feature.id} / {t.id}")
        if gene.feature.strand != t.strand or gene.feature.seqid != t.seqid:
            fail("gene/transcript strand or sequence differ", t.id)
        for o in t.others:
            if o.type == TSA:
                five = t.feature.start if t.strand == "+" else t.feature.end
                if o.start != o.end or o.start != five:
                    fail("acceptor row not at the 5' end", t.id)
                peak = t.feature.get("tss_peak")
                if not peak or int(peak.split(":")[1]) != o.start:
                    fail("acceptor row != tss_peak", t.id)
        peak = t.feature.get("tss_peak")
        if peak and t.feature.get("tss_evidence") == "CAGE":
            five = t.feature.start if t.strand == "+" else t.feature.end
            if int(peak.split(":")[1]) != five:
                fail("tss_peak != 5' end", t.id)

    # UTR rows vs exon/CDS geometry
    utr_rows = collections.defaultdict(lambda: {UTR5: [], UTR3: []})
    for f in all_rows:
        if f.type in (UTR5, UTR3):
            utr_rows[f.parent][f.type].append((f.start, f.end))
    for _gene, t in out.coding_transcripts():
        e5, e3 = t.utr_segments()
        if sorted(utr_rows[t.id][UTR5]) != sorted(e5):
            fail("five_prime_UTR rows != exons outside CDS", t.id)
        if sorted(utr_rows[t.id][UTR3]) != sorted(e3):
            fail("three_prime_UTR rows != exons outside CDS", t.id)

    print(f"{n_tx} transcripts, {len(out.genes)} genes checked")
    for check, n in sorted(fails.items()):
        print(f"FAIL {check}: {n}  e.g. {examples[check]}")
    if not fails:
        print("all checks passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
