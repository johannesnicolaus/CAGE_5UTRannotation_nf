#!/usr/bin/env python3
"""Rule tests for the bin/ scripts on small synthetic annotations.

  python3 tests/scripts/test_rules.py [-v]

Each case writes a tiny genome, GFF3, TSS table and PAS BED, runs the scripts as the
pipeline does, and checks the chosen 5' / 3' ends and operons. Coordinates are 1-based.
"""

import os
import re
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin")
sys.path.insert(0, BIN)
from gff3_utils import Annotation  # noqa: E402

SEQ = "chr1"
GENOME_LEN = 6000


def gene_rows(gid, strand, tx, cds, exons=None, biotype="protein_coding"):
    """GFF3 rows of a one-transcript coding gene; tx and cds are (start, end)."""
    exons = exons or [tx]
    a = f"{SEQ}\tT\t"
    rows = [f"{a}gene\t{tx[0]}\t{tx[1]}\t.\t{strand}\t.\tID=gene-{gid};gene_biotype={biotype}",
            f"{a}mRNA\t{tx[0]}\t{tx[1]}\t.\t{strand}\t.\tID=rna-{gid};Parent=gene-{gid}"]
    for i, (s, e) in enumerate(exons, 1):
        rows.append(f"{a}exon\t{s}\t{e}\t.\t{strand}\t.\tID=exon-{gid}-{i};Parent=rna-{gid}")
    rows.append(f"{a}CDS\t{cds[0]}\t{cds[1]}\t.\t{strand}\t0\tID=cds-{gid};Parent=rna-{gid}")
    return rows


def tss_row(pos, strand, kind, count, tid=None):
    sl = count if kind == "sl" else 0
    return f"{SEQ}\t{pos}\t{strand}\t{tid or f'{kind}{pos}'}\t{count}\t{sl}\t{1.0 if sl else 0.0:.3f}\t{kind}\tNA"


class Case:
    def __init__(self, genes, tss, pas=(), ag_before=()):
        self.dir = tempfile.mkdtemp(prefix="utr_rules_")
        self.p = lambda name: os.path.join(self.dir, name)
        seq = ["C"] * GENOME_LEN
        for pos in ag_before:  # AG just before a + strand TSS at pos
            seq[pos - 3:pos - 1] = ["A", "G"]
        with open(self.p("genome.fa"), "w") as fh:
            fh.write(f">{SEQ}\n" + "".join(seq) + "\n")
        with open(self.p("in.gff3"), "w") as fh:
            fh.write("##gff-version 3\n" + "\n".join(r for g in genes for r in g) + "\n")
        with open(self.p("tss.tsv"), "w") as fh:
            fh.write("chrom\tpeak\tstrand\ttss_id\tcount\tsl_count\tsl_fraction\ttss_type\tsl_leader\n")
            fh.write("".join(r + "\n" for r in tss))
        with open(self.p("pas.bed"), "w") as fh:  # (motif start, end, strand), 1-based
            fh.write("".join(f"{SEQ}\t{s - 1}\t{e}\tAATAAA\t0\t{st}\n" for s, e, st in pas))

    def run(self, script, *args):
        subprocess.run([sys.executable, os.path.join(BIN, script), *args], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def utr5(self, *extra):
        self.run("annotate_utr5.py", "--gff", self.p("in.gff3"), "--tss", self.p("tss.tsv"),
                 "--fasta", self.p("genome.fa"), "--output", self.p("utr5.gff3"), "--min-utr-len", "0", *extra)
        return {t.id: t for _, t in Annotation.read(self.p("utr5.gff3")).transcripts()}

    def utr3(self):
        self.utr5()
        self.run("annotate_utr3.py", "--gff", self.p("utr5.gff3"), "--pas", self.p("pas.bed"),
                 "--output", self.p("utr3.gff3"))
        return {t.id: t for _, t in Annotation.read(self.p("utr3.gff3")).transcripts()}

    def operons(self, *extra):
        self.run("annotate_operons.py", "--gff", self.p("in.gff3"), "--tss", self.p("tss.tsv"),
                 "--prefix", self.p("op"), *extra)
        with open(self.p("op.operons.gff3")) as fh:
            self.operon_rows = [line for line in fh if not line.startswith("#") and line.strip()]
        with open(self.p("op.operons.gff3")) as fh:
            return [re.search(r"(?:\t|;)genes=([^;\n]+)", line).group(1).split(",")
                    for line in fh if not line.startswith("#") and line.strip()]


def five(t):
    return t.feature.start if t.strand == "+" else t.feature.end


class FivePrime(unittest.TestCase):
    def test_strong_shadow_does_not_filter_out_acceptor_plus(self):
        # nosl shadow at -2 of the acceptor, 20x stronger: the acceptor must still win
        c = Case([gene_rows("A", "+", (801, 1500), (1001, 1300))],
                 [tss_row(899, "+", "nosl", 200), tss_row(901, "+", "sl", 10)])
        t = c.utr5()["rna-A"]
        self.assertEqual((five(t), t.feature.get("tss_type")), (901, "sl"))

    def test_strong_shadow_does_not_filter_out_acceptor_minus(self):
        c = Case([gene_rows("A", "-", (1501, 2200), (1701, 2000))],
                 [tss_row(2102, "-", "nosl", 200), tss_row(2100, "-", "sl", 10)])
        t = c.utr5()["rna-A"]
        self.assertEqual((five(t), t.feature.get("tss_type")), (2100, "sl"))

    def test_auto_excludes_nosl_upstream_of_any_acceptor(self):
        # nosl between two acceptors is upstream of the second one: it may be its outron start
        c = Case([gene_rows("A", "+", (801, 1500), (1001, 1300))],
                 [tss_row(900, "+", "sl", 50), tss_row(920, "+", "nosl", 100), tss_row(950, "+", "sl", 50)])
        t = c.utr5()["rna-A"]
        self.assertEqual((five(t), t.feature.get("tss_type")), (950, "sl"))
        self.assertEqual(t.feature.get("nosl_tss_peak"), f"{SEQ}:920:+")

    def test_stronger_nosl_downstream_of_acceptors_wins(self):
        c = Case([gene_rows("A", "+", (801, 1500), (1001, 1300))],
                 [tss_row(900, "+", "sl", 50), tss_row(960, "+", "nosl", 100)])
        t = c.utr5()["rna-A"]
        self.assertEqual((five(t), t.feature.get("tss_type")), (960, "nosl"))
        self.assertEqual(t.feature.get("sl_acceptor_peak"), f"{SEQ}:900:+")

    def test_weak_exonic_peak_ignored(self):
        c = Case([gene_rows("A", "+", (801, 1500), (1001, 1300))],
                 [tss_row(850, "+", "nosl", 100), tss_row(990, "+", "nosl", 5)])
        self.assertEqual(five(c.utr5()["rna-A"]), 850)

    def test_rerun_clears_old_attributes(self):
        c = Case([gene_rows("A", "+", (801, 1500), (1001, 1300))], [tss_row(901, "+", "sl", 10)])
        c.utr5()
        os.replace(c.p("utr5.gff3"), c.p("in.gff3"))
        with open(c.p("tss.tsv"), "w") as fh:
            fh.write("chrom\tpeak\tstrand\ttss_id\tcount\tsl_count\tsl_fraction\ttss_type\tsl_leader\n")
        t = c.utr5()["rna-A"]
        self.assertEqual(t.feature.get("tss_evidence"), "EGAPx")
        self.assertIsNone(t.feature.get("tss_type"))
        self.assertFalse(any(o.type == "trans_splice_acceptor_site" for o in t.others))


class ThreePrime(unittest.TestCase):
    def test_no_extension_into_next_gene_moved_by_5prime_step(self):
        # B's CAGE TSS (1350) lies inside A's 3' UTR (ends 1400); A's PAS would extend A to 1520
        c = Case([gene_rows("A", "+", (801, 1400), (1001, 1300)), gene_rows("B", "+", (1650, 2100), (1700, 2000))],
                 [tss_row(850, "+", "nosl", 10), tss_row(1350, "+", "nosl", 10)],
                 pas=[(1500, 1505, "+")])
        out = c.utr3()
        self.assertEqual(five(out["rna-B"]), 1350)
        self.assertEqual(out["rna-A"].feature.end, 1400)
        self.assertEqual(out["rna-A"].feature.get("utr3_evidence"), "EGAPx")

    def test_trimming_pas_still_allowed(self):
        c = Case([gene_rows("A", "+", (801, 1400), (1001, 1300)), gene_rows("B", "+", (1650, 2100), (1700, 2000))],
                 [tss_row(850, "+", "nosl", 10), tss_row(1350, "+", "nosl", 10)],
                 pas=[(1320, 1325, "+")])
        a = c.utr3()["rna-A"]
        self.assertEqual((a.feature.end, a.feature.get("utr3_evidence")), (1340, "PAS"))

    def test_pas_extends_without_neighbour(self):
        c = Case([gene_rows("A", "-", (1601, 2200), (1701, 2000))], [], pas=[(1500, 1505, "-")])
        a = c.utr3()["rna-A"]
        self.assertEqual((a.feature.start, a.feature.get("utr3_evidence")), (1485, "PAS"))


class Operons(unittest.TestCase):
    genes = [gene_rows("A", "+", (1001, 1600), (1001, 1600)), gene_rows("B", "+", (1801, 2400), (1801, 2400)),
             gene_rows("C", "+", (2601, 3200), (2601, 3200))]

    def test_operon_with_sl(self):
        c = Case(self.genes, [tss_row(1790, "+", "sl", 10), tss_row(2590, "+", "sl", 10)])
        self.assertEqual(c.operons(), [["gene-A", "gene-B", "gene-C"]])

    def test_split_at_dominant_nosl(self):
        c = Case(self.genes, [tss_row(1790, "+", "sl", 10), tss_row(2550, "+", "nosl", 100), tss_row(2590, "+", "sl", 10)])
        self.assertEqual(c.operons(), [["gene-A", "gene-B"]])

    def test_minor_promoter_of_sl_gene_does_not_split(self):
        # hybrid-operon pattern: C is mainly trans-spliced, with a weaker own promoter
        tss = [tss_row(1790, "+", "sl", 10), tss_row(2550, "+", "nosl", 10), tss_row(2590, "+", "sl", 50)]
        self.assertEqual(Case(self.genes, tss).operons(), [["gene-A", "gene-B", "gene-C"]])
        self.assertEqual(Case(self.genes, tss).operons("--split-rule", "any"), [["gene-A", "gene-B"]])

    def test_operon_labels(self):
        c = Case(self.genes, [tss_row(1790, "+", "sl", 10)])
        c.operons()
        row = c.operon_rows[0]
        self.assertIn("downstream_sl=1/2", row)
        self.assertIn("operon_evidence=partial_sl", row)
        self.assertIn("sl_upstream=true,false", row)
        self.assertIn("sl_source=SL_reads", row)

    def test_shadow_of_acceptor_does_not_split(self):
        c = Case(self.genes, [tss_row(1790, "+", "sl", 10), tss_row(2588, "+", "nosl", 100), tss_row(2590, "+", "sl", 10)])
        self.assertEqual(c.operons(), [["gene-A", "gene-B", "gene-C"]])

    def test_weak_acceptor_does_not_shield_with_min_count(self):
        c = Case(self.genes, [tss_row(1790, "+", "sl", 10), tss_row(2588, "+", "nosl", 100), tss_row(2590, "+", "sl", 1)])
        self.assertEqual(c.operons("--split-min-count", "5"), [["gene-A", "gene-B"]])

    def test_no_sl_detected_distance_only_without_genome(self):
        c = Case(self.genes, [tss_row(1790, "+", "nosl", 10)])
        self.assertEqual(c.operons(), [["gene-A", "gene-B", "gene-C"]])
        self.assertIn("sl_source=none", c.operon_rows[0])

    def test_no_sl_detected_ag_inferred_acceptors(self):
        # 1790 follows an AG: an undetected acceptor, no split; 2550 does not: a promoter, split
        c = Case(self.genes, [tss_row(1790, "+", "nosl", 10), tss_row(2550, "+", "nosl", 10)], ag_before=[1790])
        self.assertEqual(c.operons("--fasta", c.p("genome.fa")), [["gene-A", "gene-B"]])
        self.assertIn("sl_source=AG_inferred", c.operon_rows[0])
        self.assertIn("downstream_sl=1/1", c.operon_rows[0])


if __name__ == "__main__":
    unittest.main()
