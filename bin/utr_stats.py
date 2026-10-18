#!/usr/bin/env python3
"""Compare UTRs before (EGAPx) and after the pipeline, per transcript and in summary.

  utr_stats.py --before egapx.gff3 --after PREFIX.utr.gff3 --prefix PREFIX

Outputs:
  PREFIX.utr_per_transcript.tsv  one row per coding mRNA: UTR lengths before and after,
                                 evidence and TSS / PAS attributes
  PREFIX.utr_summary_mqc.tsv     MultiQC table, one row per annotation
  PREFIX.utr5_length_mqc.tsv     MultiQC bar graph of 5' UTR lengths, before and after
  PREFIX.utr3_length_mqc.tsv     the same for 3' UTRs
"""

import argparse
import collections
import statistics

from gff3_utils import Annotation

BINS = [(0, 0, "0"), (1, 50, "1-50"), (51, 100, "51-100"), (101, 200, "101-200"), (201, 500, "201-500"),
        (501, 1000, "501-1000"), (1001, None, ">1000")]

COLUMNS = ["transcript_id", "gene_id", "seqid", "strand",
           "utr5_len_egapx", "utr5_len", "tss_evidence", "tss_type", "tss_id", "tss_peak", "tss_count",
           "sl_fraction", "acceptor_AG",
           "utr3_len_egapx", "utr3_len", "utr3_evidence", "pas_motif", "pas_pos"]


def bin_label(length):
    for lo, hi, label in BINS:
        if length >= lo and (hi is None or length <= hi):
            return label
    return BINS[-1][2]


def summary(values):
    if not values:
        return "NA", "NA"
    return f"{statistics.median(values):g}", f"{statistics.mean(values):.1f}"


def pct(n, d):
    return f"{100 * n / d:.1f}" if d else "NA"


def write_bargraph(path, section_id, title, rows):
    with open(path, "w") as out:
        out.write(f"# id: '{section_id}'\n")
        out.write(f'# section_name: "{title}"\n')
        out.write("# description: 'Number of coding mRNAs per UTR length (bases), in the EGAPx input and in the output.'\n")
        out.write("# plot_type: 'bargraph'\n")
        out.write("# pconfig:\n")
        out.write(f"#     id: '{section_id}_plot'\n")
        out.write(f'#     title: "{title}"\n')
        out.write("#     ylab: 'mRNAs'\n")
        out.write("#     cpswitch_counts_label: 'mRNAs'\n")
        out.write("Sample\t" + "\t".join(label for _, _, label in BINS) + "\n")
        for name, counts in rows:
            out.write(name + "\t" + "\t".join(str(counts.get(label, 0)) for _, _, label in BINS) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--before", required=True)
    parser.add_argument("--after", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    before = {t.id: t.utr_lengths() for _, t in Annotation.read(args.before).coding_transcripts()}
    after = Annotation.read(args.after)
    noncoding = [t for _, t in after.noncoding_transcripts()]
    nc_cage = [t for t in noncoding if t.feature.get("tss_evidence") == "CAGE"]
    nc_sl = [t for t in nc_cage if t.feature.get("tss_type") == "sl"]
    rows = []
    for gene, t in after.coding_transcripts():
        f = t.feature
        utr5, utr3 = t.utr_lengths()
        utr5_0, utr3_0 = before.get(t.id, (None, None))
        rows.append({
            "transcript_id": t.id, "gene_id": gene.feature.id, "seqid": t.seqid, "strand": t.strand,
            "utr5_len_egapx": utr5_0, "utr5_len": utr5, "utr3_len_egapx": utr3_0, "utr3_len": utr3,
            **{k: f.get(k) for k in COLUMNS if k in f.attrs},
        })

    with open(f"{args.prefix}.utr_per_transcript.tsv", "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join("NA" if row.get(c) is None else str(row[c]) for c in COLUMNS) + "\n")

    n = len(rows)
    cage = [r for r in rows if r.get("tss_evidence") == "CAGE"]
    sl = [r for r in cage if r.get("tss_type") == "sl"]
    sl_ag = [r for r in sl if r.get("acceptor_AG") == "true"]
    shifted5 = [r for r in cage if r["utr5_len_egapx"] is not None]
    ext5 = sum(r["utr5_len"] > r["utr5_len_egapx"] for r in shifted5)
    trim5 = sum(r["utr5_len"] < r["utr5_len_egapx"] for r in shifted5)
    ev3 = collections.Counter(r.get("utr3_evidence") for r in rows)
    pas = [r for r in rows if r.get("utr3_evidence") == "PAS" and r["utr3_len_egapx"] is not None]
    ext3 = sum(r["utr3_len"] > r["utr3_len_egapx"] for r in pas)
    trim3 = sum(r["utr3_len"] < r["utr3_len_egapx"] for r in pas)

    med5_0, mean5_0 = summary([r["utr5_len_egapx"] for r in rows if r["utr5_len_egapx"]])
    med5, mean5 = summary([r["utr5_len"] for r in rows if r["utr5_len"]])
    med3_0, mean3_0 = summary([r["utr3_len_egapx"] for r in rows if r["utr3_len_egapx"]])
    med3, mean3 = summary([r["utr3_len"] for r in rows if r["utr3_len"]])

    header = ["Sample", "Coding mRNAs", "CAGE TSS (%)", "SL TSS (% of CAGE)", "SL with AG acceptor (%)",
              "5' UTR extended", "5' UTR trimmed",
              "5' UTR median EGAPx", "5' UTR median", "5' UTR mean EGAPx", "5' UTR mean",
              "3' polyA reads (%)", "3' PAS (%)", "3' EGAPx (%)", "3' UTR extended", "3' UTR trimmed",
              "3' UTR median EGAPx", "3' UTR median", "3' UTR mean EGAPx", "3' UTR mean",
              "Non-coding transcripts", "Non-coding CAGE TSS (%)", "Non-coding SL TSS (% of CAGE)"]
    values = [args.prefix, n, pct(len(cage), n), pct(len(sl), len(cage)), pct(len(sl_ag), len(sl)),
              ext5, trim5, med5_0, med5, mean5_0, mean5,
              pct(ev3["polyA_reads"], n), pct(ev3["PAS"], n), pct(ev3["EGAPx"], n), ext3, trim3,
              med3_0, med3, mean3_0, mean3,
              len(noncoding), pct(len(nc_cage), len(noncoding)), pct(len(nc_sl), len(nc_cage))]
    with open(f"{args.prefix}.utr_summary_mqc.tsv", "w") as out:
        out.write("# id: 'utr_summary'\n")
        out.write("# section_name: 'UTR annotation'\n")
        out.write('# description: "Coding mRNAs whose 5\' end moved to a CAGE TSS (sl: spliced-leader acceptor) '
                  'and the evidence behind each 3\' end. UTR lengths are medians and means over mRNAs with a UTR, '
                  'in the EGAPx input and in the output."\n')
        out.write("# plot_type: 'table'\n")
        out.write("# pconfig:\n")
        out.write("#     id: 'utr_summary_table'\n")
        out.write("#     title: 'UTR annotation'\n")
        out.write("\t".join(header) + "\n")
        out.write("\t".join(str(v) for v in values) + "\n")

    for kind, before_key, after_key in (("utr5", "utr5_len_egapx", "utr5_len"), ("utr3", "utr3_len_egapx", "utr3_len")):
        counts_before = collections.Counter(bin_label(r[before_key]) for r in rows if r[before_key] is not None)
        counts_after = collections.Counter(bin_label(r[after_key]) for r in rows)
        label = "5' UTR length" if kind == "utr5" else "3' UTR length"
        write_bargraph(f"{args.prefix}.{kind}_length_mqc.tsv", f"{kind}_length", label,
                       [(f"{args.prefix} EGAPx", counts_before), (f"{args.prefix} output", counts_after)])


if __name__ == "__main__":
    main()
