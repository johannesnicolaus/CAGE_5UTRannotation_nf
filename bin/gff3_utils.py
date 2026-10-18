"""Minimal GFF3 model shared by the UTR annotation scripts.

Reads an EGAPx-style GFF3 (gene -> mRNA/lncRNA/transcript -> exon/CDS) into memory,
keeps every column and the attribute order of each row, lets the scripts edit
coordinates of mRNAs and their exons, regenerates the UTR features and writes the
file back in the original gene order.

Conventions (as in EGAPx GFF3 output):
  * coordinates are 1-based, inclusive;
  * the CDS includes the stop codon, so the 3' UTR starts right after the last CDS base;
  * a CDS may be split over several rows that share one ID.
"""

import bisect
import collections
import gzip
from urllib.parse import quote, unquote

# Feature types written by this pipeline; existing ones are dropped and regenerated.
UTR5 = "five_prime_UTR"
UTR3 = "three_prime_UTR"
TSA = "trans_splice_acceptor_site"
# UTR rows are always regenerated from exons and CDS, so any in the input are dropped.
DERIVED_TYPES = {UTR5, UTR3, "5'UTR", "3'UTR", "five_prime_utr", "three_prime_utr"}

TRANSCRIPT_TYPES = {"mRNA", "transcript", "lncRNA", "ncRNA", "primary_transcript", "snRNA", "snoRNA", "rRNA", "tRNA"}

# Reserved characters that must be escaped in GFF3 attribute values (column 9).
_GFF3_RESERVED = ";=&,\t\n\r%"


def escape(value):
    """Escape an attribute value for GFF3 column 9."""
    return "".join(quote(c) if c in _GFF3_RESERVED else c for c in str(value))


class Feature:
    """One GFF3 row. `attrs` keeps the raw (still escaped) values in file order."""

    __slots__ = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attrs")

    def __init__(self, seqid, source, ftype, start, end, score, strand, phase, attrs):
        self.seqid = seqid
        self.source = source
        self.type = ftype
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attrs = attrs

    @classmethod
    def from_line(cls, line):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"GFF3 line does not have 9 columns: {line!r}")
        attrs = collections.OrderedDict()
        for item in fields[8].strip().strip(";").split(";"):
            item = item.strip()
            if not item:
                continue
            key, _, value = item.partition("=")
            attrs[key] = value
        return cls(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], attrs)

    def get(self, key, default=None):
        """Unescaped attribute value."""
        value = self.attrs.get(key)
        return default if value is None else unquote(value)

    def set(self, key, value):
        """Set an attribute (value is escaped). Existing keys keep their position."""
        self.attrs[key] = escape(value)

    @property
    def id(self):
        return self.get("ID")

    @property
    def parent(self):
        return self.get("Parent")

    def to_line(self):
        attrs = ";".join(f"{k}={v}" for k, v in self.attrs.items())
        return "\t".join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.phase, attrs]
        )


class Transcript:
    """A transcript feature with its exon and CDS rows (and any other child rows)."""

    def __init__(self, feature):
        self.feature = feature
        self.exons = []
        self.cds = []
        self.others = []  # other children, kept as they are (except DERIVED_TYPES)
        self.modified = False  # set when exon coordinates were changed

    @property
    def id(self):
        return self.feature.id

    @property
    def strand(self):
        return self.feature.strand

    @property
    def seqid(self):
        return self.feature.seqid

    @property
    def is_coding(self):
        return bool(self.cds)

    def sort(self):
        self.exons.sort(key=lambda f: f.start)
        self.cds.sort(key=lambda f: f.start)

    def cds_bounds(self):
        """(min start, max end) of the CDS rows, stop codon included."""
        return min(f.start for f in self.cds), max(f.end for f in self.cds)

    def exon_intervals(self):
        return [(e.start, e.end) for e in self.exons]

    def five_prime_partial(self):
        """True if the model is open at its 5' end (EGAPx start_range/end_range)."""
        key = "start_range" if self.strand == "+" else "end_range"
        return key in self.feature.attrs

    def three_prime_partial(self):
        key = "end_range" if self.strand == "+" else "start_range"
        return key in self.feature.attrs

    def sync_bounds(self):
        """Set the transcript row to span its exons, and mark the transcript as modified."""
        self.modified = True
        if self.exons:
            self.feature.start = min(e.start for e in self.exons)
            self.feature.end = max(e.end for e in self.exons)

    def utr_segments(self):
        """Exonic (start, end) segments 5' and 3' of the CDS, in genomic order."""
        if not self.cds or not self.exons:
            return [], []
        cds_start, cds_end = self.cds_bounds()
        left, right = [], []
        for start, end in self.exon_intervals():
            if start < cds_start:
                left.append((start, min(end, cds_start - 1)))
            if end > cds_end:
                right.append((max(start, cds_end + 1), end))
        if self.strand == "-":
            return right, left
        return left, right

    def utr_lengths(self):
        utr5, utr3 = self.utr_segments()
        return sum(e - s + 1 for s, e in utr5), sum(e - s + 1 for s, e in utr3)


class Gene:
    def __init__(self, feature):
        self.feature = feature
        self.transcripts = []
        self.others = []  # non-transcript direct children

    def sync_bounds(self):
        spans = [(t.feature.start, t.feature.end) for t in self.transcripts]
        spans += [(f.start, f.end) for f in self.others]
        if spans:
            self.feature.start = min(s for s, _ in spans)
            self.feature.end = max(e for _, e in spans)


class Annotation:
    """Whole GFF3 file: header lines, genes in file order, and orphan rows."""

    def __init__(self):
        self.header = []
        self.genes = []
        self.orphans = []  # rows whose parent chain does not reach a gene; written unchanged

    @classmethod
    def read(cls, path):
        ann = cls()
        rows = []
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    # keep the leading header and every ##sequence-region pragma (EGAPx
                    # writes one before each sequence); drop '###' and later comments
                    pragma = line.rstrip("\n")
                    if pragma.strip() == "###":
                        continue
                    if not rows or pragma.startswith("##sequence-region"):
                        ann.header.append(pragma)
                    continue
                if not line.strip():
                    continue
                rows.append(Feature.from_line(line))

        genes, transcripts = {}, {}
        for f in rows:
            if f.parent is None:
                gene = Gene(f)
                ann.genes.append(gene)
                if f.id:
                    genes[f.id] = gene
        for f in rows:
            if f.parent is not None and f.parent in genes and f.type in TRANSCRIPT_TYPES:
                t = Transcript(f)
                genes[f.parent].transcripts.append(t)
                if f.id:
                    transcripts[f.id] = t
        for f in rows:
            parent = f.parent
            if parent is None or (parent in genes and f.type in TRANSCRIPT_TYPES):
                continue
            if parent in transcripts:
                t = transcripts[parent]
                if f.type == "exon":
                    t.exons.append(f)
                elif f.type == "CDS":
                    t.cds.append(f)
                elif f.type not in DERIVED_TYPES:
                    t.others.append(f)
            elif parent in genes:
                if f.type not in DERIVED_TYPES:
                    genes[parent].others.append(f)
            else:
                ann.orphans.append(f)
        for t in transcripts.values():
            t.sort()
        return ann

    def transcripts(self):
        for gene in self.genes:
            for t in gene.transcripts:
                yield gene, t

    def coding_transcripts(self):
        for gene, t in self.transcripts():
            if t.feature.type == "mRNA" and t.is_coding and t.exons:
                yield gene, t

    def cds_index(self):
        """Per (seqid, strand): sorted CDS intervals with their gene ID."""
        index = collections.defaultdict(list)
        for gene, t in self.transcripts():
            for c in t.cds:
                index[(c.seqid, c.strand)].append((c.start, c.end, gene.feature.id))
        return IntervalIndex(index)

    def write(self, path):
        """Write the annotation, regenerating five_prime_UTR / three_prime_UTR rows."""
        with open(path, "w") as out:
            for line in self.header:
                out.write(line + "\n")
            for gene in self.genes:
                # resize only genes with an edited transcript, so other EGAPx genes stay as they are
                if any(t.modified for t in gene.transcripts):
                    gene.sync_bounds()
                out.write(gene.feature.to_line() + "\n")
                for f in gene.others:
                    out.write(f.to_line() + "\n")
                for t in gene.transcripts:
                    out.write(t.feature.to_line() + "\n")
                    # 5' -> 3' order, as in EGAPx output
                    reverse = t.strand == "-"
                    for f in sorted(t.exons, key=lambda x: x.start, reverse=reverse):
                        out.write(f.to_line() + "\n")
                    for f in sorted(t.cds, key=lambda x: x.start, reverse=reverse):
                        out.write(f.to_line() + "\n")
                    for f in utr_features(t):
                        out.write(f.to_line() + "\n")
                    for f in t.others:
                        out.write(f.to_line() + "\n")
            for f in self.orphans:
                out.write(f.to_line() + "\n")


def feature_core_id(transcript):
    """ID stem for derived features: 'rna-egapxtmp_000001-R1' -> 'egapxtmp_000001-R1'."""
    tid = transcript.id or "transcript"
    return tid[4:] if tid.startswith("rna-") else tid


def utr_features(transcript):
    """five_prime_UTR / three_prime_UTR rows of a coding mRNA, numbered 5' -> 3'.

    The rows carry only ID and Parent; the evidence attributes are on the mRNA, so they are
    also there when a UTR has length 0 (e.g. an SL acceptor right at the ATG).
    """
    if transcript.feature.type != "mRNA" or not transcript.is_coding:
        return []
    utr5, utr3 = transcript.utr_segments()
    rows = []
    core = feature_core_id(transcript)
    for kind, segments in ((UTR5, utr5), (UTR3, utr3)):
        ordered = segments if transcript.strand != "-" else list(reversed(segments))
        for n, (start, end) in enumerate(ordered, 1):
            attrs = collections.OrderedDict()
            f = Feature(transcript.seqid, transcript.feature.source, kind, start, end, ".",
                        transcript.strand, ".", attrs)
            prefix = "utr5" if kind == UTR5 else "utr3"
            f.set("ID", f"{prefix}-{core}-{n}")
            f.set("Parent", transcript.id)
            rows.append(f)
    rows.sort(key=lambda f: f.start, reverse=transcript.strand == "-")
    return rows


class IntervalIndex:
    """Sorted intervals per key, for 'does [a, b] overlap anything not owned by X' queries."""

    def __init__(self, by_key):
        self.data = {}
        for key, items in by_key.items():
            items = sorted(items)
            starts = [s for s, _, _ in items]
            # running maximum of interval ends, to stop the left scan early
            max_end, running = [], 0
            for _, e, _ in items:
                running = max(running, e)
                max_end.append(running)
            self.data[key] = (items, starts, max_end)

    def overlaps(self, key, start, end, exclude_owner=None):
        """Intervals overlapping [start, end] whose owner is not exclude_owner."""
        if key not in self.data or end < start:
            return []
        items, starts, max_end = self.data[key]
        hits = []
        i = bisect.bisect_right(starts, end) - 1
        while i >= 0 and max_end[i] >= start:
            s, e, owner = items[i]
            if e >= start and s <= end and owner != exclude_owner:
                hits.append(items[i])
            i -= 1
        return hits


class PointIndex:
    """Sorted 1-based positions per (seqid, strand) with a payload, for window queries."""

    def __init__(self, items):
        by_key = collections.defaultdict(list)
        for seqid, strand, pos, payload in items:
            by_key[(seqid, strand)].append((pos, payload))
        self.data = {}
        for key, values in by_key.items():
            values.sort(key=lambda v: v[0])
            self.data[key] = ([p for p, _ in values], [v for _, v in values])

    def window(self, seqid, strand, start, end):
        """(pos, payload) with start <= pos <= end."""
        key = (seqid, strand)
        if key not in self.data or end < start:
            return []
        positions, payloads = self.data[key]
        lo = bisect.bisect_left(positions, start)
        hi = bisect.bisect_right(positions, end)
        return list(zip(positions[lo:hi], payloads[lo:hi]))


def exonic_length_between(intervals, a, b):
    """Number of exonic bases in [a, b] (genomic, inclusive) given sorted exon intervals."""
    if b < a:
        return 0
    total = 0
    for s, e in intervals:
        lo, hi = max(s, a), min(e, b)
        if lo <= hi:
            total += hi - lo + 1
    return total


def containing_interval(intervals, pos):
    """Index of the interval containing pos, or None."""
    for i, (s, e) in enumerate(intervals):
        if s <= pos <= e:
            return i
    return None


def read_fasta(path):
    """Yield (name, sequence) from a plain or gzipped FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)
