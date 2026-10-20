# Accuracy benchmark against a trusted annotation

These scripts measure how well the pipeline restores UTRs. They need an annotation whose UTR ends are trusted. The first one used is the _C. elegans_ WormBase WS298 gene set in NCBI RefSeq format (GCF_000002985.6), whose 3' ends come largely from poly(A)-site sequencing (Mangone et al. 2010; Jan et al. 2011).

1. Make an EGAPx-like input. CDSs are kept and models open at an end (`start_range`/`end_range`) are left alone:
   - `--mode strip`: no UTRs;
   - `--mode extend`: UTRs pushed outward by 0–500 nt.

   ```bash
   tests/benchmark/degrade_annotation.py --gff reference.gff3 --mode strip  --output strip.gff3
   tests/benchmark/degrade_annotation.py --gff reference.gff3 --mode extend --output extend.gff3
   ```

2. Run the pipeline on each degraded annotation, with CAGE for the same genome.
3. Score the output against the reference. Ends that are open in the reference are skipped.

   ```bash
   tests/benchmark/compare_utrs.py --reference reference.gff3 --input strip.gff3 \
       --output results/annotation/<id>.utr.gff3 --pas results/pas/<id>.pas.bed --prefix strip
   ```

   The summary gives, for the input (baseline) and the output, the share of 5' and 3' ends within 0 / 5 / 20 / 50 nt of the reference, per evidence class. `<prefix>.pas_offset.tsv` gives where the reference ends lie relative to the polyadenylation-signal motif, which shows the best `--pas_offset`.
