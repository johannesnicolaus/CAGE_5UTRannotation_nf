//
// Uncompress the genome and annotation if needed, and collect the TSS cluster files: either
// a luscombeu/cagescan results directory, or BED files of TSS clusters from any CAGE workflow.
//

include { GUNZIP } from '../../../modules/nf-core/gunzip/main'

workflow PREPARE_INPUTS {

    take:
    ch_samplesheet // channel: [ meta, genome, annotation, tss, tss_sl ]
    cluster_set    // string: 'broad' or 'sharp'
    trans_splicing // boolean: use the spliced-leader counts when present

    main:

    //
    // Genome and annotation: gunzip the compressed ones (one GUNZIP call for both kinds)
    //
    def ch_files = ch_samplesheet
        .flatMap { meta, genome, annotation, _tss, _tss_sl ->
            [
                [ meta + [ file_type: 'genome' ], genome ],
                [ meta + [ file_type: 'annotation' ], annotation ],
            ]
        }
        .branch { _meta, f ->
            gz: f.name.endsWith('.gz')
            plain: true
        }

    GUNZIP(ch_files.gz)

    def ch_uncompressed = ch_files.plain
        .mix(GUNZIP.out.gunzip)
        .branch { meta, _f ->
            genome: meta.file_type == 'genome'
            annotation: meta.file_type == 'annotation'
        }

    def ch_fasta = ch_uncompressed.genome.map { meta, f -> [ meta - meta.subMap('file_type'), f ] }
    def ch_gff   = ch_uncompressed.annotation.map { meta, f -> [ meta - meta.subMap('file_type'), f ] }

    //
    // TSS clusters: a cagescan results directory (clusters, counts, SL counts when the dataset
    // has a spliced leader) or BED files (non-SL clusters, optional SL acceptor clusters)
    //
    def ch_cage = ch_samplesheet.map { meta, _genome, _annotation, tss, tss_sl ->
        tssInputFiles(meta, tss, tss_sl, cluster_set, trans_splicing)
    }

    emit:
    fasta = ch_fasta // channel: [ meta, fasta ]
    gff   = ch_gff   // channel: [ meta, gff ]
    cage  = ch_cage  // channel: [ meta, clusters.bed, counts.tsv or [], sl_counts.tsv or [], sl_leaders.tsv or [], sl_clusters.bed or [] ]
}

//
// TSS input files of one sample: [ meta, clusters, counts, sl_counts, sl_leaders, sl_bed ]
//
def tssInputFiles(meta, tss, tss_sl, cluster_set, trans_splicing) {
    def path = file(tss.toString())
    if (path.isDirectory()) {
        if (tss_sl) {
            error("Sample '${meta.id}': tss_sl is for BED input; a cagescan results directory carries its own SL counts.")
        }
        return cagescanTssFiles(meta, path, cluster_set, trans_splicing)
    }
    if (!tss_sl) {
        log.info("Sample '${meta.id}': no SL acceptor clusters (tss_sl); every TSS will be labelled nosl and operons are not annotated.")
    }
    def sl_bed = (tss_sl && trans_splicing) ? file(tss_sl.toString()) : []
    return [ meta, path, [], [], [], sl_bed ]
}

//
// Paths of the TSS cluster files of one cluster set in a cagescan results directory
// (cagescan docs/output.md: tss_clusters/<set>/all_samples.<set>.*).
//
def cagescanTssFiles(meta, cagescan_dir, cluster_set, trans_splicing) {
    def dir = file(cagescan_dir.toString()).resolve("tss_clusters/${cluster_set}")
    def clusters = dir.resolve("all_samples.${cluster_set}.tss_clusters.bed")
    def counts = dir.resolve("all_samples.${cluster_set}.tss_cluster_counts.tsv")
    def sl_counts = dir.resolve("all_samples.${cluster_set}.tss_cluster_sl_counts.tsv")
    def sl_leaders = dir.resolve("all_samples.${cluster_set}.tss_cluster_sl_leaders.tsv")
    [ clusters, counts ].each { required ->
        if (!required.exists()) {
            error("Sample '${meta.id}': ${required} not found. Is ${cagescan_dir} a luscombeu/cagescan results directory with the '${cluster_set}' cluster set (--cluster_set)?")
        }
    }
    if (!trans_splicing) {
        return [ meta, clusters, counts, [], [], [] ]
    }
    if (!sl_counts.exists()) {
        log.info("Sample '${meta.id}': no spliced-leader counts in ${dir}; every TSS will be labelled nosl and operons are not annotated.")
        return [ meta, clusters, counts, [], [], [] ]
    }
    return [ meta, clusters, counts, sl_counts, sl_leaders.exists() ? sl_leaders : [], [] ]
}
