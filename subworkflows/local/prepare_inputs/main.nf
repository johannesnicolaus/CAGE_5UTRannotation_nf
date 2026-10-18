//
// Uncompress the genome and annotation if needed, and find the TSS cluster files in a
// luscombeu/cagescan results directory.
//

include { GUNZIP } from '../../../modules/nf-core/gunzip/main'

workflow PREPARE_INPUTS {

    take:
    ch_samplesheet // channel: [ meta, genome, annotation, cagescan_dir ]
    cluster_set    // string: 'broad' or 'sharp'
    trans_splicing // boolean: use the spliced-leader counts when present

    main:

    //
    // Genome and annotation: gunzip the compressed ones (one GUNZIP call for both kinds)
    //
    def ch_files = ch_samplesheet
        .flatMap { meta, genome, annotation, _cagescan_dir ->
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
    // cagescan outputs: TSS clusters, counts and (only for datasets with a spliced leader) SL counts
    //
    def ch_cage = ch_samplesheet.map { meta, _genome, _annotation, cagescan_dir ->
        cagescanTssFiles(meta, cagescan_dir, cluster_set, trans_splicing)
    }

    emit:
    fasta = ch_fasta // channel: [ meta, fasta ]
    gff   = ch_gff   // channel: [ meta, gff ]
    cage  = ch_cage  // channel: [ meta, tss_clusters.bed, tss_cluster_counts.tsv, tss_cluster_sl_counts.tsv or [], tss_cluster_sl_leaders.tsv or [] ]
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
        return [ meta, clusters, counts, [], [] ]
    }
    if (!sl_counts.exists()) {
        log.info("Sample '${meta.id}': no spliced-leader counts in ${dir}; every TSS will be labelled nosl and operons are not annotated.")
        return [ meta, clusters, counts, [], [] ]
    }
    return [ meta, clusters, counts, sl_counts, sl_leaders.exists() ? sl_leaders : [] ]
}
