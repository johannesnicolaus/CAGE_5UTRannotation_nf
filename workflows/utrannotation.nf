/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUTS         } from '../subworkflows/local/prepare_inputs'
include { CAGE_TSS_TABLE         } from '../modules/local/cage_tss_table/main'
include { SCAN_PAS               } from '../modules/local/scan_pas/main'
include { ANNOTATE_UTR5          } from '../modules/local/annotate_utr5/main'
include { ANNOTATE_UTR3          } from '../modules/local/annotate_utr3/main'
include { UTR_STATS              } from '../modules/local/utr_stats/main'
include { ANNOTATE_OPERONS       } from '../modules/local/annotate_operons/main'
include { AGAT_SPSTATISTICS      } from '../modules/nf-core/agat/spstatistics/main'
include { AGAT_CONVERTSPGFF2GTF  } from '../modules/nf-core/agat/convertspgff2gtf/main'
include { GFFREAD                } from '../modules/nf-core/gffread/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_utrannotation_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTRANNOTATION {

    take:
    ch_samplesheet // channel: [ meta, genome, annotation, cagescan_dir ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: uncompress the genome and annotation, find the cagescan TSS cluster files
    //
    PREPARE_INPUTS(ch_samplesheet, params.cluster_set)
    def ch_fasta = PREPARE_INPUTS.out.fasta  // [ meta, fasta ]
    def ch_gff   = PREPARE_INPUTS.out.gff    // [ meta, gff ]

    //
    // MODULE: one TSS table per annotation, each TSS labelled sl / nosl
    //
    CAGE_TSS_TABLE(PREPARE_INPUTS.out.cage)

    //
    // MODULE: polyadenylation-signal motifs in the genome
    //
    SCAN_PAS(ch_fasta)

    //
    // MODULE: 5' ends from CAGE TSSs
    //
    ANNOTATE_UTR5(
        ch_gff
            .join(CAGE_TSS_TABLE.out.tsv, failOnMismatch: true)
            .join(ch_fasta, failOnMismatch: true)
    )

    //
    // MODULE: 3' ends from poly(A) evidence or polyadenylation signals
    //
    ANNOTATE_UTR3(
        ANNOTATE_UTR5.out.gff.join(SCAN_PAS.out.bed, failOnMismatch: true)
    )
    def ch_utr_gff = ANNOTATE_UTR3.out.gff

    //
    // MODULE: operons, only for genomes whose CAGE data have spliced-leader reads
    // (without SL every TSS is nosl and would split every operon)
    //
    def ch_has_sl = PREPARE_INPUTS.out.cage.map { meta, _clusters, _counts, sl_counts -> [ meta, sl_counts ? true : false ] }
    ANNOTATE_OPERONS(
        ch_utr_gff
            .join(CAGE_TSS_TABLE.out.tsv, failOnMismatch: true)
            .join(ch_has_sl, failOnMismatch: true)
            .filter { _meta, _gff, _tss, has_sl -> has_sl && !params.skip_operons }
            .map { meta, gff, tss, _has_sl -> [ meta, gff, tss ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_OPERONS.out.mqc.map { _meta, file -> file })

    //
    // MODULE: UTR lengths and evidence, before and after
    //
    UTR_STATS(ch_gff.join(ch_utr_gff, failOnMismatch: true))
    ch_multiqc_files = ch_multiqc_files.mix(UTR_STATS.out.mqc.flatMap { _meta, files -> files })

    //
    // MODULE: AGAT statistics of the input and output annotations
    //
    AGAT_SPSTATISTICS(
        ch_gff.map { meta, gff -> [ meta + [ stage: 'egapx' ], gff ] }
            .mix(ch_utr_gff.map { meta, gff -> [ meta + [ stage: 'utr' ], gff ] })
    )

    //
    // MODULE: GTF version of the output
    //
    AGAT_CONVERTSPGFF2GTF(ch_utr_gff)

    //
    // MODULE: transcript and protein sequences of the output annotation
    //
    def ch_gffread_in = ch_utr_gff
        .join(ch_fasta, failOnMismatch: true)
        .flatMap { meta, gff, fasta ->
            [
                [ meta + [ seqtype: 'transcripts' ], gff, fasta ],
                [ meta + [ seqtype: 'proteins' ], gff, fasta ],
            ]
        }
        .multiMap { meta, gff, fasta ->
            gff: [ meta, gff ]
            fasta: fasta
        }
    GFFREAD(ch_gffread_in.gff, ch_gffread_in.fasta)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'utrannotation_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'utrannotation'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    utr_gff        = ch_utr_gff                                               // channel: [ meta, gff3 ]
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                              // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
