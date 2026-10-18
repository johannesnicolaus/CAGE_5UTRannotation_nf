#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    luscombeu/utrannotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/luscombeu/utrannotation
----------------------------------------------------------------------------------------
*/

//
// Typed declarations of the boolean and numeric parameters (defaults: nextflow.config). From Nextflow 26.04
// command-line values arrive as strings unless declared with a type here ('--flag false' would be truthy).
//
params {
    tss_max_dist        : Integer = 500
    tss_min_rel         : Float   = 0.1
    tss_rel_whole_utr   : Boolean = true
    tss_noncoding       : Boolean = true
    tss_min_count       : Float   = 0.0
    utr5_min_len        : Integer = 0
    sl_min_fraction     : Float   = 0.5
    sl_window_up        : Integer = 5
    sl_window_down      : Integer = 5
    keep_polya_evidence : Boolean = true
    pas_window          : Integer = 500
    pas_offset          : Integer = 15
    trans_splicing      : Boolean = true
    operons             : Boolean = true
    operon_max_gap      : Integer = 500
    operon_max_gap_tx   : Integer = 100
    operon_split_min_count : Float = 0.0
    version             : Boolean = false
    plaintext_email     : Boolean = false
    monochrome_logs     : Boolean = false
    validate_params     : Boolean = true
    help_full           : Boolean = false
    show_hidden         : Boolean = false
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTRANNOTATION  } from './workflows/utrannotation'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_utrannotation_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_utrannotation_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow LUSCOMBEU_UTRANNOTATION {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    UTRANNOTATION (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = UTRANNOTATION.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    LUSCOMBEU_UTRANNOTATION (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        LUSCOMBEU_UTRANNOTATION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
