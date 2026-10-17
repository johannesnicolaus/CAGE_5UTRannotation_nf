process UTR_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(before, stageAs: 'before/*'), path(after)

    output:
    tuple val(meta), path("*.utr_per_transcript.tsv"), emit: tsv
    tuple val(meta), path("*_mqc.tsv")               , emit: mqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    utr_stats.py \\
        --before ${before} \\
        --after ${after} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.utr_per_transcript.tsv ${prefix}.utr_summary_mqc.tsv ${prefix}.utr5_length_mqc.tsv ${prefix}.utr3_length_mqc.tsv
    """
}
