process ANNOTATE_UTR5 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(gff), path(tss), path(fasta)

    output:
    tuple val(meta), path("*.utr5.gff3"), emit: gff
    tuple val(meta), path("*.utr5.log") , emit: log
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_utr5.py \\
        --gff ${gff} \\
        --tss ${tss} \\
        --fasta ${fasta} \\
        ${args} \\
        --output ${prefix}.utr5.gff3 \\
        2> ${prefix}.utr5.log \\
        || { cat ${prefix}.utr5.log >&2; exit 1; }
    cat ${prefix}.utr5.log >&2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.utr5.gff3 ${prefix}.utr5.log
    """
}
