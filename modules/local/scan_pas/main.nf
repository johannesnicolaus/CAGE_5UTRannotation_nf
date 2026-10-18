process SCAN_PAS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.pas.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    scan_pas.py \\
        --fasta ${fasta} \\
        ${args} \\
        --output ${prefix}.pas.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pas.bed
    """
}
