process ANNOTATE_OPERONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(gff), path(tss)

    output:
    tuple val(meta), path("*.operons.gff3")    , emit: gff
    tuple val(meta), path("*.operons.bed")     , emit: bed
    tuple val(meta), path("*.operon_genes.tsv"), emit: tsv
    tuple val(meta), path("*_mqc.tsv")         , emit: mqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_operons.py \\
        --gff ${gff} \\
        --tss ${tss} \\
        ${args} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.operons.gff3 ${prefix}.operons.bed ${prefix}.operon_genes.tsv ${prefix}.operon_summary_mqc.tsv
    """
}
