process CAGE_TSS_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(clusters), path(counts), path(sl_counts)

    output:
    tuple val(meta), path("*.tss.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sl_arg = sl_counts ? "--sl-counts ${sl_counts}" : ''
    """
    cage_tss_table.py \\
        --clusters ${clusters} \\
        --counts ${counts} \\
        ${sl_arg} \\
        ${args} \\
        --output ${prefix}.tss.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'chrom\\tpeak\\tstrand\\ttss_id\\tcount\\tsl_count\\tsl_fraction\\ttss_type\\n' > ${prefix}.tss.tsv
    """
}
