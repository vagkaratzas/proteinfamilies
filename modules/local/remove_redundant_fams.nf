process REMOVE_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta) , path(mapping)
    tuple val(meta2), path(domtbl)
    tuple val(meta3), path(fasta, stageAs: "fasta_folder/*")
    val(length_threshold)

    output:
    tuple val(meta), path("non_redundant/*"), emit: fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    remove_redundant_fams.py \\
        --mapping ${mapping} \\
        --domtbl ${domtbl} \\
        --fasta_folder fasta_folder \\
        --length_threshold ${length_threshold} \\
        --out_folder non_redundant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
