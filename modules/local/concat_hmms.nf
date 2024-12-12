process CONCAT_HMMS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(hmm, stageAs: "hmm/*")

    output:
    tuple val(meta), path("${meta.id}.hmm.gz"), emit: hmm
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    zcat hmm/*.hmm.gz | gzip > ${meta.id}.hmm.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
