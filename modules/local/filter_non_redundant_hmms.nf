process FILTER_NON_REDUNDANT_HMMS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.13.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta) , path(seqs, stageAs: "seqs/*")
    tuple val(meta2), path(models, stageAs: "models/*")

    output:
    tuple val(meta), path("non_redundant/*"), emit: hmm
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_non_redundant_hmms.py \\
        --seqs seqs \\
        --models models \\
        --out_folder non_redundant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
