process FILTER_NON_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta) , path(redundant_ids)
    tuple val(meta2), path(seqs  , stageAs: "input_seqs/*")
    tuple val(meta3), path(models, stageAs: "input_hmm/*")
    tuple val(meta4), path(seeds , stageAs: "input_seed_msa/*")
    tuple val(meta5), path(full  , stageAs: "input_full_msa/*")

    output:
    tuple val(meta), path("fasta/*")   , emit: fasta
    tuple val(meta), path("hmm/*")     , emit: hmm
    tuple val(meta), path("seed_msa/*"), emit: seed_msa
    tuple val(meta), path("full_msa/*"), emit: full_msa
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_non_redundant_fams.py \\
        --redundant_ids ${redundant_ids} \\
        --seqs input_seqs \\
        --models input_hmm \\
        --seeds input_seed_msa \\
        --alns input_full_msa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
