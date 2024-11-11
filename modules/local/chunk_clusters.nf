process CHUNK_CLUSTERS {
    tag "$meta.id"

    conda "conda-forge::biopython=1.84 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e25bf39c3ebfbf72a1ddcaeea8f246c03e6da84b240ec0f3b2814d9719ec66b1/data' :
        'community.wave.seqera.io/library/biopython_pandas:641c5796a40a11b9' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)
    val(size_threshold)

    output:
    tuple val(meta), path("chunked_fasta/*"), emit: fasta_chunks
    path "versions.yml"                     , emit: versions

    script:
    """
    chunk_clusters.py --clustering ${clustering} \
        --sequences ${sequences} \
        --threshold ${size_threshold} \
        --out_folder chunked_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
