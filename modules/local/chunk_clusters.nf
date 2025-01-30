process CHUNK_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::biopython=1.84 conda-forge::polars=1.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/9266f998993b43dd4ac4ee5f75aa3217cd05b0d2f38df1dc110ec9c5832627fd/data' :
        'community.wave.seqera.io/library/biopython_polars:42dc770aaa311a0c' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)
    val(size_threshold)

    output:
    tuple val(meta), path("chunked_fasta/*"), emit: fasta_chunks
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = sequences.getName().endsWith(".gz") ? true : false
    def fasta_name    = sequences.name.replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $sequences > $fasta_name
    fi

    chunk_clusters.py \\
        --clustering ${clustering} \\
        --sequences ${fasta_name} \\
        --threshold ${size_threshold} \\
        --out_folder chunked_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
        polars: \$(python -c "import importlib.metadata; print(importlib.metadata.version('polars'))")
    END_VERSIONS
    """
}
