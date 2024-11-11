process CHUNK_CLUSTERS {

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(clustering)
    val(size_threshold)

    output:
    path("chunked_clusters"), emit: chunked_clusters
    path "versions.yml"     , emit: versions

    script:
    """
    combine_tables.py --clustering ${clustering} \
        --threshold ${size_threshold} \
        --out_folder chunked_clusters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
