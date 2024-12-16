process FILTER_RECRUITED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb3700531c7ec639f59f084ab64c05e881d654dcf829db163539f2f0b095e09d/data' :
        'community.wave.seqera.io/library/biopython:1.84--3318633dad0031e7' }"

    input:
    tuple val(meta) , path(sto)
    tuple val(meta2), path(domtbl)
    val(length_threshold)

    output:
    tuple val(meta), path("*.fas.gz")  , emit: full_msa
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_recruited.py \\
        --alignment ${sto} \\
        --domtbl ${domtbl} \\
        --length_threshold ${length_threshold} \\
        --out_msa ${prefix}.fas.gz \\
        --out_fasta ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
