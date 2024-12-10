/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { EXECUTE_CLUSTERING } from '../../subworkflows/local/execute_clustering'
// include { REMOVE_REDUNDANT_SEQS } from '../../modules/local/remove_redundant_seqs.nf'
// include { REMOVE_REDUNDANT_FAMS } from '../../modules/local/remove_redundant_fams.nf'

workflow REMOVE_REDUNDANCY {
    take:
    full_msa // tuple val(meta), path(fas)
    fasta    // tuple val(meta), path(fasta)

    main:
    ch_versions = Channel.empty()

    // TODO
    // if (params.remove_family_redundancy) {
    //
    // }

    // TODO only for remaining families (after removal of redundant ones)
    if (params.remove_sequence_redundancy) {
        fasta = fasta
            .map { meta, aln ->
                def id = meta.id + "_" + meta.chunk
                return[ [id: id], aln ]
            }
        params.cluster_seq_identity = 0.95
        params.cluster_coverage = 0.95
        EXECUTE_CLUSTERING(fasta)
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )
        clustering_tsv = EXECUTE_CLUSTERING.out.clustering_tsv
        // TODO continue with producing MSAs from cluster reps only
        clustering_tsv.view()
    }

    emit:
    versions = ch_versions
    full_msa = full_msa
}
