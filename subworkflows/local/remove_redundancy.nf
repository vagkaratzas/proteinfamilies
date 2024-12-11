/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { CONCAT_HMMS           } from '../../modules/local/concat_hmms.nf'
// include { REMOVE_REDUNDANT_FAMS } from '../../modules/local/remove_redundant_fams.nf'
include { EXECUTE_CLUSTERING    } from '../../subworkflows/local/execute_clustering'
include { REMOVE_REDUNDANT_SEQS } from '../../modules/local/remove_redundant_seqs.nf'
include { ALIGN_SEQUENCES       } from '../../subworkflows/local/align_sequences'

workflow REMOVE_REDUNDANCY {
    take:
    full_msa // tuple val(meta), path(fas)
    fasta    // tuple val(meta), path(fasta)
    hmm      // tuple val(meta), path(hmm)

    main:
    ch_versions = Channel.empty()

    if (params.remove_family_redundancy) {
        hmm
            .map { meta, model -> [ [id: meta.id], model ] }
            .groupTuple(by: 0)
            .set { ch_hmm }
        CONCAT_HMMS( ch_hmm )
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

        // TODO

    }

    // TODO only for remaining families (after removal of redundant ones)..fasta = ...
    if (params.remove_sequence_redundancy) {
        EXECUTE_CLUSTERING( fasta )
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )
        clustering_tsv = EXECUTE_CLUSTERING.out.clustering_tsv

        // Join together to ensure in sync
        ch_input_for_seq_removal = fasta
            .join(clustering_tsv)
            .multiMap { meta, seqs, clusters ->
                seqs: [ meta, seqs ]
                clusters: [ meta, clusters ]
            }

        REMOVE_REDUNDANT_SEQS( ch_input_for_seq_removal.clusters, ch_input_for_seq_removal.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )

        ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta )
        ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
        full_msa = ALIGN_SEQUENCES.out.alignments
    }

    emit:
    versions = ch_versions
    full_msa = full_msa
}
