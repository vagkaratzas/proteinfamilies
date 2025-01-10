/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { EXTRACT_FAMILY_REPS       } from '../../modules/local/extract_family_reps.nf'
include { CAT_CAT                   } from '../../modules/nf-core/cat/cat'
include { HMMER_HMMSEARCH           } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { REMOVE_REDUNDANT_FAMS     } from '../../modules/local/remove_redundant_fams.nf'
include { FILTER_NON_REDUNDANT_HMMS } from '../../modules/local/filter_non_redundant_hmms.nf'
include { EXECUTE_CLUSTERING        } from '../../subworkflows/local/execute_clustering'
include { REMOVE_REDUNDANT_SEQS     } from '../../modules/local/remove_redundant_seqs.nf'
include { ALIGN_SEQUENCES           } from '../../subworkflows/local/align_sequences'

workflow REMOVE_REDUNDANCY {
    take:
    msa   // tuple val(meta), path(fas)
    fasta // tuple val(meta), path(fasta)
    hmm   // tuple val(meta), path(hmm)

    main:
    ch_versions = Channel.empty()

    if (params.remove_family_redundancy) {
        msa
            .map { meta, aln -> [[id: meta.id], aln] }
            .groupTuple(by: 0)
            .set { ch_msa }
        EXTRACT_FAMILY_REPS( ch_msa )
        ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )

        hmm
            .map { meta, model -> [[id: meta.id], model] }
            .groupTuple(by: 0)
            .set { ch_hmm }
        CAT_CAT( ch_hmm )
        ch_versions = ch_versions.mix( CAT_CAT.out.versions )

        CAT_CAT.out.file_out
            .combine(EXTRACT_FAMILY_REPS.out.fasta, by: 0)
            .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }
            .set { ch_input_for_hmmsearch }

        ch_input_for_hmmsearch.view()

    //     HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    //     ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    //     fasta
    //         .map { meta, fas -> [[id: meta.id], fas] }
    //         .groupTuple(by: 0)
    //         .set { fasta }

    //     // Join to ensure in sync
    //     ch_input_for_fam_removal = EXTRACT_FAMILY_REPS.out.map
    //         .join(HMMER_HMMSEARCH.out.domain_summary)
    //         .join(fasta)
    //         .multiMap { meta, map, domtbl, seqs ->
    //             map: [meta, map]
    //             domtbl: [meta, domtbl]
    //             seqs: [meta, seqs]
    //         }

    //     REMOVE_REDUNDANT_FAMS( ch_input_for_fam_removal.map, ch_input_for_fam_removal.domtbl, ch_input_for_fam_removal.seqs, params.hmmsearch_family_length_threshold )
    //     ch_versions = ch_versions.mix( REMOVE_REDUNDANT_FAMS.out.versions )
    //     fasta = REMOVE_REDUNDANT_FAMS.out.fasta

    //     // Join to ensure in sync
    //     ch_input_for_hmm_filtering = fasta
    //         .join(ch_hmm)
    //         .multiMap { meta, seqs, models ->
    //             seqs: [meta, seqs]
    //             models: [meta, models]
    //         }
    //     FILTER_NON_REDUNDANT_HMMS( ch_input_for_hmm_filtering.seqs, ch_input_for_hmm_filtering.models )
    //     ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_HMMS.out.versions )

    //     fasta
    //         .transpose()
    //         .map { meta, file ->
    //             [[id: meta.id, chunk: file.getSimpleName().split('_')[-1]], file]
    //         }
    //         .set { fasta }
    }

    // if (params.remove_sequence_redundancy) {
    //     EXECUTE_CLUSTERING( fasta )
    //     ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

    //     REMOVE_REDUNDANT_SEQS( EXECUTE_CLUSTERING.out.clusters, EXECUTE_CLUSTERING.out.seqs )
    //     ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )

    //     ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta )
    //     ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    //     msa = ALIGN_SEQUENCES.out.alignments
    // }

    emit:
    versions = ch_versions
    // msa      = msa
}
