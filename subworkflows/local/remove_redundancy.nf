/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { EXTRACT_FAMILY_REPS   } from '../../modules/local/extract_family_reps.nf'
include { CONCAT_HMMS           } from '../../modules/local/concat_hmms.nf'
include { HMMER_HMMSEARCH       } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { REMOVE_REDUNDANT_FAMS } from '../../modules/local/remove_redundant_fams.nf'
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
        full_msa
            .map { meta, aln -> [ [id: meta.id], aln ] }
            .groupTuple(by: 0)
            .set { ch_full_msa }
        EXTRACT_FAMILY_REPS( ch_full_msa )
        ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )

        hmm
            .map { meta, model -> [ [id: meta.id], model ] }
            .groupTuple(by: 0)
            .set { ch_hmm }
        CONCAT_HMMS( ch_hmm )
        ch_versions = ch_versions.mix( CONCAT_HMMS.out.versions )

        CONCAT_HMMS.out.hmm
            .combine(EXTRACT_FAMILY_REPS.out.fasta, by: 0)
            .map { meta, model, seqs -> [ meta, model, seqs, false, false, true ] }
            .set { ch_input_for_hmmsearch }

        HMMER_HMMSEARCH( ch_input_for_hmmsearch )
        ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

        fasta
            .map{ meta, fas -> [ [id: meta.id], fas ] }
            .groupTuple(by: 0)
            .set { fasta }

        // Join together to ensure in sync
        ch_input_for_fam_removal = EXTRACT_FAMILY_REPS.out.map
            .join(HMMER_HMMSEARCH.out.domain_summary)
            .join(fasta)
            .multiMap { meta, map, domtbl, seqs ->
                map: [ meta, map ]
                domtbl: [ meta, domtbl ]
                seqs: [ meta, seqs ]
            }

        REMOVE_REDUNDANT_FAMS( ch_input_for_fam_removal.map, ch_input_for_fam_removal.domtbl, ch_input_for_fam_removal.seqs, params.hmmsearch_family_length_threshold )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_FAMS.out.versions )
        fasta = REMOVE_REDUNDANT_FAMS.out.fasta
        fasta
            .transpose()
            .map { meta, file ->
                def baseName = file.toString().split('/')[-1].split('\\.')[0].split('_')[-1]
                [ [id: meta.id, chunk: baseName], file ]
            }
            .set { fasta }
    }

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
