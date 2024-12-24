include { CAT_CAT              } from '../../modules/nf-core/cat/cat/main'
include { HMMER_HMMSEARCH      } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { BRANCH_HITS_FASTA    } from '../../modules/local/branch_hits_fasta'
include { SEQKIT_SEQ           } from '../../modules/nf-core/seqkit/seq/main'
include { CAT_CAT as CAT_FASTA } from '../../modules/nf-core/cat/cat/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update // channel: [meta, sequences, existing_hmms_to_update, existing_msas_to_update]

    main:
    ch_versions            = Channel.empty()
    ch_updated_family_reps = Channel.empty()
    ch_no_hit_seqs         = Channel.empty()

    // TODO: check that the HMMs and the MSAs match //

    // Squeeze the HMMs into a single file
    CAT_CAT(
        ch_samplesheet_for_update.map { meta, _fasta, existing_hmms_to_update, _existing_msas_to_update ->
            [meta, file("${existing_hmms_to_update}/*.{hmm,hmm.gz}")]
        }
    )
    ch_versions = ch_versions.mix( CAT_CAT.out.versions )

    // Prep the sequences to search against the HMM concatenated model of families
    CAT_CAT.out.file_out
        .combine(ch_samplesheet_for_update, by: 0)
        .map { meta, concatenated_hmm, fasta, _existing_hmms_to_update, _existing_msas_to_update -> [meta, concatenated_hmm, fasta, false, false, true] }
        .set { ch_input_for_hmmsearch }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    ch_input_for_branch_hits = HMMER_HMMSEARCH.out.domain_summary
        .join(ch_samplesheet_for_update)
        .multiMap { meta, domtbl, fasta, _existing_hmms_to_update, _existing_msas_to_update ->
            domtbl: [ meta, domtbl ]
            fasta: [ meta, fasta ]
        }

    // Branch hit families/fasta proteins from non hit fasta proteins
    BRANCH_HITS_FASTA ( ch_input_for_branch_hits.fasta, ch_input_for_branch_hits.domtbl, params.hmmsearch_query_length_threshold )
    ch_versions = ch_versions.mix( BRANCH_HITS_FASTA.out.versions )
    ch_no_hit_seqs = BRANCH_HITS_FASTA.out.non_hit_fasta

    BRANCH_HITS_FASTA.out.hits
        .transpose()
        .map { meta, file ->
            def baseName = file.toString().split('/')[-1].split('\\.')[0]
            [[id: meta.id, family: baseName], file]
        }
        .set { hits_fasta }

    ch_samplesheet_for_update
        .map { meta, _fasta, _existing_hmms_to_update, existing_msas_to_update ->
            [meta, file("${existing_msas_to_update}/*")]
        }
        .transpose()
        .map { meta, file ->
            def baseName = file.toString().split('/')[-1].split('\\.')[0]
            [[id: meta.id, family: baseName], file]
        }
        .set { family_msas }

    // Keep fasta with family sequences by removing gaps
    SEQKIT_SEQ( family_msas )
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // Match newly recruited sequences with existing ones for each family
    hits_fasta
        .combine(SEQKIT_SEQ.out.fastx, by: 0)
        .map { meta, new_fasta, family_fasta ->
            [meta, [new_fasta, family_fasta]]
        }
        .set { ch_input_for_cat }

    // Aggregate each family's MSA sequences with the newly recruited ones
    CAT_FASTA( ch_input_for_cat )
    ch_versions = ch_versions.mix( CAT_FASTA.out.versions )
    CAT_FASTA.out.file_out.view()

    emit:
    versions            = ch_versions
    no_hit_seqs         = ch_no_hit_seqs
    updated_family_reps = ch_updated_family_reps
}
