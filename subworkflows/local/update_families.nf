include { HMMER_HMMSEARCH                                     } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { CAT_CAT as CAT_HMMS                                 } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_MSAS                                 } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_SEQS                                 } from '../../modules/nf-core/cat/cat/main'
include { SEQKIT_SEQ                                          } from '../../modules/nf-core/seqkit/seq/main'
include { EXTRACT_FAMILY_REPS as EXTRACT_FAMILY_REPS_CREATION } from '../../modules/local/extract_family_reps'
include { EXTRACT_FAMILY_REPS as EXTRACT_FAMILY_REPS_UPDATE   } from '../../modules/local/extract_family_reps'
include { GET_NONHITS_SEQS                                    } from '../../modules/local/get_nonhits_seqs'
include { FILTER_RECRUITED                                    } from '../../modules/local/filter_recruited'
include { PROCESS_FAMILIES                                    } from '../../subworkflows/local/process_families'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update // channel: [meta, sequences, existing_hmms_to_update, existing_msas_to_update]

    main:
    ch_versions            = Channel.empty()
    ch_updated_family_reps = Channel.empty()
    ch_no_hit_seqs         = Channel.empty()

    // TODO: check that the HMMs and the MSAs match //

    // Squeeze the hmm models into a single file
    CAT_HMMS(
        ch_samplesheet_for_update.map { meta, _fasta, existing_hmms_to_update, _existing_msas ->
            [meta, file("${existing_hmms_to_update}/*.{hmm,hmm.gz}")]
        }
    )

    // Prep the sequences to search them against the HMM concatenated model of families //
    ch_samplesheet_for_update = ch_samplesheet_for_update
        .join(
            CAT_HMMS.out.file_out
        )
        .map { meta, fasta, _existing_hmms, existing_msas_to_update, concatenated_hmm ->
            [meta, fasta, concatenated_hmm, existing_msas_to_update]
        }

    HMMER_HMMSEARCH(
        ch_samplesheet_for_update.map { meta, fasta, concatenated_hmm, _existing_msas ->
            [meta, concatenated_hmm, fasta, false, false, true]
        }
    )

    // We only keep those sequences that match the HMM models with the families to update, the seqs that don't match
    // are sent to the "create families mode"
    GET_NONHITS_SEQS( ch_samplesheet_for_update.map { meta, fasta, _concatenated_hmm, _existing_msas -> [meta, fasta] }.join(HMMER_HMMSEARCH.out.output) )
    ch_versions = ch_versions.mix( GET_NONHITS_SEQS.out.versions )
    ch_no_hit_seqs = GET_NONHITS_SEQS.out.fasta

    FILTER_RECRUITED( HMMER_HMMSEARCH.out.alignments, HMMER_HMMSEARCH.out.domain_summary, params.hmmsearch_query_length_threshold )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    // Get the sequences from the MSAs
    CAT_MSAS(
        ch_samplesheet_for_update.map { meta, _fasta, _concatenated_hmm, existing_msas_to_update ->
            [meta, file("${existing_msas_to_update}/*.{aln,aln.gz}")]
        }
    )
    ch_versions = ch_versions.mix(CAT_MSAS.out.versions)

    // An .aln is basically a fasta file with gaps, seqkit - seq will remove the gaps to create a traditional fasta
    SEQKIT_SEQ( CAT_MSAS.out.file_out )
    ch_versions = ch_versions.mix( SEQKIT_SEQ.out.versions )

    // Add the sequences that match with the seqs from the MSAs
    // TOOD: this is untested
    CAT_SEQS(
        SEQKIT_SEQ.out.fastx.join(
            FILTER_RECRUITED.out.fasta,
            by: [0]
        ).map { meta, msas_fasta, filtered_fastas ->
            [
                [meta, [msas_fasta, filtered_fastas]]
            ]
        }
    )
    ch_versions = ch_versions.mix(CAT_SEQS.out.versions)

    // Process-update existing families
    PROCESS_FAMILIES( CAT_SEQS.out.file_out )
    ch_versions = ch_versions.mix( PROCESS_FAMILIES.out.versions )
    ch_updated_family_reps = PROCESS_FAMILIES.out.family_reps

    emit:
    versions            = ch_versions
    no_hit_seqs         = ch_no_hit_seqs
    updated_family_reps = ch_updated_family_reps
}
