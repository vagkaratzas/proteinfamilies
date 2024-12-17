/*
    FAMILY MODEL GENERATION
*/

include { ALIGN_SEQUENCES  } from '../../subworkflows/local/align_sequences'
include { CLIPKIT          } from '../../modules/nf-core/clipkit/main'
include { CLIP_ENDS        } from '../../modules/local/clip_ends.nf'
include { HMMER_HMMBUILD   } from '../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMSEARCH  } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { FILTER_RECRUITED } from '../../modules/local/filter_recruited.nf'

workflow GENERATE_FAMILIES {
    take:
    sequences // tuple val(meta), path(fasta)
    fasta_chunks

    main:
    ch_versions = Channel.empty()
    ch_msa      = Channel.empty()
    ch_fasta    = Channel.empty()
    ch_hmm      = Channel.empty()

    fasta_chunks
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path).baseName], file_path ]
        }
        .set { ch_fasta }

    ALIGN_SEQUENCES( ch_fasta )
    ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    ch_msa = ALIGN_SEQUENCES.out.alignments

    if (params.trim_seed_msa) {
        if (params.clipping_tool == 'clipkit') {
            CLIPKIT( ch_msa )
            ch_versions = ch_versions.mix( CLIPKIT.out.versions )
            ch_msa = CLIPKIT.out.clipkit
        } else { // fallback: local module clip_ends
            CLIP_ENDS( ch_msa, params.gap_threshold )
            ch_versions = ch_versions.mix( CLIP_ENDS.out.versions )
            ch_msa = CLIP_ENDS.out.fas
        }
    }

    HMMER_HMMBUILD( ch_msa, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )
    ch_hmm = HMMER_HMMBUILD.out.hmm

    // Combine with same id to ensure in sync
    ch_hmm
        .map { meta, hmm -> [ [id: meta.id], meta, hmm ] }
        .combine(sequences, by: 0)
        .map { id, meta, hmm, seqs -> [ meta, hmm, seqs, true, params.hmmsearch_write_target, params.hmmsearch_write_domain ] } // write_align must always be true
        .set { ch_input_for_hmmsearch }

    if (params.recruit_sequences_with_models) {
        HMMER_HMMSEARCH( ch_input_for_hmmsearch )
        ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

        FILTER_RECRUITED( HMMER_HMMSEARCH.out.alignments, HMMER_HMMSEARCH.out.domain_summary, params.hmmsearch_query_length_threshold )
        ch_versions = ch_versions.mix( FILTER_RECRUITED.out.versions )
        ch_msa = FILTER_RECRUITED.out.full_msa
        ch_fasta = FILTER_RECRUITED.out.fasta
    }

    emit:
    versions = ch_versions
    msa      = ch_msa
    fasta    = ch_fasta
    hmm      = ch_hmm
}
