/*
    FAMILY MODEL GENERATION
*/

include { ALIGN_SEQUENCES  } from '../../subworkflows/local/align_sequences'
include { MAFFT_ALIGN      } from '../../modules/nf-core/mafft/align/main'
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
    ch_versions   = Channel.empty()
    ch_alignments = Channel.empty()

    fasta_chunks
        .transpose()
        .map { meta, file ->
            def baseName = file.toString().split('/')[-1].split('\\.')[0]
            [ [id: meta.id, chunk: baseName], file ]
        }
        .set { msa_input_ch }

    ALIGN_SEQUENCES( msa_input_ch )
    ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )

    if (params.trim_seed_msa) {
        if (params.clipping_tool == 'clipkit') {
            CLIPKIT( ALIGN_SEQUENCES.out.alignments )
            ch_versions = ch_versions.mix( CLIPKIT.out.versions )
            ch_alignments = CLIPKIT.out.clipkit
        } else { // fallback: local module clip_ends
            CLIP_ENDS( ALIGN_SEQUENCES.out.alignments, params.gap_threshold )
            ch_versions = ch_versions.mix( CLIP_ENDS.out.versions )
            ch_alignments = CLIP_ENDS.out.fas
        }
    }

    HMMER_HMMBUILD( ch_alignments, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )

    // Combine with same id to ensure in sync
    HMMER_HMMBUILD.out.hmm
        .map { meta, hmm -> [ [id: meta.id], meta, hmm ] }
        .combine(sequences, by: 0)
        .map { id, meta, hmm, seqs -> [ meta, hmm, seqs, true, params.hmmsearch_write_target, params.hmmsearch_write_domain ] } // write_align must always be true
        .set { ch_input_for_hmmsearch }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    FILTER_RECRUITED( HMMER_HMMSEARCH.out.alignments, HMMER_HMMSEARCH.out.domain_summary, params.hmmsearch_query_length_threshold )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )
    ch_full_msa = FILTER_RECRUITED.out.full_msa
    ch_fasta    = FILTER_RECRUITED.out.fasta

    emit:
    versions = ch_versions
    full_msa = ch_full_msa
    fasta    = ch_fasta
}
