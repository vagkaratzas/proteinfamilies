/*
    FAMILY MODEL GENERATION
*/

include { FAMSA_ALIGN      } from '../../modules/nf-core/famsa/align/main'
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

    if (params.alignment_tool == 'famsa') {
        alignment_res = FAMSA_ALIGN(msa_input_ch, [[:],[]], false)
        ch_versions = ch_versions.mix( FAMSA_ALIGN.out.versions )
        alignment_ch  = alignment_res.alignment
    } else { // fallback: mafft
        alignment_res = MAFFT_ALIGN( msa_input_ch, [[:], []], [[:], []], [[:], []], [[:], []], [[:], []], false )
        ch_versions = ch_versions.mix( MAFFT_ALIGN.out.versions )
        alignment_ch  = alignment_res.fas
    }

    if (params.trim_seed_msa) {
        if (params.clipping_tool == 'clipkit') {
            CLIPKIT(alignment_ch)
            ch_versions = ch_versions.mix( CLIPKIT.out.versions )
            alignment_ch = CLIPKIT.out.clipkit
        } else { // fallback: local module clip_ends
            CLIP_ENDS(alignment_ch, params.gap_threshold)
            ch_versions = ch_versions.mix( CLIP_ENDS.out.versions )
            alignment_ch = CLIP_ENDS.out.fas
        }
    }

    HMMER_HMMBUILD( alignment_ch, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )

    // Combine with same id to ensure in sync
    HMMER_HMMBUILD.out.hmm
        .map { meta, hmm -> [ [id: meta.id], meta, hmm ] }
        .combine(sequences, by: 0)
        .map { id, meta, hmm, seqs -> [ meta, hmm, seqs, true, params.hmmsearch_write_target, params.hmmsearch_write_domain ] } // write_align must always be true
        .set { ch_input_for_hmmsearch }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions   = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    FILTER_RECRUITED( HMMER_HMMSEARCH.out.alignments, HMMER_HMMSEARCH.out.domain_summary, params.hmmsearch_query_length_threshold )
    ch_versions   = ch_versions.mix( HMMER_HMMSEARCH.out.versions )
    ch_alignments = FILTER_RECRUITED.out.full_msa

    emit:
    versions   = ch_versions
    alignments = ch_alignments
}
