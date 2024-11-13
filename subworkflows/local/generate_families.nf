/*
    FAMILY MODEL GENERATION
*/

include { FAMSA_ALIGN     } from '../../modules/nf-core/famsa/align/main'
include { HMMER_HMMBUILD  } from '../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMSEARCH } from '../../modules/nf-core/hmmer/hmmsearch/main'

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

    FAMSA_ALIGN(msa_input_ch, [[:],[]], false)
    ch_versions = ch_versions.mix( FAMSA_ALIGN.out.versions )

    HMMER_HMMBUILD( FAMSA_ALIGN.out.alignment, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )

    // Combine with same id to ensure in sync
    HMMER_HMMBUILD.out.hmm
        .map{ meta, hmm -> [ [id: meta.id], meta, hmm ] }
        .combine(sequences, by: 0)
        .map { id, meta, hmm, seqs -> [ meta, hmm, seqs, true, params.hmmsearch_write_target, params.hmmsearch_write_domain ] } // write_align must always be true
        .set { ch_input_for_hmmsearch }

    HMMER_HMMSEARCH{ ch_input_for_hmmsearch }
    ch_versions   = ch_versions.mix( HMMER_HMMSEARCH.out.versions )
    ch_alignments = HMMER_HMMSEARCH.out.alignments

    emit:
    versions   = ch_versions
    alignments = ch_alignments
}
