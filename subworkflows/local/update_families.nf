include { UNTAR as UNTAR_HMM    } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_MSA    } from '../../modules/nf-core/untar/main'
include { CAT_CAT as CAT_HMM    } from '../../modules/nf-core/cat/cat/main'
include { HMMER_HMMSEARCH       } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { BRANCH_HITS_FASTA     } from '../../modules/local/branch_hits_fasta'
include { SEQKIT_SEQ            } from '../../modules/nf-core/seqkit/seq/main'
include { CAT_CAT as CAT_FASTA  } from '../../modules/nf-core/cat/cat/main'
include { EXECUTE_CLUSTERING    } from '../../subworkflows/local/execute_clustering'
include { REMOVE_REDUNDANT_SEQS } from '../../modules/local/remove_redundant_seqs.nf'
include { ALIGN_SEQUENCES       } from '../../subworkflows/local/align_sequences'
include { CLIPKIT               } from '../../modules/nf-core/clipkit/main'
include { CLIP_ENDS             } from '../../modules/local/clip_ends.nf'
include { HMMER_HMMBUILD        } from '../../modules/nf-core/hmmer/hmmbuild/main'
include { EXTRACT_FAMILY_REPS   } from '../../modules/local/extract_family_reps.nf'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update // channel: [meta, sequences, existing_hmms_to_update, existing_msas_to_update]

    main:
    ch_versions            = Channel.empty()
    ch_updated_family_reps = Channel.empty()
    ch_no_hit_seqs         = Channel.empty()


    ch_samplesheet_for_update
        .multiMap { meta, _fasta, existing_hmms_to_update, existing_msas_to_update ->
            hmm: [ meta, existing_hmms_to_update ]
            msa: [ meta, existing_msas_to_update ]
        }
        .set { ch_input_for_untar }

    UNTAR_HMM( ch_input_for_untar.hmm )
    ch_versions = ch_versions.mix( UNTAR_HMM.out.versions )

    UNTAR_MSA( ch_input_for_untar.msa )
    ch_versions = ch_versions.mix( UNTAR_MSA.out.versions )

    // TODO: check that the HMMs and the MSAs match //

    // Squeeze the HMMs into a single file
    CAT_HMM( UNTAR_HMM.out.untar.map { meta, folder -> [meta, file("$folder/*")] } )
    ch_versions = ch_versions.mix( CAT_HMM.out.versions )

    // Prep the sequences to search against the HMM concatenated model of families
    CAT_HMM.out.file_out
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

    UNTAR_MSA.out.untar
        .map { meta, folder ->
            [meta, file("$folder/*")]
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
    SEQKIT_SEQ.out.fastx
        .combine(hits_fasta, by: 0)
        .map { meta, family_fasta, new_fasta  ->
            [meta, [family_fasta, new_fasta]]
        }
        .set { ch_input_for_cat }

    // Aggregate each family's MSA sequences with the newly recruited ones
    CAT_FASTA( ch_input_for_cat )
    ch_versions = ch_versions.mix( CAT_FASTA.out.versions )
    fasta_ch = CAT_FASTA.out.file_out

    if (params.remove_sequence_redundancy) {
        // Strict clustering to remove redundancy
        EXECUTE_CLUSTERING( fasta_ch )
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

        // Join together to ensure in sync
        ch_input_for_seq_removal = fasta_ch
            .join(EXECUTE_CLUSTERING.out.clustering_tsv)
            .multiMap { meta, seqs, clusters ->
                seqs: [meta, seqs]
                clusters: [meta, clusters]
            }

        REMOVE_REDUNDANT_SEQS( ch_input_for_seq_removal.clusters, ch_input_for_seq_removal.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )
        fasta_ch = REMOVE_REDUNDANT_SEQS.out.fasta
    }

    ALIGN_SEQUENCES( fasta_ch )
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

    ch_msa
        .map { meta, aln -> [ [id: meta.id], aln ] }
        .groupTuple(by: 0)
        .set { ch_msa }
    EXTRACT_FAMILY_REPS( ch_msa )
    ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )
    ch_updated_family_reps = ch_updated_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    emit:
    versions            = ch_versions
    no_hit_seqs         = ch_no_hit_seqs
    updated_family_reps = ch_updated_family_reps
}
