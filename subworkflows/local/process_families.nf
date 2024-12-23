include { CHUNK_CLUSTERS      } from '../../modules/local/chunk_clusters.nf'
include { EXTRACT_FAMILY_REPS } from '../../modules/local/extract_family_reps'

include { EXECUTE_CLUSTERING  } from './execute_clustering'
include { GENERATE_FAMILIES   } from './generate_families'
include { REMOVE_REDUNDANCY   } from './remove_redundancy'


workflow PROCESS_FAMILIES {
    take:
    ch_sequences // channel: [meta, sequences]

    main:
    ch_versions = Channel.empty()

    // Avoid clustering or generating families with less sequences than params.cluster_size_threshold
    // TODO: review this
    ch_sequences = ch_sequences.filter { _map, seqs ->
        {
            seqs.countFasta() >= params.cluster_size_threshold
        }
    }

    // Clustering
    EXECUTE_CLUSTERING(ch_sequences)
    ch_versions = ch_versions.mix(EXECUTE_CLUSTERING.out.versions)

    // Join together to ensure in sync
    ch_input_for_cluster_chunking = ch_sequences
        .join(EXECUTE_CLUSTERING.out.clustering_tsv)
        .multiMap { meta, seqs, clusters ->
            seqs: [meta, seqs]
            clusters: [meta, clusters]
        }

    CHUNK_CLUSTERS(
        ch_input_for_cluster_chunking.clusters,
        ch_input_for_cluster_chunking.seqs,
        params.cluster_size_threshold,
    )
    ch_versions = ch_versions.mix(CHUNK_CLUSTERS.out.versions)

    // Multiple sequence alignment and family generation
    GENERATE_FAMILIES(
        ch_sequences,
        CHUNK_CLUSTERS.out.fasta_chunks,
    )
    ch_versions = ch_versions.mix(GENERATE_FAMILIES.out.versions)

    // Remove redundant sequences and families
    REMOVE_REDUNDANCY( GENERATE_FAMILIES.out.msa, GENERATE_FAMILIES.out.fasta, GENERATE_FAMILIES.out.hmm )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    // Post-processing
    REMOVE_REDUNDANCY.out.msa
        .map { meta, aln -> [ [id: meta.id], aln ] }
        .groupTuple(by: 0)
        .set { ch_msa }

    EXTRACT_FAMILY_REPS(ch_msa)
    ch_versions = ch_versions.mix(EXTRACT_FAMILY_REPS.out.versions)

    emit:
    family_reps = EXTRACT_FAMILY_REPS.out.map
    versions    = ch_versions
}
