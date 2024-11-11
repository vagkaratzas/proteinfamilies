/*
    SEQUENCE CLUSTERING
*/

include { MMSEQS_CREATEDB  } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER   } from '../../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_LINCLUST  } from '../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV } from '../../modules/nf-core/mmseqs/createtsv/main'
include { CHUNK_CLUSTERS   } from '../../modules/local/chunk_clusters.nf'

workflow EXECUTE_CLUSTERING {
    take:
    sequences // tuple val(meta), path(fasta)

    main:
    ch_versions       = Channel.empty()
    ch_clustering_tsv = Channel.empty()

    MMSEQS_CREATEDB( sequences )
    ch_versions = ch_versions.mix( MMSEQS_CREATEDB.out.versions )

    if (params.clustering_tool == 'cluster') {
        cluster_res = MMSEQS_CLUSTER( MMSEQS_CREATEDB.out.db )
    } else { // fallback: linclust
        cluster_res = MMSEQS_LINCLUST( MMSEQS_CREATEDB.out.db )
    }
    ch_versions = ch_versions.mix( cluster_res.versions )

    // Join together to ensure in sync
    ch_input_for_createtsv = MMSEQS_CREATEDB.out.db
        .join(cluster_res.db_cluster)
        .multiMap { meta, db, db_cluster ->
            db: [ meta, db ]
            db_cluster: [ meta, db_cluster ]
        }

    MMSEQS_CREATETSV(ch_input_for_createtsv.db_cluster, ch_input_for_createtsv.db, ch_input_for_createtsv.db)
    ch_versions       = ch_versions.mix( MMSEQS_CREATETSV.out.versions )
    ch_clustering_tsv = MMSEQS_CREATETSV.out.tsv

    CHUNK_CLUSTERS(MMSEQS_CREATETSV.out.tsv, params.cluster_size_threshold)
    ch_versions       = ch_versions.mix( CHUNK_CLUSTERS.out.versions )
    ch_cluster_chunks = CHUNK_CLUSTERS.out.chunked_clusters

    emit:
    versions       = ch_versions
    clustering_tsv = ch_clustering_tsv // channel: [ val(meta), tsv ]
    cluster_chunks = ch_cluster_chunks
}
