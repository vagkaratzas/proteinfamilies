/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { MMSEQS_CREATEDB       } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER        } from '../../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_LINCLUST       } from '../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV      } from '../../modules/nf-core/mmseqs/createtsv/main'
// include { REMOVE_REDUNDANT_SEQS } from '../../modules/local/remove_redundant_seqs.nf'
// include { REMOVE_REDUNDANT_FAMS } from '../../modules/local/remove_redundant_fams.nf'

workflow REMOVE_REDUNDANCY {
    take:
    full_msa // tuple val(meta), path(fasta)

    main:
    ch_versions = Channel.empty()

    // TODO
    full_msa.view()

    emit:
    versions = ch_versions
    full_msa = full_msa
}
