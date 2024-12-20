/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                                    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PROCESS_FAMILIES as PROCESS_FAMILIES_NEW            } from '../subworkflows/local/process_families'
include { PROCESS_FAMILIES as PROCESS_FAMILIES_UPDATE         } from '../subworkflows/local/process_families'

//
// MODULE: Local to the pipeline
//
include { HMMER_HMMSEARCH                                     } from '../modules/nf-core/hmmer/hmmsearch/main'
include { CAT_CAT as CAT_HMMS                                 } from '../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_MSAS                                 } from '../modules/nf-core/cat/cat/main'
include { SEQKIT_SEQ                                          } from '../modules/nf-core/seqkit/seq/main'
include { EXTRACT_FAMILY_REPS as EXTRACT_FAMILY_REPS_CREATION } from '../modules/local/extract_family_reps'
include { EXTRACT_FAMILY_REPS as EXTRACT_FAMILY_REPS_UPDATE   } from '../modules/local/extract_family_reps'
include { GET_NONHITS_SEQS                                    } from '../modules/local/get_nonhits_seqs'
include { FILTER_RECRUITED                                    } from '../modules/local/filter_recruited'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEINFAMILIES {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet_for_create = Channel.empty()
    ch_samplesheet_for_update = Channel.empty()


    ch_samplesheet
        .branch { _meta, _fasta, existing_hmms, existing_msas ->
            to_create: !existing_hmms?.size() && !existing_msas?.size()
            to_update: existing_hmms?.size() && existing_msas?.size()
        }
        .set {
            branch_result
        }

    ch_family_reps = Channel.empty()

    /************************************/
    /* Splitting the samplesheet into 2 */
    /* - Entries to create new families */
    /*   (they only have sequences)     */
    /* - Entries to update existing     */
    /*   families (existing HMM models  */
    /*   and MSAs)                      */
    /************************************/
    ch_samplesheet_for_create = branch_result.to_create.map { meta, fasta, _existing_hmms, _existing_msas -> [meta, fasta] }
    ch_samplesheet_for_update = branch_result.to_update

    /**********************/
    /* UPDATE SUBWORKFLOW */
    /**********************/
    if (branch_result.to_update) {

        // Squeeze the hmm models into a single file
        CAT_HMMS(
            ch_samplesheet_for_update.map { meta, _fasta, existing_hmms, _existing_msas ->
                [meta, file("${existing_hmms}/*.{hmm,hmm.gz}")]
            }
        )

        // Prep the sequences to search them against the HMM concatenated model of families //
        ch_samplesheet_for_update = ch_samplesheet_for_update
            .join(
                CAT_HMMS.out.file_out
            )
            .map { meta, fasta, _existing_hmms, existing_msas, concatenated_hmm ->
                [meta, fasta, concatenated_hmm, existing_msas]
            }

        HMMER_HMMSEARCH(
            ch_samplesheet_for_update.map { meta, fasta, concatenated_hmm, _existing_msas ->
                [meta, concatenated_hmm, fasta, false, false, true]
            }
        )

        // We only keep those sequences that match the HMM model with the families to update
        GET_NONHITS_SEQS(
            ch_samplesheet_for_update.map { meta, fasta, _concatenated_hmm, _existing_msas -> [meta, fasta] }.join(HMMER_HMMSEARCH.out.output)
        )
        ch_versions = ch_versions.mix(GET_NONHITS_SEQS.out.versions)

        FILTER_RECRUITED(
            HMMER_HMMSEARCH.out.alignments,
            HMMER_HMMSEARCH.out.domain_summary,
            params.hmmsearch_query_length_threshold,
        )
        ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

        // Get the sequences from the MSAs
        CAT_MSAS(
            ch_samplesheet_for_update.map { meta, _fasta, _concatenated_hmm, existing_msas ->
                [meta, file("${existing_msas}/*.{aln,aln.gz}")]
            }
        )
        ch_versions = ch_versions.mix(CAT_MSAS.out.versions)

        // An .aln is basically a fasta file with gaps, seqkit - seq will remove the gaps to create a traditional fasta
        SEQKIT_SEQ(CAT_MSAS.out.file_out)
        ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

        // Process existing families
        PROCESS_FAMILIES_UPDATE(
            SEQKIT_SEQ.out.fastx
        )
        ch_versions = ch_versions.mix(PROCESS_FAMILIES_UPDATE.out.versions)

        ch_family_reps = ch_family_reps.mix(PROCESS_FAMILIES_UPDATE.out.family_reps)

        // Process new families
        ch_samplesheet_for_create = ch_samplesheet_for_create.mix(GET_NONHITS_SEQS.out.fasta)
    }

    PROCESS_FAMILIES_NEW(
        ch_samplesheet_for_create
    )
    ch_versions = ch_versions.mix(PROCESS_FAMILIES_NEW.out.versions)

    ch_family_reps = ch_family_reps.mix(PROCESS_FAMILIES_NEW.out.family_reps)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'proteinfamilies_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_family_reps.collect { it[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
