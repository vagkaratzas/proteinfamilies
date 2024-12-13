# nf-core/proteinfamilies: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

Initial clustering:

- [MMseqs2](#mmseqs2) initial clustering of input amino acid sequences and filtering with membership threshold

Multiple sequence alignment:

- [FAMSA](#famsa) aligner option. Best speed and sensitivity option to build seed multiple sequence alignment for the families.
- [mafft](#mafft) aligner option. Fast but not as sensitive as FAMSA to build seed multiple sequence alignment for the families.
- [ClipKIT](#clipkit) to optionally clip gapped portions of the multiple sequence alignment (MSA)

Generating family models:

- [hmmer](#hmmer) to build the family HMM (hmmbuild) and to 'fish' additional sequences from the input fasta file, with given thresholds, into the family and also build the family full MSA (hmmsearch).

Removing redundancy:

- [hmmer](#hmmer-for-redundancy-removal) to match family representative sequences against other family models in order to keep non redundant ones
- [MMseqs2](#mmseqs2-for-redundancy-removal) to strictly cluster the sequences within each of the remaining families, in order to still capture the evolutionary diversity within a family, but without keeping all the almost identical sequences.
- [FAMSA](#famsa-for-redundancy-removal) aligner option. Re-align full MSA with non-redundant set of sequences.
- [mafft](#mafft-for-redundancy-removal) aligner option. Re-align full MSA with non-redundant set of sequences.

Reporting:

- [Extract family representatives](#extract-family-representatives) to produce the final metadata file along with a fasta of all family representative sequences (can be used downstream for structural prediction).
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### MMseqs2

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
  - `initial_clustering/`
    - `mmseqs_createtsv/`
      - `<samplename>.tsv`: tab-separated table containing 2 columns; the first one with the cluster representative sequences, and the second with the cluster members
    - `mmseqs_createdb/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format db of fasta sequences
    - `mmseqs_linclust/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
    - `mmseqs_cluster/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
    - `filtered_fasta_chunks/`
      - `<samplename>/`
        - `chunked_fasta/`
          - `*.fasta`: (optional) fasta files with amino acid sequences of each cluster above the membership threshold
</details>

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### FAMSA aligner

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `famsa_align/`
    - `<samplename>/`
      - `<samplename>_*.aln`: fasta files with aligned amino acid sequences
</details>

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft aligner

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `mafft_align/`
    - `<samplename>/`
      - `<samplename>_*.fas`: fasta files with aligned amino acid sequences
</details>

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### ClipKIT

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `clipkit/`
    - `<samplename>/`
      - `*.clipkit`: gap-clipped (start, middle, end) fasta files of aligned amino acid sequences
  - `clip_ends/`
    - `<samplename>/`
      - `*.fas`: gap-clipped (only start and end) fasta files of aligned amino acid sequences
</details>

[ClipKIT](https://github.com/JLSteenwyk/ClipKIT) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### hmmer

<details markdown="1">
<summary>Output files</summary>

- `hmmer/`
  - `hmmbuild/`
    - `<samplename>/`
      - `<samplename>_*.hmm.gz`: compressed hmm model for the family
      - `<samplename>_*.hmmbuild.txt`: (optional) hmmbuild execution log
  - `hmmsearch/`
    - `<samplename>/`
      - `<samplename>_*.sto.gz`: (optional) full multiple sequence alignment of the family
      - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results along parameters info
      - `<samplename>_*.txt.gz`: (optional) hmmsearch execution log
- `full_msa/`
  - `pre_non_redundant/`
    - `<samplename>/`
      - `<samplename>_*.fas.gz`: compressed family full MSA (before checking for redundancy)
</details>

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### hmmer for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `remove_redundancy/`
  - `hmmer/`
    - `hmmbuild/`
      - `concatenated/`
        - `<samplename>_*.hmm.gz`: (optional) concatenated compressed hmm model for all families in a given sample (pre redundancy removal)
    - `hmmsearch/`
      - `<samplename>/`
        - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results of family reps against families' HMMs
  - `family_reps/`
    - `<samplename>/`
      - `<samplename>_meta_mqc.csv`: (optional) csv with metadata (Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
      - `<samplename>_reps.fa`: (optional) fasta file of all family representative sequences (one sequence per family)
  - `non_redundant_fams/`
    - `<samplename>/`
      - `non_redundant/`
        - `<samplename>_*.fasta.gz`: (optional) compressed family full MSA (after checking for family redundancy)

</details>

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### MMseqs2 for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
  - `redundancy_clustering/`
    - `mmseqs_createtsv/`
      - `<samplename>.tsv`: tab-separated table containing 2 columns; the first one with the cluster representative sequences, and the second with the cluster members
    - `mmseqs_createdb/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format db of fasta sequences
    - `mmseqs_linclust/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
    - `mmseqs_cluster/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
- `remove_redundancy/`
    - `reps_fasta/`
      - `<samplename>/`
        - `<samplename>_reps.fa`: (optional) fasta file of all family representative sequences (one sequence per family)
</details>

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### FAMSA for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `full_msa/`
  - `non_redundant/`
    - `<samplename>/`
      - `<samplename>_*.aln`: family full MSA (after checking for sequence redundancy)
</details>

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `full_msa/`
  - `non_redundant/`
    - `<samplename>/`
      - `<samplename>_*.fas`: fasta files with aligned amino acid sequences (after checking for sequence redundancy)
</details>

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### Extract family representatives

<details markdown="1">
<summary>Output files</summary>

- `family_reps/`
  - `<samplename>/`
    - `<samplename>_meta_mqc.csv`: csv with metadata to print with MultiQC (Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
    - `<samplename>_reps.fa`: fasta file of all family representative sequences (one sequence per family)
</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

Custom output MultiQC data includes a metadata file (multiqc_data/multiqc_family_metadata.txt) with family information such as: Sample,Family Id,Size,Representative Length,Representative Id,Sequence

This custom metadata is presented as a data table on the web page.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
