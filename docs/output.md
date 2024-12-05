# nf-core/proteinfamilies: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

Clustering:

- [MMseqs2](#mmseqs2) clustering of input amino acid sequences
- [Chunked fasta clusters](#chunked-fasta-clusters) derived from the MMseqs2 clusters for parallel downstream processing.

Multiple sequence alignment:

- [FAMSA](#famsa) aligner option. Best speed and sensitivity option to build seed multiple sequence alignment for the families.
- [mafft](#mafft) aligner option. Fast but not as sensitive as FAMSA to build seed multiple sequence alignment for the families.
- [ClipKIT](#clipkit) to optionally clip gapped portions of the MSA (start, middle, end)

Generating family models:

- [hmmer](#hmmer) to build the family HMM (hmmbuild) and to 'fish' additional sequences from the input fasta file into the family and also build the full MSA (hmmsearch).
- [Extract family representatives](#extract-family-representatives) to produce the final metadata file along with a fasta of all family representative sequences (can be used downstream for structural prediction).

Reporting:

- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### MMseqs2

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
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
</details>

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### Chunked fasta clusters

<details markdown="1">
<summary>Output files</summary>

- `fasta_chunks/`
  - `<samplename>/`
    - `chunked_fasta/`
      - `*.fasta`: fasta files with amino acid sequences of each cluster
</details>

### FAMSA aligner

<details markdown="1">
<summary>Output files</summary>

- `famsa_align/`
  - `<samplename>/`
    - `*.aln`: fasta files with aligned amino acid sequences
</details>

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft aligner

<details markdown="1">
<summary>Output files</summary>

- `mafft_align/`
  - `<samplename>/`
    - `*.fas`: fasta files with aligned amino acid sequences
</details>

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### ClipKIT

<details markdown="1">
<summary>Output files</summary>

- `mafft_align/`
  - `<samplename>/`
    - `*.clipkit`: gap-clipped (start, middle, end) fasta files of aligned amino acid sequences
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
      - `<samplename>_*.sto.gz`: full multiple sequence alignment of the family
      - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results along parameters info
      - `<samplename>_*.txt.gz`: (optional) hmmsearch execution log
</details>

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

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
