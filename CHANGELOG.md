# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2025/02/05]

Initial release of nf-core/proteinfamilies, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Amino acid sequence clustering (mmseqs)
- Multiple sequence alignment (famsa, mafft, clipkit)
- Hidden Markov Model generation (hmmer)
- Between families redundancy removal (hmmer)
- In-family sequence redundancy removal (mmseqs)
- Family updating (hmmer, seqkit, mmseqs, famsa, mafft, clipkit)
- Family statistics presentation (multiqc)

## v1.0.1 - [2025/03/26]

Hackathon 2025 - JSON schema improvements

### `Changed`

- slight improvements of 'nextflow_schema.json'
- slight improtmenets of 'assets/schema_input.json'.


