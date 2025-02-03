# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2025/02/03]

Initial release of nf-core/proteinfamilies, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Amino acid sequence clustering (mmseqs)
- Multiple sequence alignment (famsa, mafft, clipkit)
- Hidden Markov Model generation (hmmer)
- Between families redundancy removal (hmmer)
- In-family sequence redundancy removal (mmseqs)
- Family updating (hmmer, seqkit, mmseqs, famsa, mafft, clipkit)
- Family statistics presentation (multiqc)
