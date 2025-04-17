# my bioinifomatic pipeline

```
mypipeline
├── ATAC-seq
│   ├── ATACseq_human.sh
│   └── z.sh
├── cellranger
│   └── scRNAseq_10X.sh
├── ChIP-seq
│   └── ChIPseq_human.sh
├── CUTTag
│   ├── runCUTRUN.sh
│   ├── runCUTTag.sh
│   └── z.sh
├── RNAseq
│   ├── RNA-seq_McSplicer
│   │   └── runRNAseq_human.sh
│   ├── RNA-seq_SpliSER
│   │   └── runRNAseq_human.sh
│   ├── runRNAseq_human.sh
│   ├── runRNAseq_human.TotolPipeline.sbatch
│   ├── script
│   └── z.sh
└── script
    ├── ATACseq
    │   ├── frag_len.R
    │   └── TSS_enrich_score.R
    ├── CUTTag
    │   └── frag_len.R
    ├── path.sh
    └── RNAseq
        ├── deseq2 -> star_DE.R
        ├── mergeTPM -> mergeTPM.py
        ├── mergeTPM.py
        └── star_DE.R
```