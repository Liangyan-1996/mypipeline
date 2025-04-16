#!/bin/sh

cellranger="/public/home/liunangroup/liangyan/software/cellranger-9.0.0/cellranger"
index="/public/home/liunangroup/liangyan/Genome/10X-genomics/refdata-gex-GRCh38-2024-A/star"

sample="$1"
thread="$2"

mkdir output

$cellranger \
    --id $sample \
    --transcriptome $index \
    --fastqs raw \
    --sample $sample \
    --create-bam=true \
    --nosecondary \
    --include-introns=true \
    --localcores $thread \
    --output-dir output
