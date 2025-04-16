#!/bin/bash

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
bowtie2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bowtie2"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/samtools"
sambamba="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
bedtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bedtools"
macs2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/macs2"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"

bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/hg38/bowtie2_index/genome"
blacklist="/public/home/liunangroup/liangyan/pipeline/CUT-RUNTools-2.0/blacklist/hg38.blacklist.bed"
GenomeSize="2913022398"

sample="$1"
thread="$2"

# ------------------------

workdir=`pwd`

cd $workdir
mkdir qc
mkdir mapping
mkdir log

if [[ -e raw/${sample}_1.fastq.gz ]] && [[ -e raw/${sample}_2.fastq.gz ]]; then

$fastp \
    -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz \
    -o qc/${sample}_1.fq.gz -O qc/${sample}_2.fq.gz \
    -w $thread

$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    -X 2000 --no-mixed --no-discordant --no-unal \
    -1 qc/${sample}_1.fq.gz \
    -2 qc/${sample}_2.fq.gz 2> $workdir/log/$sample.log |\
    awk 'BEGIN{FS="\t"}{if($0~"^@"){print $0}else if($3!~"chrM"){print $0}}' |\
    $samtools view -q30 -bS -@ $thread -o mapping/$sample.bam

elif [[ -e raw/${sample}.fastq.gz ]]; then

$fastp -i raw/${sample}.fastq.gz -o qc/${sample}.fastp.fq.gz -w ${thread} -A
$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    --no-unal -U qc/${sample}.fastp.fq.gz 2> $workdir/log/${sample}.log |\
    awk 'BEGIN{FS="\t"}{if($0~"^@"){print $0}else if($3!~"chrM"){print $0}}' |\
    $samtools view -q30 -bS -@ $thread -o mapping/$sample.bam

else

echo 'Error! There is no raw data'

fi

cd mapping
$sambamba sort -m 10G -t $thread $sample.bam
$sambamba markdup \
    -r -t ${thread} ${sample}.sorted.bam ${sample}.dedup.bam 2>> $workdir/log/${sample}.log
# -r remove duplicates instead of marking them only, -t threads
$bedtools intersect \
    -nonamecheck -v -a ${sample}.dedup.bam -b $blacklist > ${sample}.rmblacklist.bam
$sambamba sort -m 10G -t ${thread} -o ${sample}.clean.bam ${sample}.rmblacklist.bam
$bamCoverage \
    -b ${sample}.clean.bam -o ${sample}.bw -p ${thread} \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize $GenomeSize --extendReads 100
rm ${sample}.bam* ${sample}.sorted.bam* ${sample}.dedup.bam* ${sample}.rmblacklist.bam*

