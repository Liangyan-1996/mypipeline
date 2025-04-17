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
specie="hs"

sample="$1"
thread="$2"
workdir=`pwd`

cd $workdir
mkdir qc
mkdir mapping
mkdir logs
mkdir peak

$fastp \
    -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz \
    -o qc/${sample}_1.fq.gz -O qc/${sample}_2.fq.gz \
    -w $thread


$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    -X 700 -I 10 --no-mixed --no-discordant --no-unal \
    -1 qc/${sample}_1.fq.gz \
    -2 qc/${sample}_2.fq.gz 2> $workdir/log/$sample.log |\
    $samtools view -q30 -bS -@ $thread -o mapping/$sample.bam

cd mapping

$sambamba sort -t $thread -m 10G ${sample}.bam

$sambamba markdup -r -t ${thread} ${sample}.sorted.bam ${sample}.dedup.bam 2>> $workdir/log/${sample}.log

$bedtools intersect -nonamecheck -v -a ${sample}.dedup.bam -b $blacklist > ${sample}.rmblack.bam
$samtools index -@ $thread ${sample}.rmblack.bam

$alignmentSieve --numberOfProcessors $thread --ATACshift -b $sample.rmblack.bam -o $sample.shift.bam

mv $sample.shift.bam $sample.clean.bam
$samtools index -@ $thread $sample.clean.bam

$bamCoverage \
    -b ${sample}.clean.bam -o ${sample}.bw -p ${thread} \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize ${genomesize}

line_align=`$samtools view -@ $thread -c ${sample}.sorted.bam`
line_MT=`$samtools view -@ $thread -c ${sample}.sorted.bam chrM`
MT_pct=`echo -e "scale=5;${line_MT} / ${line_align} * 100"|bc`
echo -e "MT proportion: ${MT_pct}%" >> $workdir/log/${sample}.log
# Remove the mitochondrial DNA and calculate its proportion
rm ${sample}.sorted.bam* ${sample}.bam* ${sample}.dedup.bam* ${sample}.rmblack.bam*

$samtools view $sample.clean.bam |\
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' |\
    sort | uniq | cut -f2 > ${sample}.length.txt
$Rscript $ScriptDir/frag_len.R ${sample}.length.txt ${sample}.length.png

$macs2 callpeak \
    -t ${sample}.clean.bam -f BAMPE --gsize $specie \
    -n ${sample}_narrow -q 0.01 --keep-dup all \
    --outdir ${workdir}/peak; 2>> ${workdir}/logs/$sample.log
$macs2 callpeak \
    -t ${sample}.clean.bam -f BAMPE --gsize $specie \
    -n ${sample}_broad -q 0.01 --keep-dup all \
    --broad --broad-cutoff 0.05 \
    --outdir ${workdir}/peak; 2>> ${workdir}/logs/$sample.log
