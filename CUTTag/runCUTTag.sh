#!/bin/bash

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
bowtie2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bowtie2"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/samtools"
sambamba="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
bedtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bedtools"
macs2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/macs2"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"
computeMatrix="/public/home/liunangroup/liangyan/software/miniconda3/bin/computeMatrix"
plotHeatmap="/public/home/liunangroup/liangyan/software/miniconda3/bin/plotHeatmap"
Rscript="/public/home/liunangroup/liangyan/software/miniconda3/envs/ATAC/bin/Rscript"
scriptdir="/public/home/liunangroup/liangyan/pipeline/mypipeline/script"
bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/hg38/bowtie2_index/genome"
blacklist="/public/home/liunangroup/liangyan/pipeline/CUT-RUNTools-2.0/blacklist/hg38.blacklist.bed"
homer="/public/home/liunangroup/liangyan/software/homer/bin/findMotifsGenome.pl"
GenomeSize="2913022398"
specie="hs"

sample="$1"
thread="$2"
workdir=`pwd`

cd $workdir
mkdir qc
mkdir mapping
mkdir log
mkdir peak_narrow peak_broad
mkdir motif
# qc
$fastp \
    -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz \
    -o qc/${sample}_1.fq.gz -O qc/${sample}_2.fq.gz \
    -w $thread
# alignment
$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    -X 700 -I 10 --no-mixed --no-discordant --no-unal \
    -1 qc/${sample}_1.fq.gz \
    -2 qc/${sample}_2.fq.gz 2> $workdir/log/$sample.log |\
    $samtools view -q30 -bS -@ $thread -o mapping/$sample.bam

cd mapping
# sorting
$sambamba sort -t $thread -m 10G ${sample}.bam
# remove PCR duplicates
$sambamba markdup -r -t ${thread} ${sample}.sorted.bam ${sample}.dedup.bam 2>> $workdir/log/${sample}.log
# remove blacklist reads
$bedtools intersect -nonamecheck -v -a ${sample}.dedup.bam -b $blacklist > ${sample}.rmblack.bam
$samtools index -@ $thread ${sample}.rmblack.bam
# shift
$alignmentSieve --numberOfProcessors $thread --ATACshift -b $sample.rmblack.bam -o $sample.shift.bam
# final sorting
$sambamba sort -t $thread -m 10G -o $sample.clean.bam $sample.shift.bam
# coverage
$bamCoverage \
    -b ${sample}.clean.bam -o ${sample}.bw -p ${thread} \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize $GenomeSize
# log
line_align=`$samtools view -@ $thread -c ${sample}.sorted.bam`
line_MT=`$samtools view -@ $thread -c ${sample}.sorted.bam chrM`
MT_pct=`echo -e "scale=5;${line_MT} / ${line_align} * 100"|bc`
echo -e "MT proportion: ${MT_pct}%" >> $workdir/log/${sample}.log
# remove tmp
rm ${sample}.sorted.bam* ${sample}.bam* ${sample}.dedup.bam* ${sample}.rmblack.bam* ${sample}.shift.bam*

# peak calling
$macs2 callpeak \
    -t ${sample}.clean.bam -f BAMPE --gsize $specie \
    -n ${sample} -q 0.01 --keep-dup all \
    --outdir ${workdir}/peak_narrow; 2>> ${workdir}/log/$sample.log
$macs2 callpeak \
    -t ${sample}.clean.bam -f BAMPE --gsize $specie \
    -n ${sample} -q 0.01 --keep-dup all \
    --broad --broad-cutoff 0.05 \
    --outdir ${workdir}/peak_broad; 2>> ${workdir}/log/$sample.log

# motif enrichment
$homer ${workdir}/peak_narrow/${sample}_summits.bed hg38 ${workdir}/motif/$sample -size -100,100 -p $thread

# post qc
# length distribution
$samtools view ${workdir}/mapping/${sample}.clean.bam |\
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' |\
    sort | uniq | cut -f2 > ${workdir}/mapping/${sample}.length.txt
$Rscript $scriptdir/CUTTag/frag_len.R ${workdir}/mapping/${sample}.length.txt ${workdir}/mapping/${sample}.length.png

# heatmap
$computeMatrix reference-point \
    -S ${workdir}/mapping/${sample}.bw \
    -R ${workdir}/peak_narrow/${sample}_summits.bed \
    -o ${workdir}/mapping/${sample}.gz \
    -a 3000 -b 3000 -p ${thread}
$plotHeatmap \
    --matrixFile ${workdir}/mapping/${sample}.gz \
    -out ${workdir}/mapping/${sample}.png \
    --dpi 300 --plotFileFormat png --yMin 0

# FRiP
line_clean=`$samtools view -c $workdir/mapping/${sample}.clean.bam`
line_inPeak=`$bedtools sort -i $workdir/peak_narrow/${sample}_peaks.narrowPeak |\
    $bedtools intersect -u -nonamecheck -a $workdir/mapping/${sample}.clean.bam -b stdin -ubam |\
    $samtools view -c -@ ${thread}`
FRiP=`echo -e "scale=5;${line_inPeak} / ${line_clean}"|bc`
echo -e "Clean Fragments: ${line_clean}" >> $workdir/log/${sample}.log
echo -e "FRiP: $FRiP" >> $workdir/log/${sample}.log
