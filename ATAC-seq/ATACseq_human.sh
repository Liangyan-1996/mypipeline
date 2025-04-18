fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
bowtie2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bowtie2"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/samtools"
sambamba="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
bedtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bedtools"
macs2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/macs2"
alignmentSieve="/public/home/liunangroup/liangyan/software/miniconda3/bin/alignmentSieve"
computeMatrix="/public/home/liunangroup/liangyan/software/miniconda3/bin/computeMatrix"
plotHeatmap="/public/home/liunangroup/liangyan/software/miniconda3/bin/plotHeatmap"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"
ScriptDir="/public/home/liunangroup/liangyan/pipeline/mypipeline/script"
Rscript="/public/home/liunangroup/liangyan/software/miniconda3/envs/ATAC/bin/Rscript"

bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/hg38/bowtie2_index/genome"
blacklist="/public/home/liunangroup/liangyan/pipeline/CUT-RUNTools-2.0/blacklist/hg38.blacklist.bed"
GenomeSize="2913022398"
specie="hs"
ucsc_anno="TxDb.Hsapiens.UCSC.hg38.knownGene"
TSS="/public/home/liunangroup/liangyan/Genome/gencode/hg38/hg38.TSS.bed"

sample="$1"
thread="$2"

# ----------------------

workdir=`pwd`
library="0"

cd $workdir
mkdir qc
mkdir mapping
mkdir log
mkdir peak

if [[ -e raw/${sample}_1.fastq.gz ]] && [[ -e raw/${sample}_2.fastq.gz ]]; then

library="2"
echo "PE data"
$fastp \
    -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz \
    -o qc/${sample}.fastp1.fq.gz -O qc/${sample}.fastp2.fq.gz \
    -w ${thread} 2> $workdir/log/${sample}.log
$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    -X 2000 --no-mixed --no-discordant --no-unal \
    -1 qc/${sample}.fastp1.fq.gz \
    -2 qc/${sample}.fastp2.fq.gz 2> $workdir/log/${sample}.log |\
    $samtools view -@ $thread -bS -q 30 -o mapping/${sample}.bam

elif [[ -e raw/${sample}.fastq.gz ]]; then

library="1"
echo "SE data"
$fastp -i raw/${sample}.fastq.gz -o qc/${sample}.fastp.fq.gz -w ${thread} -A
$bowtie2 \
    --end-to-end --very-sensitive -x $bowtie2index -p ${thread} \
    --no-unal -U qc/${sample}.fastp.fq.gz 2> $workdir/log/${sample}.log |\
    $samtools view -@ $thread -bS -q 30 -o mapping/${sample}.bam

else

echo 'Error! There is no raw data'

fi

cd mapping

# sorting
$sambamba sort -t $thread -m 10G ${sample}.bam
# remove MT reads
$samtools view -@ $thread -h ${sample}.sorted.bam |\
    awk 'BEGIN{FS="\t"}{if($3!="chrM"){print $0}}' |\
    $samtools view -@ $thread -bS > ${sample}.rmMT.bam
# remove PCR duplicates
$sambamba markdup -r -t ${thread} ${sample}.rmMT.bam ${sample}.dedup.bam 2>> $workdir/log/${sample}.log
# remove blacklist reads
$bedtools intersect -nonamecheck -v -a ${sample}.dedup.bam -b $blacklist > ${sample}.rmblack.bam
$samtools index -@ $thread ${sample}.rmblack.bam
# shift
$alignmentSieve --numberOfProcessors $thread --ATACshift -b $sample.rmblack.bam -o $sample.shift.bam
# final sorting
$sambamba sort -t $thread -m 10G -o $sample.clean.bam $sample.shift.bam
# TSS enrichment score
$Rscript $ScriptDir/ATACseq/TSS_enrich_score.R ${sample}.clean.bam ${ucsc_anno} >> $workdir/log/${sample}.log
# length distribution
$samtools view $sample.clean.bam |\
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' |\
    sort | uniq | cut -f2 > ${sample}.length.txt
$Rscript $ScriptDir/ATACseq/frag_len.R ${sample}.length.txt ${sample}.length.png
# coverage
$bamCoverage \
    -b ${sample}.clean.bam -o ${sample}.bw -p ${thread} \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize ${GenomeSize} --extendReads 100
# visualization
$computeMatrix reference-point -S ${sample}.bw -R $TSS -a 3000 -b 3000 -p ${thread} -o ${sample}.gz
$plotHeatmap -m ${sample}.gz -out ${sample}.png --dpi 300 --plotFileFormat png --yMin 0
# log
line_align=`$samtools view -@ $thread -c ${sample}.sorted.bam`
line_MT=`$samtools view -@ $thread -c ${sample}.sorted.bam chrM`
MT_pct=`echo -e "scale=5;${line_MT} / ${line_align} * 100"|bc`
echo -e "MT proportion: ${MT_pct}%" >> $workdir/log/${sample}.log
# remove tmp
rm ${sample}.sorted.bam* ${sample}.bam* ${sample}.rmMT.bam* ${sample}.rmblack.bam* ${sample}.dedup.bam* ${sample}.shift.bam*
# peak calling
if [ "$library" == "2" ]; then
$macs2 callpeak \
    -t ${sample}.clean.bam -f BAMPE --gsize $specie -n ${sample} -q 0.01 \
    --keep-dup all --outdir $workdir/peak
elif [ "$library" == "1" ]; then
$macs2 callpeak \
    -t ${sample}.clean.bam -f BAM --gsize $specie -n ${sample} -q 0.01 \
    --shift -70 --extsize 140 --nomodel --keep-dup all --outdir $workdir/peak
fi
# FRiP
line_clean=`$samtools view -c ${sample}.clean.bam`
line_inPeak=`$bedtools sort -i $workdir/peak/${sample}_peaks.narrowPeak |\
    $bedtools intersect -u -nonamecheck -a ${sample}.clean.bam -b stdin -ubam|\
    $samtools view -c -@ ${thread}`
FRiP=`echo -e "scale=5;${line_inPeak} / ${line_clean}"|bc`
echo -e "Clean Fragments: ${line_clean}" >> $workdir/log/${sample}.log
echo -e "FRiP: $FRiP" >> $workdir/log/${sample}.log
