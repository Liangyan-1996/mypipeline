#!/bin/sh

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
STAR="/public/home/liunangroup/liangyan/software/miniconda3/envs/rnaseq/bin/STAR"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/bin/samtools"
stringtie="/public/home/liunangroup/liangyan/software/miniconda3/envs/rnaseq/bin/stringtie"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"
mergeTPM="/public/home/liunangroup/liangyan/pipeline/RNAseq/script/mergeTPM"

STAR_index="/public/home/liunangroup/liangyan/Genome/ensembl/hg38/STAR_index"
GTF="/public/home/liunangroup/liangyan/Genome/ensembl/hg38/Homo_sapiens.GRCh38.110.chr.gtf"

sample="$1"
thread="$2"
library="0"
genomesize="2913022398"
# 2913022398 for hg38
# 2652783500 for mm10

mkdir QC
mkdir mapping
mkdir quantify
mkdir log

if [[ -e raw/${sample}_1.fastq.gz ]] && [[ -e raw/${sample}_2.fastq.gz ]]; then

library="2"
$fastp \
 -i raw/${sample}_1.fastq.gz \
 -I raw/${sample}_2.fastq.gz \
 -o QC/${sample}_1.fq.gz \
 -O QC/${sample}_2.fq.gz \
 -w $thread 2> log/$sample.log

$STAR \
 --twopassMode Basic \
 --runThreadN $thread \
 --genomeDir $STAR_index \
 --readFilesCommand zcat \
 --readFilesIn QC/${sample}_1.fq.gz QC/${sample}_2.fq.gz \
 --outFilterMultimapNmax 20 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.1 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --outFilterType BySJout \
 --outFilterScoreMinOverLread 0.33 \
 --outFilterMatchNminOverLread 0.33 \
 --limitSjdbInsertNsj 1200000 \
 --outSAMstrandField intronMotif \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode GeneCounts \
 --outFileNamePrefix mapping/${sample}. \
 --outSAMattrIHstart 0 # for stringtie

elif [[ -e raw/${sample}.fastq.gz ]]; then

library="1"
$fastp \
 -i raw/${sample}.fastq.gz \
 -o QC/${sample}.fq.gz \
 -w $thread 2> log/$sample.log 

$STAR \
 --twopassMode Basic \
 --runThreadN $thread \
 --genomeDir $STAR_index \
 --readFilesCommand zcat \
 --readFilesIn QC/${sample}.fq.gz \
 --outFilterMultimapNmax 20 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.1 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --outFilterType BySJout \
 --outFilterScoreMinOverLread 0.33 \
 --outFilterMatchNminOverLread 0.33 \
 --limitSjdbInsertNsj 1200000 \
 --outSAMstrandField intronMotif \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode GeneCounts \
 --outFileNamePrefix mapping/${sample}. \
 --outSAMattrIHstart 0 # for stringtie

else

echo 'Error! There is no raw data'

fi

# coverage
$samtools index -@ $thread mapping/${sample}.Aligned.sortedByCoord.out.bam

${bamCoverage} \
 -b mapping/${sample}.Aligned.sortedByCoord.out.bam \
 -o mapping/${sample}.bw \
 -p ${thread} \
 --binSize 50 \
 --normalizeUsing RPGC \
 --effectiveGenomeSize $genomesize

# strand
tail -n+5 mapping/${sample}.ReadsPerGene.out.tab | awk 'BEGIN{FS="\t"}{sum3+=$3;sum4+=$4} END {total=sum3+sum4; print "-------"; print "--fr ratio: " sum3/total; print "--rf ratio: " sum4/total; print "------"}'

strandedness=$(tail -n+5 mapping/${sample}.ReadsPerGene.out.tab | awk 'BEGIN{FS="\t"}{sum3 += $3;sum4 += $4;} END {total = sum3 + sum4; ratio3 = sum3 / total; ratio4 = sum4 / total; if (ratio3 > 0.8) {strandedness = "--fr"} else if (ratio4 > 0.8) {strandedness = "--rf"} else {strandedness = ""} print strandedness}')

# stringtie assembly
$stringtie \
 mapping/${sample}.Aligned.sortedByCoord.out.bam \
 -G $GTF -p $thread -e ${strandedness} \
 -A quantify/${sample}.txt \
 -o quantify/${sample}.exp.gtf

# quantify
$mergeTPM $(ls quantify/*txt|tr "\n" " ")

# Splicing Quantify
mkdir splicing
/public/home/liunangroup/liangyan/software/McSplicer-master/bin/exonRefine \
 --prefix splicing/$sample quantify/$sample.exp.gtf

if [[ "$library" == "1" ]];then
/public/home/liunangroup/liangyan/software/McSplicer-master/bin/sigcount \
 --strand-specific \
 mapping/${sample}.Aligned.sortedByCoord.out.bam \
 splicing/$sample.gtf \
 splicing/$sample

else if [[ "$library" == "2" ]];then
/public/home/liunangroup/liangyan/software/McSplicer-master/bin/sigcount_pe \
 --strand-specific \
 mapping/${sample}.Aligned.sortedByCoord.out.bam \
 splicing/$sample.gtf \
 splicing/$sample

fi

read_len=$(samtools view $sample.Aligned.sortedByCoord.out.bam | head -n 10001 | cut -f10 | awk '{print length}'|sort|sed -n '5000p')
python /public/home/liunangroup/liangyan/software/McSplicer-master/python_scripts/McSplicer.py \
 -gtf splicing/$sample.gtf \
 --count_file splicing/$sample.cnt \
 --out_dir splicing \
 --read_len $read_len \
 --prefix $sample \
 --anno y
