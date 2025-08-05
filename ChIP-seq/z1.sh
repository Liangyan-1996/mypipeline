thread="8"

for sample in $(cat zzsample.txt);do
sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
sh /public/home/liunangroup/liangyan/pipeline/mypipeline/ChIP-seq/ChIPseq_mouse.sh $sample $thread
"
done
