thread="8"

for i in $(cat zzsample.txt);do
sbatch -J job -p cu -c $thread -o %j.out -e %j.err --wrap="
sh /public/home/liunangroup/liangyan/pipeline/RNAseq/runRNAseq_human.sh $i $thread
"
done
