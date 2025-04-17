# ------------------------

thread="8"

for i in $(cat zzsample.txt);do
sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
sh /public/home/liunangroup/liangyan/pipeline/ATACseq/ATACseq_human.sh $i $thread
"
done

# ------------------------

thread="8"

sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
for i in $(cat zzsample.txt);do
sh /public/home/liunangroup/liangyan/pipeline/ATACseq/ATACseq_human.sh $i $thread
done
"
