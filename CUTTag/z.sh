thread="16"

for i in `cat zzsample.txt`;do
sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
sh runCUTTag.sh $i $thread 
"
done
