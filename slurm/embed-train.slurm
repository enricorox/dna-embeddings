#!/bin/bash

#SBATCH --job-name embed-train
#SBATCH --output logs/embed_%j.log
#SBATCH --error logs/embed_%j.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 35-00:00:00
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 500G
phase="training"
num_lines=1000

source /nfsd/bcb/bcbg/rossigno/anaconda3/bin/activate
conda activate dvf

cd ..

virus_file="temp-virus-$phase.fa"
host_file="temp-host-$phase.fa"
csv_file="$phase-data.csv"

echo "Starting encoding virus sequences..."
head -n $num_lines ~/bio/DeepVirFinder/datasets/"$phase"/virus-training.fa > "$virus_file"
srun python embed.py -k 3 --file "$virus_file" --label 1

echo "Starting encoding host sequences..."
head -n $num_lines ~/bio/DeepVirFinder/datasets/"$phase"/host-training.fa > "$host_file"
srun python embed.py -k 3 --file "$host_file" --label 0

header=""
for i in $(seq 1 $((4**3))); do
	header+="kmer$i,"
done
header+="label"
echo $header > "$csv_file"

echo "Joining files..."
cat "$virus_file".csv "$host_file".csv >> $csv_file
rm "$virus_file"{,.csv} "$host_file"{,.csv}
echo "Done."
