#!/bin/bash

#SBATCH --job-name embed
#SBATCH --output logs/embed_%j.log
#SBATCH --error logs/embed_%j.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 35-00:00:00
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 500G

source /nfsd/bcb/bcbg/rossigno/anaconda3/bin/activate
conda activate dvf

cd ..

echo "Starting encoding virus sequences..."
head -n 100 ~/bio/DeepVirFinder/datasets/training/virus-training.fa > virus.fa
srun python embed.py -k 3 --file virus.fa --type 1

echo "Starting encoding host sequences..."
head -n 100 ~/bio/DeepVirFinder/datasets/training/host-training.fa > host.fa
srun python embed.py -k 3 --file host.fa --type 0

header=""
for i in $(seq 1 $((4**3))); do
	header+="kmer$i,"
done
header+="type"
echo $header > data.csv

echo "Joining files..."
cat virus.fa.csv host.fa.csv >> data.csv
rm virus.fa{,.csv} host.fa{,.csv}
echo "Done."
