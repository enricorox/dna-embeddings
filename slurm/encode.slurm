#!/bin/bash

#SBATCH --job-name encode
#SBATCH --output encode.log
#SBATCH --error encode.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 35-00:00:00
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 500G

source /nfsd/bcb/bcbg/rossigno/anaconda3/bin/activate
conda activate dvf

echo "Starting encoding virus sequences..."
head -n 100 ../DeepVirFinder/datasets/training/virus-training.fa > virus.fa
srun python embed.py -k 3 --file virus.fa --classn 1

echo "Starting encoding host sequences..."
head -n 100 ../DeepVirFinder/datasets/training/host-training.fa > host.fa
srun python embed.py -k 3 --file host.fa --classn 0

echo "Joining files..."
header=""
for i in `seq 1 $((4**3))`; do
	header+=kmer$i,
done
echo $header type > data.csv
cat virus.fa.csv host.fa.csv >> data.csv
rm virus.fa{,.csv} host.fa{,.csv}
echo "Done."
