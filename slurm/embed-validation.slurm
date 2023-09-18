#!/bin/bash

#SBATCH --job-name embed-validation
#SBATCH --output logs/embed_validation_%j.log
#SBATCH --error logs/embed_validation_%j.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 35-00:00:00
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 500G
phase="validation"
num_lines=100000
whole=true

source /nfsd/bcb/bcbg/rossigno/anaconda3/bin/activate
conda activate dvf

cd ..

virus_file="temp-virus-$phase.fa"
host_file="temp-host-$phase.fa"
csv_file="$phase-data.csv"

echo "Encoding virus sequences for $phase..."
if [ "$whole" == "true" ]; then
  virus_file=~/bio/DeepVirFinder/datasets/"$phase"/virus-"$phase".fa
else
  head -n $num_lines ~/bio/DeepVirFinder/datasets/"$phase"/virus-"$phase".fa > "$virus_file"
fi
srun python embed.py -k 3 --file "$virus_file" --label 1

echo "Encoding host sequences for $phase..."
if [ "$num_lines" == "-1" ]; then
  host_file=~/bio/DeepVirFinder/datasets/"$phase"/host-"$phase".fa
else
  head -n $num_lines ~/bio/DeepVirFinder/datasets/"$phase"/host-"$phase".fa > "$host_file"
fi
srun python embed.py -k 3 --file "$host_file" --label 0

header=""
for i in $(seq 1 $((4**3))); do
	header+="kmer$i,"
done
header+="label"
echo $header > "$csv_file"

echo "Joining files..."
cat "$virus_file".csv "$host_file".csv >> $csv_file
rm "$virus_file".csv "$host_file".csv
[ "$whole" == "true" ] || rm "$virus_file" "$host_file"
echo "Done."
