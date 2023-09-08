import argparse
import csv

from Bio import SeqIO

alphabet = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
def kmer_id(kmer):
    k_id = 0
    base = 1
    for c in kmer[::-1]:
        k_id += alphabet[c] * base
        base *= 4
    return k_id

def test_kmer_id():
    for c1 in alphabet:
        for c2 in alphabet:
            for c3 in alphabet:
                mer = c1+c2+c3
                print(f"{mer} --> {kmer_id(mer)}")

def compute_frequencies(seq, k):
    freqs = [0] * (4**k)
    for i in range(len(seq) - k + 1):
        codon = seq[i:i+k]
        # skip unknown nucleotides
        if 'N' not in codon:
            freqs[kmer_id(codon)] += 1
    return freqs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DNA sequence embedder')
    parser.add_argument("--file", type=str, help="fasta file to embed")
    parser.add_argument("-k", type=int, default=3, help="k-mer size")
    parser.add_argument("-L", type=int, default=500, help="sequence length")
    parser.add_argument("--emb", type=str, default="spike2vec", help="embedding name")
    parser.add_argument("--type", type=int, help="an integer identifying the class")

    args = parser.parse_args()

    with open(args.file) as fasta_file:
        with open(args.file + ".csv", 'w') as csv_file:
            csv_writer = csv.writer(csv_file)
            # header = [f"kmer{i}," for i in range(4 ** args.k)] + ["type"]
            # csv_writer.writerow(header)
            for record in SeqIO.parse(fasta_file, "fasta"):
                freq = compute_frequencies(record.seq, args.k)
                csv_writer.writerow(compute_frequencies(record.seq, args.k) + [args.classn])

