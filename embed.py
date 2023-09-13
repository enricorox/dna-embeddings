import argparse
import csv

from Bio import SeqIO

alphabet = {'A': 0, 'C': 1, 'T': 2, 'G': 3, 'U': 2}
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
    seq = seq.upper()

    for j in range(len(seq) - k + 1):
        codon = seq[j:j+k]
        freqs[kmer_id(codon)] += 1
    return freqs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DNA sequence embedder')
    parser.add_argument("--file", type=str, help="fasta file to embed")
    parser.add_argument("-k", type=int, default=3, help="k-mer size")
    parser.add_argument("-L", type=int, default=500, help="sequence length")
    parser.add_argument("--emb", type=str, default="sec2vec", help="embedding name")
    parser.add_argument("--label", type=int, help="an integer identifying the class")

    args = parser.parse_args()

    print(f"Start embedding with {args.emb}...")

    tot_count = 0
    uncertains = 0
    unknowns = 0
    if args.emb == "sec2vec-overlapping":
        with open(args.file) as fasta_file:
            with open(args.file + ".csv", 'w') as csv_file:
                csv_writer = csv.writer(csv_file)
                # header = [f"kmer{i}," for i in range(4 ** args.k)] + ["label"]
                # csv_writer.writerow(header)
                for record in SeqIO.parse(fasta_file, "fasta"):
                    for i in range(len(record.seq) - args.L + 1):
                        lmer = record.seq[i:i+args.L]

                        # skip unknown nucleotides
                        if 'N' in lmer:
                            uncertains += 1
                            continue
                        if not ('A' in lmer or 'C' in lmer or 'G' in lmer or 'T' in lmer or 'U' in lmer):
                            unknowns += 1
                            continue

                        freq = compute_frequencies(lmer, args.k)
                        csv_writer.writerow(freq + [args.label])
                        tot_count += 1
    elif args.emb == "sec2vec":
        with open(args.file) as fasta_file:
            with open(args.file + ".csv", 'w') as csv_file:
                csv_writer = csv.writer(csv_file)
                # header = [f"kmer{i}," for i in range(4 ** args.k)] + ["label"]
                # csv_writer.writerow(header)
                for record in SeqIO.parse(fasta_file, "fasta"):
                    for l in range(len(record.seq)//args.L):
                        i = l * args.L
                        lmer = record.seq[i:i+args.L]

                        # skip unknown nucleotides
                        if 'N' in lmer:
                            uncertains += 1
                            continue
                        if not ('A' in lmer or 'C' in lmer or 'G' in lmer or 'T' in lmer or 'U' in lmer):
                            unknowns += 1
                            continue

                        freq = compute_frequencies(lmer, args.k)
                        csv_writer.writerow(freq + [args.label])
                        tot_count += 1

    print(f"Embedded sequences with label {args.label}:")
    print(f"\tnumber of vectors: {tot_count}")
    print(f"\tnumber of uncertain nucleotides: {uncertains}")
    print(f"\tnumber of unknown nucleotides: {unknowns}")
    print("Done!")
