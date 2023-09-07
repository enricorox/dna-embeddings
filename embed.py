import argparse

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DNA sequence embedder')
    parser.add_argument("--file", type=str, help="fasta file to embed")
    parser.add_argument("-k", type=int, default=31, help="k-mer size")
    parser.add_argument("-N", type=int, default=500, help="sequence length")
    parser.add_argument("--emb", type=str, default="spike2vec", help="embedding name")

    args = parser.parse_args()

