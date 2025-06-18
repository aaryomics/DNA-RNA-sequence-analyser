from Bio import SeqIO
from Bio.Seq import Seq

def analyze_sequence(seq, is_dna=True):
    s = Seq(seq)
    print("Length:", len(seq))
    print("GC Content: {:.2f}%".format((seq.count("G") + seq.count("C")) / len(seq) * 100))
    
    print("Base Composition:")
    for base in "ATGC" if is_dna else "AUGC":
        print(f"  {base}: {seq.upper().count(base)}")
    
    if is_dna:
        print("Reverse Complement:", str(s.reverse_complement()))
        print("Transcription (DNA → RNA):", s.transcribe())
    else:
        print("Transcription: Already RNA")

    print("Translation (to protein):", s.translate(to_stop=True))
    print("Longest ORF:", find_longest_orf(seq))

def find_longest_orf(seq):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for i in range(len(seq)):
        if seq[i:i+3] == start_codon:
            for j in range(i + 3, len(seq), 3):
                codon = seq[j:j+3]
                if codon in stop_codons:
                    orfs.append(seq[i:j+3])
                    break
    return max(orfs, key=len, default="None found")

def main():
    print("DNA/RNA Sequence Analyzer")
    file_path = input("Enter FASTA or FASTQ file path: ")
    file_format = "fastq" if file_path.lower().endswith(".fastq") else "fasta"
    is_dna = input("Are these DNA sequences? (y/n): ").lower() == 'y'

    records = list(SeqIO.parse(file_path, file_format))

    for record in records:
        print("\n===============================")
        print(f"ID: {record.id}")
        analyze_sequence(str(record.seq), is_dna)

if __name__ == "__main__":
    main()
