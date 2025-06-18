import streamlit as st
from Bio import SeqIO, Seq
from io import StringIO

def analyze_sequence(seq, is_dna=True):
    s = Seq.Seq(seq)
    result = {}
    result["Length"] = len(seq)
    result["GC Content (%)"] = round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
    
    base_counts = {}
    for base in ("ATGC" if is_dna else "AUGC"):
        base_counts[base] = seq.upper().count(base)
    result["Base Composition"] = base_counts

    if is_dna:
        result["Reverse Complement"] = str(s.reverse_complement())
        result["Transcription"] = str(s.transcribe())
    else:
        result["Reverse Complement"] = "N/A"
        result["Transcription"] = "Already RNA"

    result["Translation (to stop)"] = str(s.translate(to_stop=True))
    result["Full Translation"] = str(s.translate())
    result["Longest ORF"] = find_longest_orf(seq)
    return result

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

# --- STREAMLIT UI ---

st.title("🧬 DNA/RNA Sequence Analyzer")

uploaded_file = st.file_uploader("Upload a FASTA or FASTQ file", type=["fasta", "fa", "fastq"])
is_dna = st.radio("Is the input DNA?", ["Yes", "No"]) == "Yes"

if uploaded_file:
    file_format = "fastq" if uploaded_file.name.lower().endswith("fastq") else "fasta"
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    records = list(SeqIO.parse(stringio, file_format))

    for record in records:
        st.subheader(f"🧾 {record.id}")
        result = analyze_sequence(str(record.seq), is_dna)

        for k, v in result.items():
            if isinstance(v, dict):
                st.write(k)
                st.json(v)
            else:
                st.write(f"**{k}:**", v)
