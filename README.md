# DNA/RNA Sequence Analyser
This is a simple web app for analysing DNA or RNA sequences. It supports FASTA and FASTQ file formats and provides useful sequence statistics such as GC content, reverse complement, transcription, translation, and longest ORF detection. Built using Python, Biopython, and Streamlit.

## Live Demo
You can try the app here:  
https://dna-rna-sequence-analyser-7ua8nxjjnnb2hrsh4jzqt8.streamlit.app/
![image](https://github.com/user-attachments/assets/af2be196-ed78-4d6e-ad10-4782c8cb14c6)


## Features
- Upload DNA or RNA sequences in FASTA or FASTQ format
- Automatic detection of file format
- GC content calculation
- Base composition (A/T/G/C or A/U/G/C)
- Reverse complement and transcription
- Translation to protein (until first stop codon and full sequence)
- Longest ORF (Open Reading Frame) finder on forward strand

## How to Run Locally
1. Clone the repository:
   git clone https://github.com/aaryomics/dna-rna-sequence-analyser.git
   cd dna-rna-sequence-analyzer

2. Install dependencies:
   pip install -r requirements.txt

3. Start the Streamlit app:
   streamlit run streamlit_app.py
   
   The app will open in your browser.

## Example Input
Example (FASTA format):
<pre> ``` >Example_Sequence ATGCGTAGCTAGCTACGATCGATCGTAGCTAGCTAGCTACGATCG ``` </pre>

## Technologies Used
- Python
- Biopython
- Streamlit

## Author
Created by Aarya Kumar  
GitHub: https://github.com/aaryomics 
