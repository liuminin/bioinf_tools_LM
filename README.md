# bioinf_tools_LM

This is a small bioinformatics project in Python.  
Includes tools for DNA/RNA strings, FASTQ filtering and file processors (FASTA/BLAST).


## Project Structure

bioinf_tools_LM/
├─ bioinf_tools_LM.py 
├─ bio_files_processor.py 
├─ modules/
│ ├─ init.py 
│ ├─ fastq.py 
│ └─ run_dna_rna_tools.py 
└─ README.md 
---

## How to Run

Run the main script from the terminal:

python bioinf_tools_LM.py


The main file imports:
`run_dna_rna_tools` from `modules/run_dna_rna_tools.py`
`filter_fastq` from `modules/fastq.py`

---

## Main Functions

### run_dna_rna_tools
This file has basic DNA/RNA functions (e.g. reverse, transcribe, translate).

### filter_fastq
This function filters reads by  GC content, Length and Average quality (Phred33 scale)

### convert_multiline_fasta_to_oneline
Reads a multi-line FASTA and writes a new file where each sequence is on one line

### parse_blast_output
Parses text BLAST output and collects (per query) the FIRST hit from the section


## Example:

EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')


    }  
from bio_files_processor import convert_multiline_fasta_to_oneline
convert_multiline_fasta_to_oneline("example_data/example_multiline_fasta.fasta")
# -> creates "example_multiline_fasta_oneline.fasta"

from bio_files_processor import parse_blast_output
parse_blast_output("example_data/example_blast_results.txt", "blast_hits.txt")
# -> writes unique first-hit descriptions into ./blast_hits.txt

## Author:
Liubov Minina
