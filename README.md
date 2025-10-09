# bioinf_tools_LM

This is a small bioinformatics project in Python.  
It has tools for DNA/RNA sequences and for filtering FASTQ reads.


## Project Structure

├── bioinf_tools_LM.py
├── modules
│   ├── __init__.py
│   ├── fastq.py
│   └── run_dna_rna_tools.py
└── README.md

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

## Example:

EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')

    }  

## Author:
Liubov Minina