
# bioinf_tools_LM (HW16)

Small bioinformatics utilities in Python.

This homework implementing the previous semester’s functional code into OOP classes for biological sequences
and re-implements FASTQ filtering using Biopython.

---

## Project Structure

bioinf_tools_LM/  
├─ bioinf_tools_LM.py  
└─ README.md  

Note: The `modules/` subpackage from the previous semester is no longer used in this project.

---

## Requirements

- Python 3.x
- Biopython
---

## How to Run

This repository provides classes inside bioinf_tools_LM.py.
You can import and use them from Python.

python -c "from bioinf_tools_LM import DNASequence; print(DNASequence('ATGC').transcribe())"
---
## Task 1: Abstract Sequences (OOP)

Implemented classes:
 ### BiologicalSequence (abstract)
- supports len(seq)
- indexing and slicing (seq[i], seq[a:b])
 ### NucleicAcidSequence (abstract base for nucleic acids)
- complement(), reverse(), reverse_complement()
- complement logic is polymorphic (no if DNA/RNA inside)
 ### DNASequence
- transcribe() -> returns an RNASequence
### RNASequence
---
## Task 2: FASTQ Filter (Biopython)
Function:

### filter_fastq
- filters FASTQ reads by length bounds, GC% bounds, average Phred quality threshold
- implemented using Biopython (SeqIO, SeqUtils, SeqRecord)
---
## Author:
Liubov Minina