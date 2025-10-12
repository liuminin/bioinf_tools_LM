"""
bio_files_processor.py

Bio tilities for working with biological text files:
FASTA and BLAST output
"""

from typing import Optional, List


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: Optional[str] = None) -> str:
    """
    Convert a multiline FASTA file into a one line per sequence FASTA file.

    FASTA files often split long sequences into multiple lines.
    convert_multiline_fasta_to_oneline function joins all sequence lines belonging to the same header.

    Parameters:

    input_fasta : str
        Path to the input FASTA file (can have multi-line sequences)
    output_fasta : str, optional
        Path for the resulting one-line FASTA file.
        If not provided, '_oneline.fasta' will be created next to the input file.

    Returns:
   
    str
        Path to the written output file.

    """
    # Determine output file path
    if output_fasta is None:
        if input_fasta.endswith(".fasta"):
            output_fasta = input_fasta.replace(".fasta", "_oneline.fasta")
        else:
            output_fasta = input_fasta + "_oneline.fasta"

    # Prepare variables for parsing
    sequences: List[str] = []  # lines for sequence
    header: Optional[str] = None

    with open(input_fasta, "r", encoding="utf-8") as inp, open(output_fasta, "w", encoding="utf-8") as out:
        for line in inp:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith(">"):
                # If we have a header, write its sequence before switching
                if header is not None:
                    out.write(f"{header}\n")
                    out.write("".join(sequences) + "\n")
                    sequences = []
                # Start new header
                header = line
            else:
                # Sequence line (append to the current one)
                sequences.append(line)

        # Write the last sequence (if file doesn’t end with '>')
        if header is not None:
            out.write(f"{header}\n")
            out.write("".join(sequences) + "\n")

    return output_fasta


def parse_blast_output(input_file: str, output_file: str) -> str:
    """
    Parse a plain-text BLAST output file and extract top hits for each query.

    For each QUERY block, find the first hit under the section:
        'Sequences producing significant alignments:'
    and extract its description.

    Save all unique descriptions into a new file, sorted alphabetically.

    Parameters:
    
    input_file : str
        Path to the BLAST text output file.
    output_file : str
        Path for saving the extracted protein/gene descriptions.

    Returns:
   
    str
        Path to the written output file.
    """
    hits: List[str] = []
    inside_query = False
    inside_hits = False

    with open(input_file, "r", encoding="utf-8") as inp:
        for line in inp:
            line = line.strip()

            # detect start of a query
            if line.startswith("Query="):
                inside_query = True
                inside_hits = False
                continue

            # detect section with hits
            if "Sequences producing significant alignments:" in line:
                inside_hits = True
                continue

            # if inside the hits section
            if inside_hits:
                # stop if it reaches an empty line or next section
                if line == "" or line.startswith(">"):
                    inside_hits = False
                    continue

                
                # take the full description (until the end of line)
                description = line.split()[0]
                hits.append(description)

    # Deduplicate and sort
    unique_hits = sorted(set(hits))

    with open(output_file, "w", encoding="utf-8") as out:
        for hit in unique_hits:
            out.write(hit + "\n")

    return output_file