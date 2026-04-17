
from abc import ABC, abstractmethod
from typing import Tuple
import os
import logging
import argparse

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    Requirements:
    - len()
    - indexing and slicing
    - alphabet validation
    """
    def __init__(self, sequence: str):
        self.sequence = str(sequence).upper()

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, item):
        # For indexing return a single character (str)
        # For slicing return an object of the class
        if isinstance(item, slice):
            return self.__class__(self.sequence[item])
        return self.sequence[item]

    def __str__(self) -> str:
        return f"{self.__class__.__name__}('{self.sequence}')"

    def __repr__(self) -> str:
        return str(self)

    @property
    @abstractmethod
    def alphabet(self) -> set:
        """Allowed symbols for the sequence"""
        raise NotImplementedError

    def check_alphabet(self) -> bool:
        """Return True if all characters belong to the allowed alphabet"""
        for ch in self.sequence:
            if ch not in self.alphabet:
                return False
        return True


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for DNA/RNA sequences.
    Implements complement/reverse/reverse_complement

    """
    @property
    def alphabet(self) -> set:
        # Common nucleic alphabet 
        return {"A", "C", "G", "T", "U"}

    @property
    @abstractmethod
    def _complement_map(self) -> dict:
        """Mapping for complement"""
        raise NotImplementedError

    def complement(self):
        # No if DNA/RNA inside this method
        # Uses subclass-provided _complement_ma
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("Use DNASequence or RNASequence")
        if not self.check_alphabet():
            raise ValueError("Sequence contains invalid symbols for nucleic acids")
        comp_chars = []
        for ch in self.sequence:
            if ch not in self._complement_map:
                raise ValueError(f"Cannot complement symbol {ch}")
            comp_chars.append(self._complement_map[ch])
        return self.__class__("".join(comp_chars))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    @property
    def alphabet(self) -> set:
        return {"A", "C", "G", "T"}

    @property
    def _complement_map(self) -> dict:
        return {"A": "T", "T": "A", "C": "G", "G": "C"}

    def transcribe(self):
        # DNA --> RNA
        if not self.check_alphabet():
            raise ValueError("Invalid DNA alphabet")
        return RNASequence(self.sequence.replace("T", "U"))
    
   


class RNASequence(NucleicAcidSequence):
    @property
    def alphabet(self) -> set:
        return {"A", "C", "G", "U"}

    @property
    def _complement_map(self) -> dict:
        return {"A": "U", "U": "A", "C": "G", "G": "C"}
   
   

class AminoAcidSequence(BiologicalSequence):
    @property
    def alphabet(self) -> set:
        return set("ACDEFGHIKLMNPQRSTVWY")

    def molecular_weight_like(self) -> int:
        """ 
        Returns sequence length 
        """
        if not self.check_alphabet():
            raise ValueError("Invalid amino acid alphabet")
        return len(self)


def normalize_bounds(bounds, default_low, default_high):
    """Normally bounds are a pair (low, high). If bounds is a single number, it's a highest bound"""
    # turn single number into a pair
    if isinstance(bounds, (int, float)):
        low, high = float(default_low), float(bounds)
    else:
        low, high = float(bounds[0]), float(bounds[1])
    # swap if reversed
    if high < low:
        low, high = high, low
    return low, high

def setup_logger(log_file="fastq_filter.log"):
    logger = logging.getLogger("fastq_filter")
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        file_handler = logging.FileHandler(log_file)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Tuple[float, float] | float = (0, 100),
    length_bounds: Tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
    log_file: str = "fastq_filter.log",
) -> int:
    """
    Filter a FASTQ file using Biopython:
    - length bounds
    - GC bounds
    - average Phred quality threshold

    Writes filtered reads into filtered/output_fastq
    Returns number of reads written
    """
    logger = setup_logger(log_file)
    low_gc, high_gc = normalize_bounds(gc_bounds, 0, 100)
    low_len, high_len = normalize_bounds(length_bounds, 0, 2**32)
    min_q = float(quality_threshold)

    os.makedirs("filtered", exist_ok=True)
    out_path = os.path.join("filtered", output_fastq)

    if os.path.exists(out_path):
        logger.error(f"Output file already exists {out_path}")
        raise FileExistsError(f"File already exists {out_path}. Choose other name")

    passed_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        seq_len = len(record.seq)
        if not (low_len <= seq_len <= high_len):
            continue

        # GC fraction returns 0..1, convert to percent 0..100
        gc_percent = gc_fraction(record.seq) * 100.0
        if not (low_gc <= gc_percent <= high_gc):
            continue

        phred_scores = record.letter_annotations.get("phred_quality", [])
        if not phred_scores:
            continue
        avg_q = sum(phred_scores) / len(phred_scores)
        if avg_q < min_q:
            continue

        passed_records.append(record)

    # write output
    written = SeqIO.write(passed_records, out_path, "fastq")
    logger.info(f"Filtering finishd successfully. Written reads {written}")
    return int(written)

from modules.fastq import filter_fastq_file as filter_fastq
from modules.run_dna_rna_tools import run_dna_rna_tools

def parse_args():
    parser = argparse.ArgumentParser(description="FASTQ filter tool")

    parser.add_argument("input_fastq", help="Path to inpit FASTQ file")
    parser.add_argument("output_fastq", help="Name of output FASTQ file")

    parser.add_argument("--gc_min", type=float, default=0)
    parser.add_argument("--gc_max", type=float, default=100)

    parser.add_argument("--len_min", type=int, default=0)
    parser.add_argument("--len_max", type=int, default=2**32)

    parser.add_argument("--quality", type=float, default=0)
    parser.add_argument("--log_file", default="fastq_filter.log")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    written_reads = filter_fastq(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=(args.gc_min, args.gc_max),
        length_bounds=(args.len_min, args.len_max),
        quality_threshold=args.quality,
        log_file=args.log_file,
    )

    print(f"Written reads {written_reads}")

