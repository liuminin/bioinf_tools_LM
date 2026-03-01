from abc import ABC, abstractmethod
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
    Implements complement/reverse/reverse_complement.

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
                raise ValueError(f"Cannot complement symbol: {ch}")
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
        if not NucleicAcidSequence.check_alphabet(self):
            raise ValueError("Invalid DNA alphabet")
        return RNASequence(self.sequence.replace("T", "U"))
    


    # DNASequence shouldn't expose these methods publicly
    def complement(self):
        raise AttributeError("Use NucleicAcidSequence.complement()")

    def reverse(self):
        raise AttributeError("Use NucleicAcidSequence.reverse()")

    def reverse_complement(self):
        raise AttributeError("Use NucleicAcidSequence.reverse_complement()")

   


class RNASequence(NucleicAcidSequence):
    @property
    def alphabet(self) -> set:
        return {"A", "C", "G", "U"}

    @property
    def _complement_map(self) -> dict:
        return {"A": "U", "U": "A", "C": "G", "G": "C"}
    
    #  RNASequence shouldn't expose these methods publicly
    def complement(self):
        raise AttributeError("Use NucleicAcidSequence.complement()")

    def reverse(self):
        raise AttributeError("Use NucleicAcidSequence.reverse()")

    def reverse_complement(self):
        raise AttributeError("Use NucleicAcidSequence.reverse_complement()")

   
   

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






from typing import Tuple
import os

def normalize_bounds(bounds, default_low, default_high):
    """Normally bounds are a pair (low, high). If bounds is a single number, it is a highest bound"""
    # turn single number into a pair
    if isinstance(bounds, (int, float)):
        low, high = float(default_low), float(bounds)
    else:
        low, high = float(bounds[0]), float(bounds[1])
    # swap if reversed
    if high < low:
        low, high = high, low
    return low, high



from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Tuple[float, float] | float = (0, 100),
    length_bounds: Tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> int:
    """
    Filter a FASTQ file using Biopython:
    - length bounds
    - GC bounds
    - average Phred quality threshold

    Writes filtered reads into filtered/output_fastq
    Returns number of reads written
    """
    low_gc, high_gc = normalize_bounds(gc_bounds, 0, 100)
    low_len, high_len = normalize_bounds(length_bounds, 0, 2**32)
    min_q = float(quality_threshold)

    os.makedirs("filtered", exist_ok=True)
    out_path = os.path.join("filtered", output_fastq)

    if os.path.exists(out_path):
        raise FileExistsError(f"File already exists: {out_path}. Choose another name")

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
    return int(written)


