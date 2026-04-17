import os
import sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from bioinf_tools_LM import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    normalize_bounds,
    filter_fastq,
)


def test_len_dna():
    seq = DNASequence("ATGC")
    assert len(seq) == 4


def test_index_dna():
    seq = DNASequence("ATGC")
    assert seq[1] == "T"


def test_slice_dna():
    seq = DNASequence("ATGC")
    result = seq[1:3]
    assert isinstance(result, DNASequence)
    assert result.sequence == "TG"


def test_transcribe_dna():
    seq = DNASequence("ATGC")
    rna = seq.transcribe()
    assert isinstance(rna, RNASequence)
    assert rna.sequence == "AUGC"


def test_reverse_complement_dna():
    seq = DNASequence("ATGC")
    result = seq.reverse_complement()
    assert isinstance(result, DNASequence)
    assert result.sequence == "GCAT"


def test_amino_acid_invalid_alphabet():
    seq = AminoAcidSequence("ACDZ")
    assert seq.check_alphabet() is False


def test_normalize_bounds_single_value():
    low, high = normalize_bounds(50, 0, 100)
    assert low == 0
    assert high == 50


def test_filter_fastq_file_exists(tmp_path):
    input_fastq = tmp_path / "input.fastq"
    output_dir = tmp_path / "filtered"
    output_dir.mkdir()
    output_fastq = output_dir / "result.fastq"

    input_fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
    )

    output_fastq.write_text("already exists")

    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with pytest.raises(FileExistsError):
            filter_fastq(
                input_fastq=str(input_fastq),
                output_fastq="result.fastq",
            )
    finally:
        os.chdir(old_cwd)