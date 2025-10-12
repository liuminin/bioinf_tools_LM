from typing import Dict, Tuple
import os


def normalize_bounds(bounds, default_low, default_high):
    """
    Normally bounds are a pair (low, high).
    If bounds is a single number, it is a highest bound
    """
    # turn single number into a pair
    if isinstance(bounds, (int, float)):
        low, high = float(default_low), float(bounds)
    else:
        low, high = float(bounds[0]), float(bounds[1])

    # swap if reversed
    if high < low:
        low, high = high, low
    return low, high


def get_gc_percent(seq: str) -> float:
    """
    Calculate GC content in percent (0–100%).
    Only 'G' and 'C' are used
    """
    if not seq:
        return 0.0

    seq_lower = seq.lower()
    gc_count = seq_lower.count("g") + seq_lower.count("c")
    return (gc_count / len(seq_lower)) * 100.0


def get_average_quality_phred33(quality: str) -> float:
    """
    Compute average Phred33 quality for a read.

    Given formula:
        score = ord(char) - 33
    Then take the mean of all the scores in the quality string
    """
    if not quality:
        return 0.0

    total_score = 0
    for symbol in quality:
        total_score += (ord(symbol) - 33)
    return total_score / len(quality)


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Tuple[float, float] | float = (0, 100),
    length_bounds: Tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Keep only reads that pass all filters:

      1. GC % is within gc_bounds
      2. Length is within length_bounds
      3. Average Phred33 quality >= quality_threshold

    Input:
        seqs: dict like {"read_id": (sequence, quality_string)}
    Output:
        dict with the same structure, but only with passed reads
    """
    # Turn possibly single-number inputs into (low, high) pairs
    low_gc, high_gc = normalize_bounds(gc_bounds, 0, 100)
    low_length, high_length = normalize_bounds(length_bounds, 0, 2**32)
    min_quality = float(quality_threshold)

    result: Dict[str, Tuple[str, str]] = {}

    # Go through each read
    for read_id, pair in seqs.items():
        # Pair must be a pair (sequence, quality)
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            # If the data is malformed, then skip it
            continue

        sequence, quality = pair

        # Length filter
        seq_length = len(sequence)
        if not (low_length <= seq_length <= high_length):
            continue

        # GC filter
        gc_percent_value = get_gc_percent(sequence)
        if not (low_gc <= gc_percent_value <= high_gc):
            continue

        # Average quality filter (>= threshold)
        avg_quality = get_average_quality_phred33(quality)
        if avg_quality < min_quality:
            continue

        # If all checks passed, then keep this read
        result[read_id] = (sequence, quality)

    return result



# HW 5


def read_fastq(path: str) -> Dict[str, Tuple[str, str]]:
    """
    Read a FASTQ file and convert it into a dict format:
    {read_id: (sequence, quality)}

    Parameters:
   
    path : str
        Path to the input FASTQ file

    Returns:
    
    Dict[str, Tuple[str, str]]
        Mapping read_id -> (sequence, quality)
    """
    result: Dict[str, Tuple[str, str]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                break

            sequence = fh.readline()
            plus_line = fh.readline()
            quality = fh.readline()

            if not (header and sequence and plus_line and quality):
                # File ended unexpectedly, stop parsing
                break

            header = header.strip()
            sequence = sequence.strip()
            plus_line = plus_line.strip()
            quality = quality.strip()

            # Expect '@' at the beginning of header and '+' line
            if not header.startswith("@"):
                # Skip this record
                continue
            if not plus_line.startswith("+"):
                # Skip this record
                continue

            read_id = header  # keep full header as a key
            if len(sequence) != len(quality):
                # Malformed record: lengths must match
                continue

            result[read_id] = (sequence, quality)

    return result


def save_fastq(output_fastq: str, seqs: Dict[str, Tuple[str, str]]) -> str:
    """
    Save reads from the dict into a FASTQ file inside the filtered folder.

    To avoid accidental overwrite, if 'filtered/output_fastq' already exists, then
    raise FileExistsError

    Parameters:
    
    output_fastq : str
        File name without directory
    seqs : Dict[str, Tuple[str, str]]
        Mapping read_id -> (sequence, quality)

    Returns:
    str
        The full path to the written file.
    """
    os.makedirs("filtered", exist_ok=True)
    out_path = os.path.join("filtered", output_fastq)

    if os.path.exists(out_path):
        # Prevent accidental overwrite
        raise FileExistsError(
            f"File already exists: {out_path}. Choose another name"
        )

    with open(out_path, "w", encoding="utf-8") as out:
        for read_id, (sequence, quality) in seqs.items():
            out.write(f"{read_id}\n")
            out.write(f"{sequence}\n")
            out.write("+\n")
            out.write(f"{quality}\n")

    return out_path


def filter_fastq_file(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Tuple[float, float] | float = (0, 100),
    length_bounds: Tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> int:
    """
    Read a FASTQ file, filter it with the existing dict-based 'filter_fastq',
    and write filtered read into 'filtered/output_fastq'

    Parameters:
    
    input_fastq : str
        Path to the input FASTQ
    output_fastq : str
        File name for output (written into the 'filtered/' folder)
    gc_bounds : (low, high) or float
        Allowed GC bounds. If a single number is given, it's the upper bound.
    length_bounds : (low, high) or single int
        Allowed length bounds. If a single number is given, it's the upper bound
    quality_threshold : float
        Minimal average read quality Phred33 to keep.

    Returns:
    
    int
        Number of reads written to the output file.
    """
    reads = read_fastq(input_fastq)

    # filter using the earlier function
    filtered = filter_fastq(
        reads,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )

    save_fastq(output_fastq, filtered)

    return len(filtered)
