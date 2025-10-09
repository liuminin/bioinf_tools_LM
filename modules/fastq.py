
def normalize_bounds(bounds, default_low, default_high):
    """
    Make sure bounds are a pair (low, high).
    If bounds is a single number, treat it as a highest bound (e.g. (default_low, 80)).

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

def gc_percent(seq: str) -> float:
    """
    Calculate GC content in percent (0-100%).
    Only 'G' and 'C' are used.
    """
    if not seq:
        return 0.0
    s = seq.upper()
    gc = 0
    for ch in s:
        if ch == 'G' or ch == 'C':
            gc += 1
    return 100 * gc / len(s)


def average_quality_phred33(qual: str) -> float:
    """
    Compute average Phred33 quality for a read.
    Given formula: score = ord(char) - 33
    Then take the mean of all scores in the quality string.
    """
    if not qual:
        return 0.0
    total = 0
    for ch in qual:
        total += (ord(ch)-33)
    return total / len(qual)


def filter_fastq(
    seqs,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
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
    # Turn possibly single numbers inputs into (low, high) pairs
    low_gc, high_gc = normalize_bounds(gc_bounds, 0, 100)
    low_len, high_len = normalize_bounds(length_bounds, 0, 2**32)
    q_min = float(quality_threshold)

    result = {}

    # Go through each read
    for read_id, pair in seqs.items():
        # Pair must be a pair (sequence, quality)
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            # If the data is malformed, skip it
            continue

        seq, qual = pair

        # 1. Length filter 
        L = len(seq)
        if not (low_len <= L <= high_len):
            continue

        # 2. GC filter 
        gc = gc_percent(seq)
        if not (low_gc <= gc <= high_gc):
            continue

        # 3. Average quality filter (>= threshold)
        avg_q = average_quality_phred33(qual)
        if avg_q < q_min:
            continue

        # If all checks passed, then keep this read
        result[read_id] = (seq, qual)

    return result
