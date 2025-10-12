allowed = set("ATUGCatugc")


def is_nucleic_acid(seq: str) -> bool:
    """Return TRUE if the string contains only A,T, U, G, C (in any case)
    and doesn't contain simultaneously T and U
    """
    
    if not seq:
        return False  # only allowed letters

    for ch in seq:
        if ch not in allowed:
            return False

    # exclude mix T and U
    has_t = ("T" in seq) or ("t" in seq)
    has_u = ("U" in seq) or ("u" in seq)
    if has_t and has_u:
        return False
    return True


def transcribe(seq: str) -> str:
    """Transcription DNA in RNA: replace T/t with U/u"""
    seq = seq.replace("T", "U")
    seq = seq.replace("t", "u")
    return seq


def reverse(seq: str) -> str:
    """Reverse string"""
    return seq[::-1]


def is_rna_pairs(seq: str) -> bool:
    """If U/u exists and T/t doesn't, use RNA pairing rules"""
    has_t = ("T" in seq) or ("t" in seq)
    has_u = ("U" in seq) or ("u" in seq)
    return (not has_t) and has_u


def complement(seq: str) -> str:
    """Return the complementary sequence.
    If U is present, use RNA pairs (A & U, C & G),
    otherwise use DNA pairs (A & T, C & G).
    """
    use_rna = is_rna_pairs(seq)

   
    dna_map = {
        "A": "T", "a": "t",
        "T": "A", "t": "a",
        "G": "C", "g": "c",
        "C": "G", "c": "g",
    }
    rna_map = {
        "A": "U", "a": "u",
        "U": "A", "u": "a",
        "G": "C", "g": "c",
        "C": "G", "c": "g",
    }
    cmap = rna_map if use_rna else dna_map
    return "".join(cmap.get(ch, ch) for ch in seq)


def reverse_complement(seq: str) -> str:
    return reverse(complement(seq))


def run_dna_rna_tools(*args):
    """Accepts any number of sequences and the 
    last argument is the procedure name.
    Allowed names: "is_nucleic_acid", "transcribe", "reverse",
    "complement", "reverse_complement"
    Returns a single value for one input or a list for multiple inputs.

    Examples:
        run_dna_rna_tools('TTUU', 'is_nucleic_acid') -> False
        run_dna_rna_tools('ATG', 'reverse') -> 'GTA'
        run_dna_rna_tools('ATG', 'aT', 'reverse') -> ['GTA', 'Ta']
    """
   
    if len(args) < 2:
        raise ValueError(
            "Provide at least one sequence and a name for the procedure"
        )

    *seqs, procedure = args

    funcs = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    if procedure not in funcs:
        raise ValueError(f"Unknown procedure: {procedure}")

    fn = funcs[procedure]
    results = [fn(s) for s in seqs]
    return results[0] if len(results) == 1 else results
