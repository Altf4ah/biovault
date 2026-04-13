"""
BioVault Encoder
Converts binary data → DNA base sequences with GC balancing and index headers.
"""

import os
import hashlib
from dataclasses import dataclass
from typing import List

# 2-bit encoding map
BASE_MAP = {"00": "A", "01": "T", "10": "C", "11": "G"}
BITS_MAP = {v: k for k, v in BASE_MAP.items()}

# GC content target range (synthesis sweet spot)
GC_MIN = 0.40
GC_MAX = 0.60

CHUNK_SIZE = 150  # bases per oligo (synthesis-friendly length)


@dataclass
class DNAFragment:
    index: int
    sequence: str
    checksum: str
    original_length: int  # bytes in this chunk

    def to_fasta(self) -> str:
        return f">fragment_{self.index:06d}|len={self.original_length}|chk={self.checksum}\n{self.sequence}\n"

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "sequence": self.sequence,
            "checksum": self.checksum,
            "original_length": self.original_length,
        }


def bytes_to_bits(data: bytes) -> str:
    return "".join(f"{byte:08b}" for byte in data)


def bits_to_dna(bits: str) -> str:
    # Pad to even length
    if len(bits) % 2:
        bits += "0"
    return "".join(BASE_MAP[bits[i : i + 2]] for i in range(0, len(bits), 2))


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def complement_base(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}[base]


def balance_gc(seq: str) -> str:
    """
    If GC content is out of range, flip bases strategically.
    This is a heuristic — real systems use constraint-aware synthesis.
    """
    gc = gc_content(seq)
    seq = list(seq)

    if gc < GC_MIN:
        # Too many A/T — flip some to G/C
        for i, base in enumerate(seq):
            if base in ("A", "T") and gc_content("".join(seq)) < GC_MIN:
                seq[i] = complement_base(base) if base in ("G", "C") else ("C" if base == "A" else "G")
    elif gc > GC_MAX:
        # Too many G/C — flip some to A/T
        for i, base in enumerate(seq):
            if base in ("G", "C") and gc_content("".join(seq)) > GC_MAX:
                seq[i] = "A" if base == "G" else "T"

    return "".join(seq)


def fragment_checksum(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()[:8]


def encode(data: bytes, chunk_bytes: int = 75) -> List[DNAFragment]:
    """
    Encode raw bytes into a list of DNA fragments.

    Args:
        data: Raw bytes to encode
        chunk_bytes: Bytes per fragment (75 bytes → 150 bases per oligo)

    Returns:
        List of DNAFragment objects ready for synthesis or storage
    """
    fragments = []
    chunks = [data[i : i + chunk_bytes] for i in range(0, len(data), chunk_bytes)]

    for idx, chunk in enumerate(chunks):
        bits = bytes_to_bits(chunk)
        seq = bits_to_dna(bits)
        chk = fragment_checksum(chunk)

        fragments.append(
            DNAFragment(
                index=idx,
                sequence=seq,
                checksum=chk,
                original_length=len(chunk),
            )
        )

    return fragments


def to_fasta(fragments: List[DNAFragment], filename: str = "biovault_output") -> str:
    """Export fragments as FASTA format (standard bioinformatics format)."""
    header = f"; BioVault encoded file: {filename}\n"
    header += f"; Total fragments: {len(fragments)}\n"
    header += f"; Total bases: {sum(len(f.sequence) for f in fragments)}\n\n"
    return header + "".join(f.to_fasta() for f in fragments)
