"""
BioVault Primer Design
Generates PCR primer sequences for selective retrieval of stored data.

In wet-lab DNA storage, you don't read ALL stored sequences every time.
Instead, you use primers — short (~20bp) flanking sequences — to selectively
amplify only the fragments you want, like a biological search query.

This module:
- Designs forward/reverse primers for each fragment or shard group
- Validates Tm (melting temperature), GC%, hairpin potential
- Outputs primer pairs ready for synthesis
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import hashlib
import math


# Constants
PRIMER_LENGTH = 20          # Standard PCR primer length
TM_MIN = 55.0               # Min melting temperature (°C)
TM_MAX = 65.0               # Max melting temperature (°C)
GC_MIN = 0.40
GC_MAX = 0.60
MAX_HOMOPOLYMER = 4         # Max consecutive same base (e.g., AAAA bad)


@dataclass
class Primer:
    sequence: str
    direction: str          # "forward" or "reverse"
    tm: float
    gc_content: float
    warnings: List[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return len(self.warnings) == 0

    def __str__(self):
        flag = "✓" if self.is_valid else "⚠"
        return (
            f"{flag} [{self.direction:>7}] 5'-{self.sequence}-3' "
            f"Tm={self.tm:.1f}°C GC={self.gc_content:.0%}"
            + (f" WARN: {'; '.join(self.warnings)}" if self.warnings else "")
        )


@dataclass
class PrimerPair:
    fragment_index: int
    forward: Primer
    reverse: Primer
    amplicon_length: int    # Expected PCR product size

    def to_dict(self) -> dict:
        return {
            "fragment_index": self.fragment_index,
            "forward_primer": self.forward.sequence,
            "reverse_primer": self.reverse.sequence,
            "forward_tm": round(self.forward.tm, 2),
            "reverse_tm": round(self.reverse.tm, 2),
            "amplicon_length": self.amplicon_length,
            "valid": self.forward.is_valid and self.reverse.is_valid,
        }

    def __str__(self):
        return (
            f"Fragment {self.fragment_index:>6} | "
            f"Amplicon: {self.amplicon_length}bp\n"
            f"  {self.forward}\n"
            f"  {self.reverse}"
        )


# ── Thermodynamics ─────────────────────────────────────────────────────────────

def calc_tm_wallace(seq: str) -> float:
    """
    Wallace rule: Tm = 2(A+T) + 4(G+C)
    Valid for sequences 14–20 bp in ~50mM salt.
    """
    seq = seq.upper()
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc


def calc_gc(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq.upper()))


# ── Validation ─────────────────────────────────────────────────────────────────

def validate_primer(seq: str, direction: str) -> Primer:
    warnings = []
    seq = seq.upper()

    tm = calc_tm_wallace(seq)
    gc = calc_gc(seq)

    if tm < TM_MIN:
        warnings.append(f"Tm too low ({tm:.1f}°C < {TM_MIN}°C)")
    if tm > TM_MAX:
        warnings.append(f"Tm too high ({tm:.1f}°C > {TM_MAX}°C)")
    if gc < GC_MIN:
        warnings.append(f"GC too low ({gc:.0%})")
    if gc > GC_MAX:
        warnings.append(f"GC too high ({gc:.0%})")

    # Check for homopolymer runs
    for base in "ATCG":
        run = base * MAX_HOMOPOLYMER
        if run in seq:
            warnings.append(f"Homopolymer run ({run})")
            break

    # Check 3' end — should end in G or C (GC clamp)
    if seq[-1] not in ("G", "C"):
        warnings.append("No 3' GC clamp")

    # Simple self-complementarity check (last 4 bases vs complement of first 4)
    if seq[-4:] == reverse_complement(seq[:4]):
        warnings.append("Potential 3' hairpin")

    return Primer(sequence=seq, direction=direction, tm=tm, gc_content=gc, warnings=warnings)


# ── Primer Generation ──────────────────────────────────────────────────────────

def _derive_primer_from_seed(seed: str, length: int = PRIMER_LENGTH) -> str:
    """
    Derive a deterministic primer-like sequence from a seed string.
    Uses SHA-256 to generate a reproducible base sequence.
    In real systems, primers are designed against flanking sequences.
    Here we generate address-tag primers that can be prepended to oligos.
    """
    h = hashlib.sha256(seed.encode()).hexdigest()

    # Map hex chars to bases: 0-3→A, 4-7→T, 8-b→C, c-f→G
    base_map = {
        "0": "A", "1": "A", "2": "A", "3": "A",
        "4": "T", "5": "T", "6": "T", "7": "T",
        "8": "C", "9": "C", "a": "C", "b": "C",
        "c": "G", "d": "G", "e": "G", "f": "G",
    }
    seq = "".join(base_map[c] for c in h)[:length]

    # Ensure GC clamp (3' end)
    if seq[-1] not in ("G", "C"):
        seq = seq[:-1] + "C"

    return seq


def design_primer_pair(
    fragment_index: int,
    sequence: str,
    namespace: str = "biovault",
) -> PrimerPair:
    """
    Design a forward/reverse primer pair for a DNA fragment.

    In a real synthesis workflow:
    - Forward primer flanks the 5' end of the fragment
    - Reverse primer is the reverse complement of the 3' end
    - Together they amplify only this specific fragment via PCR

    Here we use deterministic address-tag primers derived from fragment identity.
    """
    fwd_seed = f"{namespace}:fwd:{fragment_index}"
    rev_seed = f"{namespace}:rev:{fragment_index}"

    fwd_seq = _derive_primer_from_seed(fwd_seed)
    rev_seq = _derive_primer_from_seed(rev_seed)

    fwd = validate_primer(fwd_seq, "forward")
    rev = validate_primer(rev_seq, "reverse")

    # Amplicon = primer + sequence + primer (with overhang)
    amplicon_len = PRIMER_LENGTH + len(sequence) + PRIMER_LENGTH

    return PrimerPair(
        fragment_index=fragment_index,
        forward=fwd,
        reverse=rev,
        amplicon_length=amplicon_len,
    )


def design_all_primers(fragments, namespace: str = "biovault") -> List[PrimerPair]:
    """Design primer pairs for all fragments."""
    return [
        design_primer_pair(f.index, f.sequence, namespace)
        for f in fragments
    ]


def primers_to_csv(pairs: List[PrimerPair]) -> str:
    """Export primer pairs as CSV (ready for synthesis order)."""
    lines = ["fragment_index,forward_primer,reverse_primer,fwd_tm,rev_tm,amplicon_bp,valid"]
    for p in pairs:
        d = p.to_dict()
        lines.append(
            f"{d['fragment_index']},"
            f"{d['forward_primer']},"
            f"{d['reverse_primer']},"
            f"{d['forward_tm']},"
            f"{d['reverse_tm']},"
            f"{d['amplicon_length']},"
            f"{d['valid']}"
        )
    return "\n".join(lines)


def primer_stats(pairs: List[PrimerPair]) -> dict:
    valid = sum(1 for p in pairs if p.forward.is_valid and p.reverse.is_valid)
    tms = [p.forward.tm for p in pairs] + [p.reverse.tm for p in pairs]
    return {
        "total_pairs": len(pairs),
        "valid_pairs": valid,
        "invalid_pairs": len(pairs) - valid,
        "avg_tm": sum(tms) / len(tms) if tms else 0,
        "min_tm": min(tms) if tms else 0,
        "max_tm": max(tms) if tms else 0,
    }
