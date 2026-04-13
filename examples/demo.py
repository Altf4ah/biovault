"""
BioVault Example: Encode → Simulate Errors → Decode
Shows the full pipeline including error simulation.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biovault.encoder import encode, to_fasta, gc_content
from biovault.decoder import decode, parse_fasta, simulate_sequencing_errors


def main():
    # ── Step 1: Original data ──────────────────────────────
    original = b"BioVault: storing data in the fabric of life itself. " * 5
    print(f"Original data: {len(original)} bytes")

    # ── Step 2: Encode to DNA ──────────────────────────────
    fragments = encode(original)
    total_bases = sum(len(f.sequence) for f in fragments)
    avg_gc = sum(gc_content(f.sequence) for f in fragments) / len(fragments)

    print(f"\n[ENCODE]")
    print(f"  Fragments  : {len(fragments)}")
    print(f"  Total bases: {total_bases:,}")
    print(f"  Avg GC     : {avg_gc:.1%}")
    print(f"  First seq  : {fragments[0].sequence[:40]}...")

    # ── Step 3: Export as FASTA ────────────────────────────
    fasta = to_fasta(fragments, filename="example.txt")
    print(f"\n[FASTA PREVIEW]\n{fasta[:300]}...")

    # ── Step 4: Simulate sequencing errors (1% rate) ───────
    noisy = simulate_sequencing_errors(fragments, error_rate=0.01)
    diffs = sum(
        sum(a != b for a, b in zip(o.sequence, n.sequence))
        for o, n in zip(fragments, noisy)
    )
    print(f"\n[SIMULATE ERRORS]")
    print(f"  Error rate : 1%")
    print(f"  Total mutations introduced: {diffs}")

    # ── Step 5: Decode clean version ──────────────────────
    parsed = parse_fasta(fasta)
    recovered, failed = decode(parsed)

    print(f"\n[DECODE]")
    print(f"  Failed fragments: {len(failed)}")
    print(f"  Recovered bytes : {len(recovered)}")
    print(f"  Match           : {'✓ PERFECT' if recovered == original else '✗ MISMATCH'}")

    print(f"\n[SAMPLE OUTPUT]")
    print(f"  {recovered[:80].decode()}")


if __name__ == "__main__":
    main()
