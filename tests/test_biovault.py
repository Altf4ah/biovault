"""
BioVault Test Suite
"""

import pytest
import os
import tempfile
from biovault.encoder import encode, to_fasta, gc_content
from biovault.decoder import decode, parse_fasta, simulate_sequencing_errors


# ── Encoder Tests ──────────────────────────────────────────

def test_encode_small():
    data = b"Hello, BioVault!"
    fragments = encode(data)
    assert len(fragments) > 0
    assert all(set(f.sequence).issubset(set("ATCG")) for f in fragments)


def test_roundtrip_text():
    original = b"BioVault is a DNA storage system. This is a test of the encoder and decoder."
    fragments = encode(original)
    recovered, failed = decode(fragments)
    assert failed == []
    assert recovered == original


def test_roundtrip_binary():
    original = os.urandom(512)
    fragments = encode(original)
    recovered, failed = decode(fragments)
    assert failed == []
    assert recovered == original


def test_roundtrip_large():
    original = os.urandom(4096)
    fragments = encode(original)
    recovered, failed = decode(fragments)
    assert failed == []
    assert recovered == original


def test_gc_content_balanced():
    data = os.urandom(256)
    fragments = encode(data)
    for f in fragments:
        gc = gc_content(f.sequence)
        assert 0.30 <= gc <= 0.70, f"GC content {gc:.2f} out of acceptable range"


def test_fragment_indices_sequential():
    data = os.urandom(500)
    fragments = encode(data)
    for i, f in enumerate(fragments):
        assert f.index == i


# ── FASTA Tests ────────────────────────────────────────────

def test_fasta_roundtrip():
    original = b"FASTA roundtrip test data!"
    fragments = encode(original)
    fasta = to_fasta(fragments)

    assert ">fragment_" in fasta
    assert "len=" in fasta
    assert "chk=" in fasta

    parsed = parse_fasta(fasta)
    recovered, failed = decode(parsed)
    assert failed == []
    assert recovered == original


def test_fasta_file_roundtrip():
    original = b"Testing file I/O for BioVault FASTA export."
    fragments = encode(original)
    fasta = to_fasta(fragments, filename="test.txt")

    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as f:
        f.write(fasta)
        tmp_path = f.name

    try:
        with open(tmp_path) as f:
            loaded = f.read()
        parsed = parse_fasta(loaded)
        recovered, failed = decode(parsed)
        assert failed == []
        assert recovered == original
    finally:
        os.unlink(tmp_path)


# ── Simulation Tests ───────────────────────────────────────

def test_simulate_errors_changes_sequence():
    data = b"Test error simulation input data for BioVault."
    fragments = encode(data)
    noisy = simulate_sequencing_errors(fragments, error_rate=0.05)

    # At 5% error rate, sequences should differ
    original_seqs = [f.sequence for f in fragments]
    noisy_seqs = [f.sequence for f in noisy]
    assert original_seqs != noisy_seqs


def test_simulate_zero_error_rate():
    data = b"No errors should appear."
    fragments = encode(data)
    noisy = simulate_sequencing_errors(fragments, error_rate=0.0)
    for orig, noisy_frag in zip(fragments, noisy):
        assert orig.sequence == noisy_frag.sequence


# ── Compression Tests ──────────────────────────────────────

def test_compression_roundtrip_text():
    from biovault.compression import compress, decompress
    data = b"Repetitive text for compression testing. " * 50
    result = compress(data)
    recovered = decompress(result.data)
    assert recovered == data

def test_compression_roundtrip_binary():
    from biovault.compression import compress, decompress
    data = os.urandom(512)
    result = compress(data)
    recovered = decompress(result.data)
    assert recovered == data

def test_compression_savings_on_repetitive():
    from biovault.compression import compress
    data = b"AAAA" * 1000
    result = compress(data)
    assert result.compressed_size < result.original_size

def test_compression_no_inflate_on_random():
    from biovault.compression import compress
    data = os.urandom(256)
    result = compress(data)
    # Random data may not compress — should never inflate more than magic header
    assert result.compressed_size <= result.original_size + 10


# ── Primer Tests ───────────────────────────────────────────

def test_primer_design_returns_pairs():
    from biovault.primers import design_all_primers
    data = os.urandom(300)
    fragments = encode(data)
    pairs = design_all_primers(fragments)
    assert len(pairs) == len(fragments)

def test_primer_sequences_are_dna():
    from biovault.primers import design_all_primers
    data = b"Primer sequence validation test data here."
    fragments = encode(data)
    pairs = design_all_primers(fragments)
    for pair in pairs:
        assert set(pair.forward.sequence).issubset(set("ATCG"))
        assert set(pair.reverse.sequence).issubset(set("ATCG"))

def test_primer_gc_clamp():
    from biovault.primers import design_all_primers
    data = b"GC clamp test for primer 3-prime end validation."
    fragments = encode(data)
    pairs = design_all_primers(fragments)
    for pair in pairs:
        assert pair.forward.sequence[-1] in ("G", "C")
        assert pair.reverse.sequence[-1] in ("G", "C")

def test_primer_csv_export():
    from biovault.primers import design_all_primers, primers_to_csv
    data = b"CSV export test for primer pairs in BioVault."
    fragments = encode(data)
    pairs = design_all_primers(fragments)
    csv = primers_to_csv(pairs)
    assert "fragment_index" in csv
    assert "forward_primer" in csv
    lines = csv.strip().split("\n")
    assert len(lines) == len(pairs) + 1  # header + data rows

def test_full_pipeline_with_compression():
    from biovault.compression import compress, decompress
    original = b"Full pipeline test with compression enabled. " * 20
    compressed = compress(original)
    fragments = encode(compressed.data)
    recovered_compressed, failed = decode(fragments)
    assert failed == []
    recovered = decompress(recovered_compressed)
    assert recovered == original


def test_checksum_detects_corruption():
    from biovault.decoder import ChecksumError, decode_fragment
    from biovault.encoder import DNAFragment

    data = b"Integrity check test."
    fragments = encode(data)
    frag = fragments[0]

    # Corrupt the sequence
    corrupted_seq = "A" * len(frag.sequence)
    bad_frag = DNAFragment(
        index=frag.index,
        sequence=corrupted_seq,
        checksum=frag.checksum,
        original_length=frag.original_length,
    )

    with pytest.raises((ChecksumError, KeyError, ValueError)):
        decode_fragment(bad_frag, verify=True)
