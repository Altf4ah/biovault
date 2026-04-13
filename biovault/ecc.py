"""
BioVault Error Correction
Reed-Solomon ECC layer for DNA storage — corrects synthesis and sequencing errors.
DNA sequencing typically has ~0.1–1% error rate per base.
"""

from typing import List, Tuple
from .encoder import DNAFragment


try:
    from reedsolo import RSCodec, ReedSolomonError
    RS_AVAILABLE = True
except ImportError:
    RS_AVAILABLE = False


ECC_SYMBOLS = 10  # Can correct up to 5 byte errors per block


class ECCNotAvailable(Exception):
    pass


def _check_rs():
    if not RS_AVAILABLE:
        raise ECCNotAvailable(
            "reedsolo not installed. Run: pip install reedsolo"
        )


def encode_with_ecc(data: bytes, ecc_symbols: int = ECC_SYMBOLS) -> bytes:
    """Add Reed-Solomon error correction codes to raw bytes."""
    _check_rs()
    rs = RSCodec(ecc_symbols)
    return bytes(rs.encode(data))


def decode_with_ecc(data: bytes, ecc_symbols: int = ECC_SYMBOLS) -> Tuple[bytes, bool]:
    """
    Decode Reed-Solomon protected bytes.

    Returns:
        (decoded_bytes, was_corrected: bool)
    """
    _check_rs()
    rs = RSCodec(ecc_symbols)
    try:
        decoded, _, errata = rs.decode(data)
        corrected = len(errata) > 0
        return bytes(decoded), corrected
    except ReedSolomonError as e:
        raise ValueError(f"Uncorrectable errors in fragment: {e}")


def apply_ecc_to_fragments(
    fragments: List[DNAFragment],
) -> List[DNAFragment]:
    """
    Re-encode fragment sequences with ECC protection.
    Converts: sequence → bytes → RS-encoded bytes → new DNA sequence
    """
    from .encoder import bytes_to_bits, bits_to_dna, balance_gc, fragment_checksum

    _check_rs()
    protected = []

    for frag in fragments:
        # Convert DNA back to raw bytes
        from .decoder import dna_to_bits, bits_to_bytes
        bits = dna_to_bits(frag.sequence)
        raw = bits_to_bytes(bits, frag.original_length)

        # Add ECC
        ecc_data = encode_with_ecc(raw)

        # Re-encode as DNA
        new_bits = bytes_to_bits(ecc_data)
        new_seq = balance_gc(bits_to_dna(new_bits))

        protected.append(
            DNAFragment(
                index=frag.index,
                sequence=new_seq,
                checksum=fragment_checksum(raw),  # checksum on original data
                original_length=frag.original_length,
            )
        )

    return protected


def ecc_summary(ecc_symbols: int = ECC_SYMBOLS) -> dict:
    """Return ECC capability stats."""
    return {
        "ecc_symbols": ecc_symbols,
        "max_correctable_errors": ecc_symbols // 2,
        "overhead_bytes": ecc_symbols,
        "overhead_bases": ecc_symbols * 4,  # 1 byte = 4 bases
    }
