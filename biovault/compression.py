"""
BioVault Compression
Pre-compression before DNA encoding to minimize synthesis cost.
Fewer bytes → fewer bases → cheaper synthesis.

Strategy: Try zstd, gzip, and raw — pick smallest output.
"""

import gzip
import zlib
from enum import Enum
from dataclasses import dataclass
from typing import Tuple

try:
    import zstandard as zstd
    ZSTD_AVAILABLE = True
except ImportError:
    ZSTD_AVAILABLE = False


MAGIC_ZSTD = b"\xBB\x01"   # BioVault zstd marker
MAGIC_GZIP = b"\xBB\x02"   # BioVault gzip marker
MAGIC_RAW  = b"\xBB\x00"   # BioVault raw (no compression)
MAGIC_LEN  = 2


class Codec(Enum):
    RAW  = "raw"
    GZIP = "gzip"
    ZSTD = "zstd"


@dataclass
class CompressionResult:
    data: bytes
    codec: Codec
    original_size: int
    compressed_size: int

    @property
    def ratio(self) -> float:
        return self.compressed_size / self.original_size if self.original_size else 1.0

    @property
    def savings_pct(self) -> float:
        return (1 - self.ratio) * 100

    def __str__(self):
        return (
            f"Codec: {self.codec.value} | "
            f"{self.original_size:,}B → {self.compressed_size:,}B "
            f"({self.savings_pct:.1f}% saved)"
        )


def _try_zstd(data: bytes, level: int = 3) -> bytes:
    if not ZSTD_AVAILABLE:
        return data
    cctx = zstd.ZstdCompressor(level=level)
    return cctx.compress(data)


def _try_gzip(data: bytes, level: int = 6) -> bytes:
    return gzip.compress(data, compresslevel=level)


def compress(data: bytes) -> CompressionResult:
    """
    Compress data using the best available codec.
    Returns CompressionResult with the smallest output and a magic header.
    """
    candidates = [(Codec.RAW, data)]

    gz = _try_gzip(data)
    candidates.append((Codec.GZIP, gz))

    if ZSTD_AVAILABLE:
        zd = _try_zstd(data)
        candidates.append((Codec.ZSTD, zd))

    # Pick winner (smallest)
    best_codec, best_data = min(candidates, key=lambda x: len(x[1]))

    # Prepend magic header
    magic = {Codec.RAW: MAGIC_RAW, Codec.GZIP: MAGIC_GZIP, Codec.ZSTD: MAGIC_ZSTD}
    framed = magic[best_codec] + best_data

    return CompressionResult(
        data=framed,
        codec=best_codec,
        original_size=len(data),
        compressed_size=len(framed),
    )


def decompress(data: bytes) -> bytes:
    """
    Decompress BioVault-framed data. Auto-detects codec from magic header.
    """
    if len(data) < MAGIC_LEN:
        raise ValueError("Data too short to contain BioVault magic header")

    magic = data[:MAGIC_LEN]
    payload = data[MAGIC_LEN:]

    if magic == MAGIC_RAW:
        return payload
    elif magic == MAGIC_GZIP:
        return gzip.decompress(payload)
    elif magic == MAGIC_ZSTD:
        if not ZSTD_AVAILABLE:
            raise RuntimeError("zstandard not installed — cannot decompress zstd payload")
        dctx = zstd.ZstdDecompressor()
        return dctx.decompress(payload)
    else:
        # No magic header — assume raw (legacy support)
        return data


def compression_report(data: bytes) -> dict:
    """Benchmark all codecs on the given data."""
    report = {"original_bytes": len(data), "codecs": {}}

    for codec, fn in [("gzip", _try_gzip), ("zlib", zlib.compress)]:
        compressed = fn(data)
        report["codecs"][codec] = {
            "size": len(compressed),
            "ratio": len(compressed) / len(data),
            "savings_pct": (1 - len(compressed) / len(data)) * 100,
        }

    if ZSTD_AVAILABLE:
        for level in [1, 3, 9, 19]:
            compressed = _try_zstd(data, level=level)
            report["codecs"][f"zstd_l{level}"] = {
                "size": len(compressed),
                "ratio": len(compressed) / len(data),
                "savings_pct": (1 - len(compressed) / len(data)) * 100,
            }

    return report
