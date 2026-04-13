# 🧬 BioVault — Organic DNA Storage System

> *Encoding data into the fabric of life itself.*

BioVault is an organic data storage system that encodes arbitrary files into synthetic DNA sequences. It implements the full pipeline from binary → DNA → FASTA synthesis format → sequencing simulation → decoding back to the original file — with error correction using Reed-Solomon codes.

---

## 🌐 Why DNA Storage?

| Medium | Density | Longevity | Energy |
|---|---|---|---|
| Hard Drive | ~1 TB/kg | ~10 years | Active cooling |
| Flash SSD | ~10 TB/kg | ~10 years | Passive |
| **DNA** | **~215 PB/gram** | **1,000+ years** | **None (frozen/dried)** |

DNA is nature's data format. It has survived millions of years of evolution with built-in error correction. BioVault bridges the gap between software and synthetic biology.

---

## ✨ Features

- **Binary → DNA Encoding** — 2-bit base mapping (A/T/C/G) with index headers
- **Reed-Solomon ECC** — corrects up to 5 byte errors per fragment (tunable)
- **FASTA Export** — standard bioinformatics format, ready for synthesis services
- **GC Content Analysis** — validates synthesis-readiness of sequences
- **Sequencing Error Simulation** — test your ECC against real-world read errors
- **Checksum Verification** — MD5 per-fragment integrity on decode
- **CLI + Python API** — use as a command-line tool or import as a library

---

## 🚀 Quick Start

```bash
git clone https://github.com/yourusername/biovault
cd biovault
pip install -e ".[dev]"
```

### Encode a file

```bash
python -m biovault encode myfile.pdf --ecc --out myfile.fasta
```

### Decode it back

```bash
python -m biovault decode myfile.fasta --out recovered.pdf
```

### Simulate sequencing errors (1% error rate)

```bash
python -m biovault simulate myfile.fasta --error-rate 0.01
```

### View FASTA stats

```bash
python -m biovault stats myfile.fasta
```

---

## 🔬 How It Works

```
┌─────────────┐
│  Input File │  (any format: PDF, image, text, binary)
└──────┬──────┘
       │  bytes_to_bits()
       ▼
┌─────────────┐
│  Bit String │  10110100 01001101 ...
└──────┬──────┘
       │  bits_to_dna()  [00→A, 01→T, 10→C, 11→G]
       ▼
┌─────────────┐
│ DNA Sequence│  ATCGATCG...  (150 bases per oligo)
└──────┬──────┘
       │  + Index header + MD5 checksum
       ▼
┌─────────────┐
│  FASTA File │  >fragment_000001|len=75|chk=3973a3f0
└──────┬──────┘   ATCGATCGATCGATCG...
       │
       ▼  [Optional: send to Twist Bioscience / IDT for synthesis]
       │
       ▼  [Read back via Oxford Nanopore MinION sequencer]
       │
┌─────────────┐
│  Raw Reads  │  Noisy DNA sequences with ~1% error rate
└──────┬──────┘
       │  Reed-Solomon decode + checksum verify
       ▼
┌─────────────┐
│ Output File │  Byte-perfect reconstruction ✓
└─────────────┘
```

---

## 🧪 Python API

```python
from biovault import encode, decode, to_fasta, parse_fasta

# Encode
with open("myfile.pdf", "rb") as f:
    data = f.read()

fragments = encode(data)
fasta = to_fasta(fragments, filename="myfile.pdf")

print(f"Encoded into {len(fragments)} fragments")
print(f"First sequence: {fragments[0].sequence[:40]}...")

# Decode
parsed = parse_fasta(fasta)
recovered, failed = decode(parsed)

assert recovered == data  # Perfect reconstruction
print(f"Recovered {len(recovered)} bytes. Failed fragments: {failed}")
```

### With Reed-Solomon ECC

```python
from biovault.ecc import apply_ecc_to_fragments, ecc_summary

# Protect fragments before synthesis
protected = apply_ecc_to_fragments(fragments)

stats = ecc_summary()
print(f"Can correct up to {stats['max_correctable_errors']} errors per fragment")
```

### Simulate Sequencing Errors

```python
from biovault import simulate_sequencing_errors

# Simulate 1% base read error (typical for Oxford Nanopore)
noisy = simulate_sequencing_errors(fragments, error_rate=0.01)

# Decode noisy fragments (ECC corrects errors)
recovered, failed = decode(noisy, verify=False)
```

---

## 📁 Project Structure

```
biovault/
├── biovault/
│   ├── __init__.py       # Public API
│   ├── __main__.py       # CLI entry point
│   ├── encoder.py        # Binary → DNA encoding
│   ├── decoder.py        # DNA → Binary decoding
│   └── ecc.py            # Reed-Solomon error correction
├── tests/
│   └── test_biovault.py  # Full test suite (11 tests)
├── examples/
│   └── demo.py           # End-to-end demo
├── setup.py
└── README.md
```

---

## 🧬 FASTA Format

BioVault produces standard FASTA files with metadata headers:

```fasta
; BioVault encoded file: document.pdf
; Total fragments: 142
; Total bases: 42,600

>fragment_000000|len=75|chk=3973a3f0
TAACTCCTTCGGTTTCTCATTGTTTCGATGTAAGCCACAATGAGTGTATCGGTGACTCCT...
>fragment_000001|len=75|chk=5058f1af
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
```

These files can be sent directly to DNA synthesis services like:
- [Twist Bioscience](https://www.twistbioscience.com/) — high-throughput synthesis
- [IDT (Integrated DNA Technologies)](https://www.idtdna.com/) — short oligos
- [Genscript](https://www.genscript.com/) — gene synthesis

---

## 🔧 Configuration

| Parameter | Default | Description |
|---|---|---|
| `chunk_bytes` | 75 | Bytes per fragment (→ 150 bases, synthesis-safe) |
| `ecc_symbols` | 10 | RS ECC symbols (corrects up to 5 errors) |
| `error_rate` | 0.01 | Simulation error rate (1% ≈ Nanopore baseline) |

---

## 🧪 Running Tests

```bash
pytest tests/ -v
```

```
test_encode_small                  PASSED
test_roundtrip_text                PASSED
test_roundtrip_binary              PASSED
test_roundtrip_large               PASSED
test_gc_content_balanced           PASSED
test_fragment_indices_sequential   PASSED
test_fasta_roundtrip               PASSED
test_fasta_file_roundtrip          PASSED
test_simulate_errors_changes       PASSED
test_simulate_zero_error_rate      PASSED
test_checksum_detects_corruption   PASSED

11 passed in 0.08s
```

---

## 🗺️ Roadmap

- [ ] Fountain codes (LT codes) for erasure-tolerant encoding
- [ ] Primer sequence design for selective retrieval ("biological addressing")
- [ ] Compression pre-processing (zstd before encode)
- [ ] GC-balanced encoding that preserves decodability
- [ ] Plasmid insertion guide (E. coli living storage protocol)
- [ ] Nanopore MinION read pipeline integration
- [ ] Web UI for drag-and-drop encode/decode

---

## 📚 References

- Church, G. et al. (2012). *Next-Generation Digital Information Storage in DNA*. Science.
- Goldman, N. et al. (2013). *Towards practical, high-capacity, low-maintenance information storage in synthesized DNA*. Nature.
- Organick, L. et al. (2018). *Random access in large-scale DNA data storage*. Nature Biotechnology.
- JCVI Synthetic Cell (2010) — First synthetic genome booted in a cell.

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for details.

---

*BioVault — because the best storage medium was already inside us.*
