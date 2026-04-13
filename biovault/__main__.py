"""
BioVault CLI
Usage:
    python -m biovault encode <file> [--ecc] [--out output.fasta]
    python -m biovault decode <file.fasta> [--out output_file]
    python -m biovault simulate <file.fasta> [--error-rate 0.01]
    python -m biovault stats <file.fasta>
"""

import sys
import os
import json
import argparse

from .encoder import encode, to_fasta
from .decoder import decode, parse_fasta, simulate_sequencing_errors


def cmd_encode(args):
    with open(args.file, "rb") as f:
        data = f.read()

    print(f"[+] Encoding {len(data):,} bytes...")
    fragments = encode(data)

    if args.ecc:
        try:
            from .ecc import apply_ecc_to_fragments, ecc_summary
            fragments = apply_ecc_to_fragments(fragments)
            stats = ecc_summary()
            print(f"[+] ECC applied — can correct up to {stats['max_correctable_errors']} errors/fragment")
        except Exception as e:
            print(f"[WARN] ECC skipped: {e}")

    fasta_out = to_fasta(fragments, filename=os.path.basename(args.file))

    out_path = args.out or (args.file + ".fasta")
    with open(out_path, "w") as f:
        f.write(fasta_out)

    total_bases = sum(len(fr.sequence) for fr in fragments)
    print(f"[+] Encoded into {len(fragments)} fragments ({total_bases:,} bases)")
    print(f"[+] Saved to: {out_path}")

    if args.json:
        meta = {
            "source_file": args.file,
            "source_bytes": len(data),
            "fragments": len(fragments),
            "total_bases": total_bases,
            "ecc_enabled": args.ecc,
        }
        print(json.dumps(meta, indent=2))


def cmd_decode(args):
    with open(args.file, "r") as f:
        fasta_text = f.read()

    fragments = parse_fasta(fasta_text)
    print(f"[+] Loaded {len(fragments)} fragments from {args.file}")

    decoded, failed = decode(fragments, verify=not args.no_verify)

    if failed:
        print(f"[WARN] {len(failed)} fragments failed checksum: {failed}")
    else:
        print(f"[+] All fragments decoded successfully ✓")

    out_path = args.out or "biovault_decoded.bin"
    with open(out_path, "wb") as f:
        f.write(decoded)

    print(f"[+] Decoded {len(decoded):,} bytes → {out_path}")


def cmd_simulate(args):
    with open(args.file, "r") as f:
        fasta_text = f.read()

    fragments = parse_fasta(fasta_text)
    noisy = simulate_sequencing_errors(fragments, error_rate=args.error_rate)

    out_path = args.out or args.file.replace(".fasta", f"_noisy_{args.error_rate}.fasta")

    from .encoder import DNAFragment
    noisy_fasta = "".join(f.to_fasta() for f in noisy)
    with open(out_path, "w") as f:
        f.write(noisy_fasta)

    print(f"[+] Simulated {args.error_rate*100:.1f}% error rate on {len(noisy)} fragments")
    print(f"[+] Noisy FASTA saved to: {out_path}")


def cmd_stats(args):
    with open(args.file, "r") as f:
        fasta_text = f.read()

    from .encoder import gc_content
    fragments = parse_fasta(fasta_text)
    total_bases = sum(len(f.sequence) for f in fragments)
    gc_values = [gc_content(f.sequence) for f in fragments]
    avg_gc = sum(gc_values) / len(gc_values)

    print(f"\n{'='*40}")
    print(f"  BioVault Stats: {args.file}")
    print(f"{'='*40}")
    print(f"  Fragments    : {len(fragments)}")
    print(f"  Total bases  : {total_bases:,}")
    print(f"  Avg GC       : {avg_gc:.1%}")
    print(f"  Min GC       : {min(gc_values):.1%}")
    print(f"  Max GC       : {max(gc_values):.1%}")
    print(f"  Avg length   : {total_bases // len(fragments)} bases/fragment")
    print(f"{'='*40}\n")


def main():
    parser = argparse.ArgumentParser(
        prog="biovault",
        description="🧬 BioVault — Organic DNA Storage System"
    )
    sub = parser.add_subparsers(dest="command")

    # encode
    p_enc = sub.add_parser("encode", help="Encode a file to DNA FASTA")
    p_enc.add_argument("file", help="Input file to encode")
    p_enc.add_argument("--out", help="Output FASTA path")
    p_enc.add_argument("--ecc", action="store_true", help="Apply Reed-Solomon ECC")
    p_enc.add_argument("--json", action="store_true", help="Output metadata as JSON")

    # decode
    p_dec = sub.add_parser("decode", help="Decode a FASTA back to original file")
    p_dec.add_argument("file", help="Input FASTA file")
    p_dec.add_argument("--out", help="Output file path")
    p_dec.add_argument("--no-verify", action="store_true", help="Skip checksum verification")

    # simulate
    p_sim = sub.add_parser("simulate", help="Simulate sequencing errors on a FASTA")
    p_sim.add_argument("file", help="Input FASTA file")
    p_sim.add_argument("--error-rate", type=float, default=0.01, help="Error rate (default: 0.01)")
    p_sim.add_argument("--out", help="Output path for noisy FASTA")

    # stats
    p_stats = sub.add_parser("stats", help="Show stats for a FASTA file")
    p_stats.add_argument("file", help="Input FASTA file")

    args = parser.parse_args()

    if args.command == "encode":
        cmd_encode(args)
    elif args.command == "decode":
        cmd_decode(args)
    elif args.command == "simulate":
        cmd_simulate(args)
    elif args.command == "stats":
        cmd_stats(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
