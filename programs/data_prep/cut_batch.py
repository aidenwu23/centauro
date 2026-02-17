#!/usr/bin/env python3
"""Batch run cut_events over a range of run numbers.

Usage:
  python3 scripts/data/cut_batch.py -l 0 -u 9
  python3 scripts/data/cut_batch.py -l 0 -u 9 --dry-run
  python3 scripts/data/cut_batch.py -l 0 -u 9 --drop-inputs

Assumptions:
  Inputs: data/events/filtered/NNNN_skim.root
  Outputs: data/events/filtered/NNNN_skim_q2y.root
  cut_events binary: build/cut_events relative to repo root
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
  # CLI: define run range, file layout, binary location, and dry-run toggle.
  parser = argparse.ArgumentParser(description="Batch cut EDM4eic files with cut_events.")
  parser.add_argument("-l", "--lower", type=int, required=True, help="First run index (inclusive).")
  parser.add_argument("-u", "--upper", type=int, required=True, help="Last run index (inclusive).")
  parser.add_argument("-w", "--width", type=int, default=4, help="Zero-pad width for run numbers.")
  parser.add_argument(
      "--input-dir",
      default="data/events/filtered",
      help="Directory holding input skims (default: data/events/filtered).")
  parser.add_argument(
      "--cut-binary",
      default="build/cut_events",
      help="Path to cut_events binary relative to repo root (default: build/cut_events).")
  parser.add_argument(
      "--drop-inputs",
      action="store_true",
      help="Delete the input skim after a successful cut (default: keep inputs).")
  parser.add_argument(
      "--dry-run",
      action="store_true",
      help="Print planned commands without running them.")
  parser.add_argument(
      "--cut-args",
      default="",
      help="Extra args to pass to cut_events (quoted string).")
  return parser.parse_args()


def main() -> int:
  args = parse_args()
  if args.upper < args.lower:
    print("[cut_batch] upper must be >= lower", file=sys.stderr)
    return 1

  # Resolve repo root and paths so the loop can stay simple.
  root = Path(__file__).resolve().parent.parent.parent
  cut_bin = root / args.cut_binary
  if not cut_bin.is_file():
    print(f"[cut_batch] cut_events not found at {cut_bin}", file=sys.stderr)
    return 1

  in_dir = root / args.input_dir
  if not in_dir.is_dir():
    print(f"[cut_batch] input directory not found: {in_dir}", file=sys.stderr)
    return 1

  processed = 0
  skipped = 0
  failed = 0

  for idx in range(args.lower, args.upper + 1):
    # Build input/output names for this run index.
    tag = f"{idx:0{args.width}d}"
    in_file = in_dir / f"{tag}_skim.root"
    out_file = in_dir / f"{tag}_skim_q2y.root"

    if out_file.exists():
      print(f"[skip] {out_file.name} exists; skipping")
      skipped += 1
      continue

    if not in_file.is_file():
      print(f"[warn] input missing: {in_file.name}; skipping", file=sys.stderr)
      failed += 1
      continue

    cmd = [str(cut_bin), "-i", str(in_file), "-o", str(out_file)]
    if args.cut_args:
      cmd.extend(args.cut_args.split())
    if args.dry_run:
      print(f"[plan] {' '.join(cmd)}")
      skipped += 1
      continue

    print(f"[cut] {' '.join(cmd)}")
    try:
      subprocess.run(cmd, check=True)
      processed += 1
      if args.drop_inputs:
        try:
          in_file.unlink()
          print(f"[cleanup] deleted {in_file.name}")
        except OSError as exc:
          print(f"[cleanup] warning: failed to delete {in_file.name}: {exc}", file=sys.stderr)
    except subprocess.CalledProcessError as exc:
      print(f"[cut] failed (exit {exc.returncode}) for {in_file.name}", file=sys.stderr)
      failed += 1

  print(f"[cut_batch] done processed={processed} skipped={skipped} failed={failed} dry_run={args.dry_run}")
  return 0 if failed == 0 else 1


if __name__ == "__main__":
  sys.exit(main())
