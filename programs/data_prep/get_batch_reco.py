#!/usr/bin/env python3
"""Download and reduce EDM4eic ROOT files in numbered batches.

Features:
  - Takes a numeric lower/upper range and processes in chunks (default 5 files).
  - Downloads inputs into data/events/raw preserving the remote filename.
  - Reduces to data/events/filtered/NNNN_skim.root using the build/reducer binary.
  - Skips reduction when the processed file already exists; downloads a missing
    original even if the skim is present.
  - Optional removal of originals when skims are available.

The remote template may contain {padded} (zero-padded index) and {idx} (raw int).
Default template matches the JLab numbered files (5.0000, 5.0001, ...):
  root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.0/epic_craterlake/DIS/NC/18x275/minQ2=100/pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_5.{padded}.eicrecon.edm4eic.root

python scripts/data/get_batch_reco.py -l 0 -u 9 -c 1 --drop-inputs
python scripts/data/get_batch_reco.py -l 2 -u 2 -c 1


# 25.08.0
python scripts/data/get_batch_reco.py -l 0 -u 0 -c 1 --remote-template "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.08.0/epic_craterlake/DIS/NC/18x275/minQ2=100/pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_5.{padded}.eicrecon.edm4eic.root"

python scripts/data/get_batch_reco.py --lower 0 -upper 19 --chunk-size 5 \
    --drop-inputs --dry-run
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from itertools import islice
from pathlib import Path
from typing import Iterable, List


# Yield input values in fixed-size groups.
def chunked(seq: Iterable[int], size: int) -> Iterable[List[int]]:
  iterator = iter(seq)
  while True:
    block = list(islice(iterator, size))
    if not block:
      break
    yield block


# Collect CLI settings for the batch run.
def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description="Batch download + reduce EDM4eic ROOT files.")
  parser.add_argument(
      "-l",
      "--lower",
      type=int,
      required=True,
      help="First index to process (inclusive).")
  parser.add_argument(
      "-u",
      "--upper",
      type=int,
      required=True,
      help="Last index to process (inclusive).")
  parser.add_argument(
      "-c",
      "--chunk-size",
      type=int,
      default=5,
      help="How many files to process per chunk (download + reduce). Default: 5.")
  parser.add_argument(
      "-w",
      "--width",
      type=int,
      default=4,
      help="Zero-padding width for file indices. Default: 4 (0000, 0001, ...).")
  parser.add_argument(
      "--remote-template",
      default=("root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.0/epic_craterlake/"
               "DIS/NC/18x275/minQ2=100/pythia8NCDIS_18x275_minQ2=100_"
               "beamEffects_xAngle=-0.025_hiDiv_5.{padded}.eicrecon.edm4eic.root"),
      help="Template for remote files; supports {padded} and {idx}.")
  parser.add_argument(
      "--drop-inputs",
      action="store_true",
      help="Do not keep original inputs; remove them when a reduced file is available.")
  parser.add_argument(
      "--dry-run",
      action="store_true",
      help="Print planned actions without downloading or reducing.")
  return parser.parse_args()


# Pad file indices for consistent filenames.
def format_index(idx: int, width: int) -> str:
  return f"{idx:0{width}d}"


# Download a single remote file; respect dry-run.
def download_one(remote: str, dest: Path, dry_run: bool) -> bool:
  if dry_run:
    print(f"[download] {remote} -> {dest} (dry-run)")
    return True
  print(f"[download] {remote} -> {dest}")
  try:
    subprocess.run(["xrdcp", remote, str(dest)], check=True)
    return True
  except subprocess.CalledProcessError as exc:
    print(f"[download] xrdcp failed (exit {exc.returncode}) for {remote}", file=sys.stderr)
    return False


# Run reducer binary for one input/output pair.
def reduce_one(src: Path, dst: Path, reducer: Path, dry_run: bool) -> bool:
  cmd = [str(reducer), "-i", str(src), "-o", str(dst)]
  if dry_run:
    print(f"[reduce] {' '.join(cmd)} (dry-run)")
    return True
  print(f"[reduce] {' '.join(cmd)}")
  try:
    subprocess.run(cmd, check=True)
    return True
  except subprocess.CalledProcessError as exc:
    print(f"[reduce] reducer failed (exit {exc.returncode}) for {src}", file=sys.stderr)
    return False


def trim_inputs(orig_dir: Path, processed_dir: Path, width: int) -> None:
  # Remove any originals that already have a processed partner.
  tags = set()
  for processed in processed_dir.glob("*_skim.root"):
    name = processed.name
    if not name.endswith("_skim.root"):
      continue
    tag = name[:-len("_skim.root")]
    if tag and (width == 0 or len(tag) == width):
      tags.add(tag)

  for orig in orig_dir.iterdir():
    if not orig.is_file():
      continue
    for tag in tags:
      if tag in orig.name:
        try:
          orig.unlink()
          print(f"[cleanup] Deleted stale input {orig.name}")
        except OSError as exc:
          print(f"[cleanup] Warning: failed to delete {orig}: {exc}", file=sys.stderr)
        break


def main() -> int:
  args = parse_args()
  # Validate numeric inputs early.
  if args.upper < args.lower:
    print("[batch] upper must be >= lower", file=sys.stderr)
    return 1
  if args.chunk_size <= 0:
    print("[batch] chunk-size must be positive", file=sys.stderr)
    return 1

  # Resolve repository paths for binaries and data.
  chunk = args.chunk_size
  tmpl = args.remote_template
  root_dir = Path(__file__).resolve().parent.parent.parent
  orig_dir = root_dir / "data" / "events" / "raw"
  processed_dir = root_dir / "data" / "events" / "filtered"
  reducer = root_dir / "build" / "reducer"
  if not reducer.is_file():
    print("[batch] No reducer binary found at build/reducer. Build it first.", file=sys.stderr)
    return 1

  orig_dir.mkdir(parents=True, exist_ok=True)
  processed_dir.mkdir(parents=True, exist_ok=True)

  indices = range(args.lower, args.upper + 1)
  # Track work for the final summary.
  processed = 0
  skipped = 0

  has_idx = "{padded}" in tmpl or "{idx}" in tmpl
  # Warn when template reuses a single source for all indices.
  if not has_idx and args.upper > args.lower:
    print("[batch] Warning: remote-template has no {padded}/{idx}; same source will be used for all indices.")

  # Process downloads and reductions in manageable chunks.
  for block in chunked(indices, chunk):
    for idx in block:
      tag = format_index(idx, args.width)
      print(f" ")
      print(f"[batch] Processing {tag}")
      remote = tmpl.format(padded=tag, idx=idx) if has_idx else tmpl
      remote_name = Path(remote).name
      in_file = orig_dir / remote_name
      out_file = processed_dir / f"{tag}_skim.root"
      out_exists = out_file.exists()
      in_exists = in_file.exists()

      # In drop-inputs mode we avoid keeping originals whenever a skim exists.
      if args.drop_inputs:
        if out_exists:
          if in_exists and not args.dry_run:
            try:
              in_file.unlink()
              print(f"[cleanup] Deleted {in_file.name}")
            except OSError as exc:
              print(f"[cleanup] Warning: failed to delete {in_file}: {exc}", file=sys.stderr)
          skipped += 1
          continue
        if not in_exists:
          if not download_one(remote, in_file, args.dry_run):
            print(f"[batch] Aborting due to download failure for index {idx}", file=sys.stderr)
            return 1
        if not reduce_one(in_file, out_file, reducer, args.dry_run):
          print(f"[batch] Aborting due to reduce failure for index {idx}", file=sys.stderr)
          return 1
        processed += 1
        if not args.dry_run and in_file.exists():
          try:
            in_file.unlink()
            print(f"[cleanup] Deleted {in_file.name}")
          except OSError as exc:
            print(f"[cleanup] Warning: failed to delete {in_file}: {exc}", file=sys.stderr)
        continue

      if out_exists:
        print(f"[skip] Reduced exists: {out_file.name}; skipping reduction")
        if not in_exists:
          if not download_one(remote, in_file, args.dry_run):
            print(f"[batch] Aborting due to download failure for index {idx}", file=sys.stderr)
            return 1
        skipped += 1
        continue

      if not in_exists:
        if not download_one(remote, in_file, args.dry_run):
          print(f"[batch] Aborting due to download failure for index {idx}", file=sys.stderr)
          return 1

      if not reduce_one(in_file, out_file, reducer, args.dry_run):
        print(f"[batch] Aborting due to reduce failure for index {idx}", file=sys.stderr)
        return 1

      processed += 1

  if args.drop_inputs and not args.dry_run:
    # Sweep any other inputs that already have matching skims.
    trim_inputs(orig_dir, processed_dir, args.width)

  print(f"")
  print(f"[batch] Done. processed={processed} skipped={skipped} dry_run={args.dry_run}")
  return 0


if __name__ == "__main__":
  sys.exit(main())
