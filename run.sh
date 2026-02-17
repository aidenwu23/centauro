#!/usr/bin/env bash
# Run a previously built executable (typically under ./build) with optional arguments.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"

die() { echo "[error] $*" >&2; exit 1; }
info() { echo "[info] $*"; }

[[ $# -ge 1 ]] || die "Usage: $0 <path-to-binary-or-source> [-- program args]"

INPUT="$1"
shift || true

resolve_exec() {
  local target="$1"

  # Special-case cross-frame splitting (name collision with scripts/splitting.cc)
  if [[ "$target" == */scripts/cross_frame/splitting.cc || "$target" == "scripts/cross_frame/splitting.cc" || "$target" == "cross_frame/splitting.cc" ]]; then
    local candidate="$BUILD_DIR/cf_splitting"
    if [[ -x "$candidate" ]]; then
      echo "$candidate"; return 0
    fi
  fi

  # Special-case same-frame matcher since executable name differs from source basename.
  if [[ "$target" == */programs/same_frame/match_jets.cc || "$target" == "programs/same_frame/match_jets.cc" || "$target" == "same_frame/match_jets.cc" ]]; then
    local candidate="$BUILD_DIR/sf_match_jets"
    if [[ -x "$candidate" ]]; then
      echo "$candidate"; return 0
    fi
  fi

  # Special-case diff-frame matcher since executable name differs from source basename.
  if [[ "$target" == */programs/diff_frame/match_jets.cc || "$target" == "programs/diff_frame/match_jets.cc" || "$target" == "diff_frame/match_jets.cc" ]]; then
    local candidate="$BUILD_DIR/df_match_jets"
    if [[ -x "$candidate" ]]; then
      echo "$candidate"; return 0
    fi
  fi

  # Direct executable path as given
  if [[ -x "$target" && -f "$target" ]]; then
    echo "$target"; return 0
  fi

  # Executable relative to repo root
  if [[ -x "$ROOT_DIR/$target" && -f "$ROOT_DIR/$target" ]]; then
    echo "$ROOT_DIR/$target"; return 0
  fi

  # If user passed a source-like path, try the build output (basename, no extension)
  if [[ -f "$target" ]]; then
    local base; base="$(basename "${target%.*}")"
    local candidate="$BUILD_DIR/$base"
    if [[ -x "$candidate" ]]; then
      echo "$candidate"; return 0
    fi
  fi

  # Last resort: look for basename (with or without extension) inside build dir
  local base_input; base_input="$(basename "${target%.*}")"
  if [[ -x "$BUILD_DIR/$base_input" ]]; then
    echo "$BUILD_DIR/$base_input"; return 0
  fi

  return 1
}

if ! EXEC_PATH="$(resolve_exec "$INPUT")"; then
  die "Could not find executable for '$INPUT'. Build it first with ./build.sh"
fi

info "Running: $EXEC_PATH $*"
"$EXEC_PATH" "$@"
rc=$?
info "Process exited with status $rc"
exit $rc
