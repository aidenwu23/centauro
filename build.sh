#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
MODULE="centauro"

usage() {
  cat <<EOF
Usage: $0 [clean] [fast] [cmake-opts...]

Options:
  clean   Remove the build directory before configuring
  fast    Skip configure, just build/install

Env:
  EIC_SHELL_PREFIX       Hint for dependencies (EIC env)
  BUILD_NPROC            Parallel build jobs (default: nproc or 1)
  CENTAURO_BUILD_TYPE    CMAKE_BUILD_TYPE (default: RelWithDebInfo)
  CENTAURO_WITH_EDM      ON/OFF EDM helpers (default: ON)
  CENTAURO_INSTALL_PREFIX Install prefix (default: build/install)
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

clean=0
fast=0
extra_opts=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    clean) clean=1 ;;
    fast) fast=1 ;;
    *) extra_opts+=("$1") ;;
  esac
  shift
done

install_prefix="${CENTAURO_INSTALL_PREFIX:-${BUILD_DIR}/install}"
nproc="${BUILD_NPROC:-}"
if [[ -z "$nproc" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    nproc="$(nproc)"
  else
    nproc=1
  fi
fi

# Build a CMake-friendly prefix path (semicolon-separated).
existing_prefix_path="${CMAKE_PREFIX_PATH:-}"
existing_prefix_path="${existing_prefix_path//:/;}"
cmake_prefix_entries=()
[[ -n "${EIC_SHELL_PREFIX:-}" ]] && cmake_prefix_entries+=("$EIC_SHELL_PREFIX")
# Keep local install ahead of system hints so we stay isolated.
cmake_prefix_entries+=("$install_prefix" "/opt/local")
[[ -n "$existing_prefix_path" ]] && cmake_prefix_entries+=("$existing_prefix_path")
export CMAKE_PREFIX_PATH="$(IFS=';'; echo "${cmake_prefix_entries[*]}")"

build_type="${CENTAURO_BUILD_TYPE:-RelWithDebInfo}"
with_edm="${CENTAURO_WITH_EDM:-ON}"

cmake_gen=(
  cmake
  -S "$ROOT_DIR"
  -B "$BUILD_DIR"
  -DCMAKE_BUILD_TYPE="$build_type"
  -DCMAKE_INSTALL_PREFIX="$install_prefix"
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"
  -DCMAKE_FIND_DEBUG_MODE=OFF
  -DCENTAURO_WITH_EDM="$with_edm"
)
cmake_build=(cmake --build "$BUILD_DIR" -j"$nproc")
cmake_install=(cmake --install "$BUILD_DIR")

cat <<EOF

BUILDING:
module    = $MODULE
clean     = $clean
fast      = $fast
prefix    = $install_prefix
nproc     = $nproc
extraOpts = ${extra_opts[*]:-(none)}

CMAKE COMMANDS:
========================================
[ generate buildsystem ]
${cmake_gen[*]} ${extra_opts[*]}

[ build ]
${cmake_build[*]}

[ install ]
${cmake_install[*]}
========================================
EOF

if [[ $clean -eq 1 ]]; then
  echo "--- removing $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

[[ $fast -eq 0 ]] && "${cmake_gen[@]}" "${extra_opts[@]}"
"${cmake_build[@]}"
"${cmake_install[@]}"
