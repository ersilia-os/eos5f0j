#!/usr/bin/env bash
#
# Make the checkmol binary available to the Ersilia model eos5f0j.
#
# Prebuilt, statically linked, FP-exception-masked checkmol binaries for the two
# Linux architectures Ersilia targets are committed next to this script:
#     checkmol-linux-x86_64   checkmol-linux-aarch64
# They are built from the committed Pascal sources:
#     checkmol.pas   (upstream, GPL-3.0; https://homepage.univie.ac.at/norbert.haider/cheminf/)
#     fpufix.pas     (masks FP exceptions so checkmol does not trap with EInvalidOp
#                     on Linux, where Free Pascal unmasks them by default)
#
# main.py resolves the architecture-specific binary directly, so this script is
# not strictly required at runtime. It is kept because install.yml invokes it: on
# Linux it exposes a plain `checkmol` (and a copy on the interpreter's PATH); on
# other platforms, or if no matching prebuilt binary exists, it compiles from the
# committed sources when a Free Pascal compiler (fpc) is available.
#
# To reproduce the vendored binaries, run this script on a Debian 11 (bullseye,
# glibc 2.31) image for each architecture with `fpc` + `binutils` installed.
set -eo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$CODE_DIR/checkmol"
OS="$(uname -s)"
MACH="$(uname -m)"

# 1. Prefer the committed prebuilt binary for this Linux architecture.
PREBUILT=""
if [ "$OS" = "Linux" ]; then
  case "$MACH" in
    x86_64|amd64)  PREBUILT="$CODE_DIR/checkmol-linux-x86_64" ;;
    aarch64|arm64) PREBUILT="$CODE_DIR/checkmol-linux-aarch64" ;;
  esac
fi

install_binary() {  # install_binary <src>
  cp "$1" "$OUT"
  chmod +x "$OUT"
  echo "Installed checkmol -> $OUT"
  # Also expose it on the interpreter's bin so source runs (conda run bash run.sh)
  # find it via os.path.dirname(sys.executable) / PATH.
  local pybin=""
  if [ -n "${CONDA_PREFIX:-}" ] && [ -d "$CONDA_PREFIX/bin" ]; then
    pybin="$CONDA_PREFIX/bin"
  elif command -v python3 >/dev/null 2>&1; then
    pybin="$(cd "$(dirname "$(command -v python3)")" && pwd)"
  fi
  if [ -n "$pybin" ] && [ -w "$pybin" ]; then
    cp "$1" "$pybin/checkmol" && chmod +x "$pybin/checkmol"
    echo "checkmol also installed -> $pybin/checkmol"
  fi
}

if [ -n "$PREBUILT" ] && [ -f "$PREBUILT" ]; then
  install_binary "$PREBUILT"
  exit 0
fi

# 2. Fallback: compile from the committed sources if a compiler is available
#    (e.g. on a developer machine whose platform has no prebuilt binary).
echo "No prebuilt checkmol for $OS-$MACH; attempting to compile from source." >&2
if ! command -v fpc >/dev/null 2>&1; then
  echo "WARNING: 'fpc' (Free Pascal) not found and no prebuilt binary for this" >&2
  echo "         platform. checkmol will be resolved from PATH at runtime." >&2
  exit 0
fi

WORK="$(mktemp -d "${TMPDIR:-$HOME}/checkmol-build.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp "$CODE_DIR/checkmol.pas" "$CODE_DIR/fpufix.pas" "$WORK/"
# Inject the FP-exception-mask unit into checkmol's uses clause, then compile.
sed -i.bak -E 's/^([[:space:]]*SYSUTILS,[[:space:]]*MATH)[[:space:]]*;/\1, fpufix;/' \
  "$WORK/checkmol.pas"
( cd "$WORK" && fpc -Sd -O2 checkmol.pas )
install_binary "$WORK/checkmol"
