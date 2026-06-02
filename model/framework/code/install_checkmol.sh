#!/usr/bin/env bash
#
# Build the checkmol binary for the Ersilia model eos5f0j.
#
# This script is self-contained and resolves its own location, so it can be
# invoked from anywhere. It:
#   1. uses an existing Free Pascal compiler (fpc) if one is on PATH;
#   2. otherwise downloads the official Free Pascal release for the current
#      architecture and installs it into a private prefix (no root required);
#   3. downloads the checkmol Pascal source and compiles it next to main.py.
#
# checkmol is GPL-3.0; see https://homepage.univie.ac.at/norbert.haider/cheminf/
set -eo pipefail

FPC_VERSION="3.2.2"
CHECKMOL_URL="https://homepage.univie.ac.at/norbert.haider/download/chemistry/checkmol/checkmol.pas"

# 1. Locate this script -> the model code directory (the compile target). main.py
#    looks for the checkmol binary next to itself, so we build it right here.
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$CODE_DIR/checkmol"

# 2. Scratch space. Do not assume /tmp exists: honour $TMPDIR, else fall back to
#    $HOME (which is always present), and let mktemp create a unique directory.
WORK="$(mktemp -d "${TMPDIR:-$HOME}/checkmol-build.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# Small download helper: prefer wget, fall back to curl.
fetch() {  # fetch <url> <dest>
  if command -v wget >/dev/null 2>&1; then
    wget -q -O "$2" "$1"
  elif command -v curl >/dev/null 2>&1; then
    curl -fsSL -o "$2" "$1"
  else
    echo "ERROR: neither wget nor curl is available to download files" >&2
    exit 1
  fi
}

# 3. Resolve a Free Pascal compiler.
if command -v fpc >/dev/null 2>&1; then
  FPC_BIN="$(command -v fpc)"
else
  OS="$(uname -s)"
  MACH="$(uname -m)"
  if [ "$OS" != "Linux" ]; then
    echo "ERROR: 'fpc' is not on PATH and automatic Free Pascal install is only" >&2
    echo "       supported on Linux (detected '$OS'). Please install Free Pascal." >&2
    exit 1
  fi
  case "$MACH" in
    x86_64|amd64)  FPC_ARCH="x86_64";  PPC="ppcx64" ;;
    aarch64|arm64) FPC_ARCH="aarch64"; PPC="ppca64" ;;
    *) echo "ERROR: unsupported architecture '$MACH' for automatic FPC install" >&2; exit 1 ;;
  esac

  TARBALL="fpc-${FPC_VERSION}.${FPC_ARCH}-linux.tar"
  URL="https://sourceforge.net/projects/freepascal/files/Linux/${FPC_VERSION}/${TARBALL}/download"
  echo "No fpc on PATH; downloading Free Pascal ${FPC_VERSION} for ${FPC_ARCH}-linux ..."
  fetch "$URL" "$WORK/$TARBALL"
  tar -xf "$WORK/$TARBALL" -C "$WORK"

  SRC="$WORK/fpc-${FPC_VERSION}.${FPC_ARCH}-linux"
  BIN_TAR="$SRC/binary.${FPC_ARCH}-linux.tar"
  PREFIX="$HOME/.local/share/eos5f0j-fpc"
  FPCDIR="$PREFIX/lib/fpc/${FPC_VERSION}"
  mkdir -p "$PREFIX/etc"

  # Extract only what we need to compile: the compiler base and the config tool.
  tar -xOf "$BIN_TAR" "base.${FPC_ARCH}-linux.tar.gz"           | tar -xzf - -C "$PREFIX"
  tar -xOf "$BIN_TAR" "utils-fpcmkcfg.${FPC_ARCH}-linux.tar.gz" | tar -xzf - -C "$PREFIX"

  # The release archive ships the compiler at lib/fpc/<ver>/ppcXXX but not the
  # bin/ symlink the interactive installer would create. The fpc driver locates
  # the compiler next to itself, so link it into bin/ (otherwise: "ppcx64 can't
  # be executed, error code: 127").
  ln -sf "$FPCDIR/$PPC" "$PREFIX/bin/$PPC"

  # Generate a config so the driver can locate its RTL units (SysUtils, Math).
  "$FPCDIR/samplecfg" "$FPCDIR" "$PREFIX/etc" >/dev/null 2>&1 || true
  export PPC_CONFIG_PATH="$PREFIX/etc"
  FPC_BIN="$PREFIX/bin/fpc"
fi

echo "Using Free Pascal compiler: $FPC_BIN"

# 4. Download the checkmol source.
fetch "$CHECKMOL_URL" "$WORK/checkmol.pas"

# 5. Compile checkmol. Delphi mode (-S2) is REQUIRED by the source. Build inside
#    the scratch dir (keeps .o/.ppu artifacts out of the model) then move the
#    executable into the model code directory.
( cd "$WORK" && "$FPC_BIN" -S2 -O2 checkmol.pas )
mv "$WORK/checkmol" "$OUT"
chmod +x "$OUT"
echo "checkmol compiled -> $OUT"
