#!/bin/bash
# Build the element-level tangent assembly mex file for Octave.
# Run from the slope_stability/ directory.
#
# Usage:
#   bash +ASSEMBLY/mex/build_assemble_K_tangent.sh
#
# The compiled .mex file is placed in +ASSEMBLY/ so that
# Octave can find it via the ASSEMBLY package namespace.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_3D="$SCRIPT_DIR/assemble_K_tangent_vals.c"
SRC_2D="$SCRIPT_DIR/assemble_K_tangent_vals_2D.c"
OUT_DIR="$SCRIPT_DIR/.."

echo "Compiling element-level tangent assembly mex files ..."

# Detect the right mkoctfile
if command -v mkoctfile &>/dev/null; then
    MKOCTFILE=mkoctfile
elif [ -x "${OCTAVE_BIN%/*}/mkoctfile" ]; then
    MKOCTFILE="${OCTAVE_BIN%/*}/mkoctfile"
else
    echo "ERROR: mkoctfile not found. Source activate_optimized_octave.sh first."
    exit 1
fi

compile_one () {
    local src="$1"
    local out="$2"
    echo "  -> $(basename "$out")"
    $MKOCTFILE --mex -O2 \
        -DHAVE_OPENMP \
        "$src" \
        -o "$out" \
        -fopenmp -lgomp
}

compile_one "$SRC_3D" "$OUT_DIR/assemble_K_tangent_vals.mex"
compile_one "$SRC_2D" "$OUT_DIR/assemble_K_tangent_vals_2D.mex"

echo "Done.  Outputs:"
echo "  $OUT_DIR/assemble_K_tangent_vals.mex"
echo "  $OUT_DIR/assemble_K_tangent_vals_2D.mex"
