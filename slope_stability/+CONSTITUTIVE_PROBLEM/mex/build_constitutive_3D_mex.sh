#!/bin/bash
# Build the constitutive_problem 2D/3D mex files for Octave.
# Run from the slope_stability/ directory.
#
# Usage:
#   bash +CONSTITUTIVE_PROBLEM/mex/build_constitutive_3D_mex.sh
#
# Produces four .mex files in +CONSTITUTIVE_PROBLEM/:
#   constitutive_problem_2D_S_mex.mex     (stress only)
#   constitutive_problem_2D_SDS_mex.mex   (stress + consistent tangent)
#   constitutive_problem_3D_S_mex.mex     (stress only)
#   constitutive_problem_3D_SDS_mex.mex   (stress + consistent tangent)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$SCRIPT_DIR/.."

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
    echo "Compiling $(basename "$src") ..."
    $MKOCTFILE --mex -O2 \
        -DHAVE_OPENMP \
        -I"$SCRIPT_DIR" \
        "$src" \
        -o "$out" \
        -fopenmp -lgomp
}

compile_one "$SCRIPT_DIR/constitutive_problem_2D_S_mex.c" "$OUT_DIR/constitutive_problem_2D_S_mex.mex"
compile_one "$SCRIPT_DIR/constitutive_problem_2D_SDS_mex.c" "$OUT_DIR/constitutive_problem_2D_SDS_mex.mex"
compile_one "$SCRIPT_DIR/constitutive_problem_3D_S_mex.c" "$OUT_DIR/constitutive_problem_3D_S_mex.mex"
compile_one "$SCRIPT_DIR/constitutive_problem_3D_SDS_mex.c" "$OUT_DIR/constitutive_problem_3D_SDS_mex.mex"

echo "Done.  Outputs:"
echo "  $OUT_DIR/constitutive_problem_2D_S_mex.mex"
echo "  $OUT_DIR/constitutive_problem_2D_SDS_mex.mex"
echo "  $OUT_DIR/constitutive_problem_3D_S_mex.mex"
echo "  $OUT_DIR/constitutive_problem_3D_SDS_mex.mex"
