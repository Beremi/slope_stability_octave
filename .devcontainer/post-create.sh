#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

required_paths=(
  ".octave_all/env.sh"
  ".octave_all/bin/octave-rsb"
  ".venv/bin/python"
  ".venv/bin/activate"
  "setup/activate_optimized_octave.sh"
  "slope_stability/+LINEAR_SOLVERS/hypre_boomeramg_mex.mex"
  "slope_stability/+CONSTITUTIVE_PROBLEM/constitutive_problem_3D_SDS_mex.mex"
  "slope_stability/+ASSEMBLY/assemble_K_tangent_vals.mex"
)

needs_clean=0

relocatable_path_matches() {
  local file_path="$1"
  local expected="$2"
  [[ -f "${file_path}" ]] && grep -Fq "${expected}" "${file_path}"
}

if [[ -e "${ROOT_DIR}/.octave_all/env.sh" ]] && ! relocatable_path_matches "${ROOT_DIR}/.octave_all/env.sh" "${ROOT_DIR}/.octave_all/install/"; then
  needs_clean=1
fi

if [[ -e "${ROOT_DIR}/.octave_all/bin/octave-rsb" ]] && ! relocatable_path_matches "${ROOT_DIR}/.octave_all/bin/octave-rsb" "${ROOT_DIR}/.octave_all/install/"; then
  needs_clean=1
fi

if [[ -e "${ROOT_DIR}/.venv/bin/activate" ]] && ! relocatable_path_matches "${ROOT_DIR}/.venv/bin/activate" "${ROOT_DIR}/.venv"; then
  needs_clean=1
fi

needs_bootstrap=0
for rel_path in "${required_paths[@]}"; do
  if [[ ! -e "${ROOT_DIR}/${rel_path}" ]]; then
    needs_bootstrap=1
    break
  fi
done

if [[ "${needs_clean}" == "1" ]]; then
  cd "${ROOT_DIR}"
  setup/clean_local_builds.sh
  needs_bootstrap=1
fi

if [[ "${needs_bootstrap}" == "1" ]]; then
  cd "${ROOT_DIR}"
  ./bootstrap_all.sh --no-clean --no-verify
  exit 0
fi

echo "Local Octave stack already present and path-correct; skipping bootstrap."
