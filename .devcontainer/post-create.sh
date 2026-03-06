#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

required_paths=(
  ".octave_all/env.sh"
  ".octave_all/bin/octave-rsb"
  ".octave_all/bin/octave-rsb-jupyter"
  ".venv/bin/python"
  ".venv/bin/activate"
  ".venv/share/jupyter/kernels/octave-local-rsb/kernel.json"
  "setup/activate_optimized_octave.sh"
  "slope_stability/+LINEAR_SOLVERS/hypre_boomeramg_mex.mex"
  "slope_stability/+CONSTITUTIVE_PROBLEM/constitutive_problem_3D_SDS_mex.mex"
  "slope_stability/+ASSEMBLY/assemble_K_tangent_vals.mex"
)

needs_clean=0
needs_bootstrap=0

relocatable_path_matches() {
  local file_path="$1"
  local expected="$2"
  [[ -f "${file_path}" ]] && grep -Fq "${expected}" "${file_path}"
}

find_octave_config_header() {
  local build_root="${ROOT_DIR}/.octave_all/build"
  if [[ ! -d "${build_root}" ]]; then
    return 0
  fi

  find "${build_root}" -path '*/octave-*-build/config.h' | sort | head -n 1
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

if [[ -e "${ROOT_DIR}/setup/activate_optimized_octave.sh" ]] && ! relocatable_path_matches "${ROOT_DIR}/setup/activate_optimized_octave.sh" ".octave_all/env.sh"; then
  needs_bootstrap=1
fi

if [[ "${REQUIRE_QT_TOOLKIT:-0}" == "1" ]]; then
  octave_config_h="$(find_octave_config_header)"
  if [[ -z "${octave_config_h}" ]]; then
    needs_clean=1
  elif ! grep -q 'HAVE_QT 1' "${octave_config_h}" || ! grep -q 'HAVE_OPENGL 1' "${octave_config_h}"; then
    needs_clean=1
  fi
fi

if [[ -e "${ROOT_DIR}/.octave_all/env.sh" ]] && grep -Fq 'OMP_NUM_THREADS:-16' "${ROOT_DIR}/.octave_all/env.sh"; then
  needs_bootstrap=1
fi

if [[ -e "${ROOT_DIR}/.venv/share/jupyter/kernels/octave-local-rsb/kernel.json" ]] && grep -Fq '"OMP_NUM_THREADS"' "${ROOT_DIR}/.venv/share/jupyter/kernels/octave-local-rsb/kernel.json"; then
  needs_bootstrap=1
fi

octave_prefix=""
if [[ -d "${ROOT_DIR}/.octave_all/install" ]]; then
  octave_prefix="$(find "${ROOT_DIR}/.octave_all/install" -maxdepth 1 -mindepth 1 -type d -name 'octave-*' | sort | head -n 1)"
fi
if [[ -x "${octave_prefix}/bin/octave" ]]; then
  if ! (
    cd "${ROOT_DIR}" && \
    bash -lc 'source .octave_all/env.sh && octave --quiet --eval "assert (exist (\"delaunay\", \"file\") == 2);"' \
  ) >/dev/null 2>&1; then
    needs_bootstrap=1
  fi
fi

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
