# shellcheck shell=bash

if [[ $- != *i* ]]; then
  return 0
fi

if [[ "${SLOPE_STABILITY_ENV_ACTIVATED:-0}" == "1" ]]; then
  return 0
fi

repo_root="${SLOPE_STABILITY_REPO:-}"
if [[ -z "${repo_root}" || ! -f "${repo_root}/bootstrap_all.sh" || ! -d "${repo_root}/setup" ]]; then
  detected_root="$(git -C "${PWD}" rev-parse --show-toplevel 2>/dev/null || true)"
  if [[ -n "${detected_root}" && -f "${detected_root}/bootstrap_all.sh" && -d "${detected_root}/setup" ]]; then
    repo_root="${detected_root}"
  else
    return 0
  fi
fi

activated=0

install_root="${repo_root}/.octave_all/install"
openblas_prefix=""
octave_prefix=""
librsb_prefix=""
octave_libexec_dir=""

if [[ -d "${install_root}" ]]; then
  openblas_prefix="$(find "${install_root}" -maxdepth 1 -mindepth 1 -type d -name 'openblas-*' | sort | head -n 1)"
  octave_prefix="$(find "${install_root}" -maxdepth 1 -mindepth 1 -type d -name 'octave-*' | sort | head -n 1)"
  librsb_prefix="$(find "${install_root}" -maxdepth 1 -mindepth 1 -type d -name 'librsb-*' | sort | head -n 1)"
fi

if [[ -n "${octave_prefix}" && -d "${octave_prefix}/lib/octave" ]]; then
  octave_libexec_dir="$(find "${octave_prefix}/lib/octave" -maxdepth 1 -mindepth 1 -type d | sort | head -n 1)"
fi

if [[ -n "${openblas_prefix}" && -n "${octave_prefix}" && -n "${librsb_prefix}" && -n "${octave_libexec_dir}" ]]; then
  octave_lib_path="${octave_libexec_dir}:${octave_prefix}/lib:${librsb_prefix}/lib:${openblas_prefix}/lib"
  if env LD_LIBRARY_PATH="${octave_lib_path}:${LD_LIBRARY_PATH:-}" "${octave_prefix}/bin/octave-cli" --version >/dev/null 2>&1; then
    export OPENBLAS_PREFIX="${openblas_prefix}"
    export OCTAVE_PREFIX="${octave_prefix}"
    export LIBRSB_PREFIX="${librsb_prefix}"
    export OCTAVE_BIN="${octave_prefix}/bin/octave"
    export PATH="${octave_prefix}/bin:${librsb_prefix}/bin:${PATH}"
    export LD_LIBRARY_PATH="${octave_lib_path}:${LD_LIBRARY_PATH:-}"
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
    activated=1
  fi
fi

if [[ -x "${repo_root}/.venv/bin/python" ]] && [[ "${VIRTUAL_ENV:-}" != "${repo_root}/.venv" ]] && "${repo_root}/.venv/bin/python" -V >/dev/null 2>&1; then
  export VIRTUAL_ENV="${repo_root}/.venv"
  export PATH="${repo_root}/.venv/bin:${PATH}"
  unset PYTHONHOME
  activated=1
fi

if [[ "${activated}" == "1" ]]; then
  export SLOPE_STABILITY_ENV_ACTIVATED=1
fi
