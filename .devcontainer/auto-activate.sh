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
    if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
      _slope_threads=""
      _slope_quota=""
      _slope_period=""
      _slope_quota_threads=""
      if command -v getconf >/dev/null 2>&1; then
        _slope_threads="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
      fi
      if [[ -z "${_slope_threads:-}" ]] && command -v nproc >/dev/null 2>&1; then
        _slope_threads="$(nproc 2>/dev/null || true)"
      fi
      if [[ -r /sys/fs/cgroup/cpu.max ]]; then
        read -r _slope_quota _slope_period < /sys/fs/cgroup/cpu.max || true
        if [[ "${_slope_quota:-}" =~ ^[0-9]+$ && "${_slope_period:-}" =~ ^[1-9][0-9]*$ ]]; then
          _slope_quota_threads="$(( (_slope_quota + _slope_period - 1) / _slope_period ))"
        fi
      elif [[ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us && -r /sys/fs/cgroup/cpu/cpu.cfs_period_us ]]; then
        _slope_quota="$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us 2>/dev/null || true)"
        _slope_period="$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us 2>/dev/null || true)"
        if [[ "${_slope_quota:-}" =~ ^[0-9]+$ && "${_slope_period:-}" =~ ^[1-9][0-9]*$ && "${_slope_quota}" -gt 0 ]]; then
          _slope_quota_threads="$(( (_slope_quota + _slope_period - 1) / _slope_period ))"
        fi
      fi
      if [[ "${_slope_threads:-}" =~ ^[1-9][0-9]*$ && "${_slope_quota_threads:-}" =~ ^[1-9][0-9]*$ && "${_slope_quota_threads}" -lt "${_slope_threads}" ]]; then
        _slope_threads="${_slope_quota_threads}"
      elif [[ ! "${_slope_threads:-}" =~ ^[1-9][0-9]*$ && "${_slope_quota_threads:-}" =~ ^[1-9][0-9]*$ ]]; then
        _slope_threads="${_slope_quota_threads}"
      fi
      export OMP_NUM_THREADS="${_slope_threads:-1}"
      unset _slope_threads _slope_quota _slope_period _slope_quota_threads
    fi
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
