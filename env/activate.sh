#!/usr/bin/env bash
set -euo pipefail

# BASH_SOURCE doesn't exist when this is sourced from zsh (macOS's default
# login shell) - fall back to zsh's own "%N" script-path expansion so
# `source env/activate.sh` works from either shell.
if [[ -n "${ZSH_VERSION:-}" ]]; then
    # shellcheck disable=SC2296,SC2298
    this_script="${(%):-%N}"
else
    this_script="${BASH_SOURCE[0]}"
fi
repo_root="$(cd "$(dirname "${this_script}")/.." && pwd)"
venv="${repo_root}/.venv"
nxf_dir="${repo_root}/env/bin"

cd "${repo_root}"

if [[ ! -x "${venv}/bin/python" || ! -x "${nxf_dir}/nextflow" ]]; then
    echo "[error] drugmr environment not found at ${venv} or ${nxf_dir}." >&2
    echo "[error] run once: ./env/bootstrap.sh" >&2
    return 1 2>/dev/null || exit 1
fi

if ! "${venv}/bin/python" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)'; then
    echo "[error] ${venv} does not use python 3.12 or newer." >&2
    echo "[error] rename it and rerun: mv ${venv} ${venv}.old" >&2
    return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC1091
source "${venv}/bin/activate"

export DRUGMR_ROOT="${repo_root}"
export PATH="${nxf_dir}:${PATH}"

nxf_version_line="$(grep 'nextflowVersion' "${repo_root}/nextflow.config")"
export NXF_VER="$(grep -oE '[0-9]+\.[0-9]+\.[0-9]+' <<< "${nxf_version_line}")"

# Nextflow needs java on PATH, and unlike the venv/nextflow launcher itself,
# a `module load` from bootstrap.sh doesn't carry over into a new login
# shell - so redo it here too, every session, not just once at bootstrap time.
if ! command -v java >/dev/null 2>&1 && type module >/dev/null 2>&1; then
    java_module="$(module -t avail 2>&1 | awk -F/ '/(^|\/)Java\// && $2+0 >= 17 {print; exit}')"
    [[ -n "${java_module}" ]] && module load "${java_module}" >/dev/null 2>&1
fi

if ! command -v java >/dev/null 2>&1; then
    echo "[error] java is required to run nextflow but was not found on PATH." >&2
    return 1 2>/dev/null || exit 1
fi

if ! python -c 'import drugmr' >/dev/null 2>&1; then
    echo "[error] drugmr is not installed in ${venv}." >&2
    echo "[error] run once: ./env/bootstrap.sh" >&2
    return 1 2>/dev/null || exit 1
fi

python -c 'import drugmr; print("[done] drugmr environment active")'
echo "[done] nextflow ${NXF_VER} on PATH"
