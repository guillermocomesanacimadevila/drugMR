#!/usr/bin/env bash
# No "-u": Lmod's own `module` bash function isn't nounset-safe internally,
# and this script calls it - "-u" here would risk killing the whole (sourced)
# shell on an unrelated unset variable deep inside Lmod's implementation.
set -eo pipefail

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

NXF_VER="$(grep 'nextflowVersion' "${repo_root}/nextflow.config" 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
if [[ -z "${NXF_VER}" ]]; then
    echo "[error] could not find a nextflowVersion in ${repo_root}/nextflow.config." >&2
    return 1 2>/dev/null || exit 1
fi
export NXF_VER

# Nextflow needs java on PATH, and unlike the venv/nextflow launcher itself,
# a `module load` from bootstrap.sh doesn't carry over into a new login
# shell - so redo it here too, every session, not just once at bootstrap time.
if ! command -v java >/dev/null 2>&1 && type module >/dev/null 2>&1; then
    # Lmod's terse listing renders an aliased module as e.g.
    # "Java/17(@Java/17.0.15)" on one line - strip the "(@...)" part so we
    # pass module load a real name, not that whole alias annotation.
    java_module="$(module -t avail 2>&1 | awk -F/ '/(^|\/)Java\// && $2+0 >= 17 {sub(/\(.*/, ""); print; exit}')"
    if [[ -n "${java_module}" ]]; then
        module load "${java_module}" >/dev/null 2>&1 || true
    fi
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
