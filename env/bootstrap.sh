#!/usr/bin/env bash
# No "-u": Lmod's own `module` bash function isn't nounset-safe internally,
# and this script calls it - "-u" here would risk killing the whole script
# on an unrelated unset variable deep inside Lmod's implementation.
set -eo pipefail

# BASH_SOURCE doesn't exist when this is sourced from zsh (macOS's default
# login shell) - fall back to zsh's own "%N" script-path expansion so this
# still works if someone sources it out of habit, not just `./env/bootstrap.sh`.
if [[ -n "${ZSH_VERSION:-}" ]]; then
    # shellcheck disable=SC2296,SC2298
    this_script="${(%):-%N}"
else
    this_script="${BASH_SOURCE[0]}"
fi
repo_root="$(cd "$(dirname "${this_script}")/.." && pwd)"
venv="${repo_root}/.venv"

cd "${repo_root}"

# Make the module command available in non-interactive HPC shells.
if ! type module >/dev/null 2>&1; then
    for modules_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r "${modules_init}" ]]; then
            # shellcheck disable=SC1090
            source "${modules_init}"
            break
        fi
    done
fi

python_bin=""
if type module >/dev/null 2>&1; then
    module_candidates=()
    if [[ -n "${DRUGMR_PYTHON_MODULE:-}" ]]; then
        module_candidates+=("${DRUGMR_PYTHON_MODULE}")
    else
        while IFS= read -r module_name; do
            [[ -n "${module_name}" ]] && module_candidates+=("${module_name}")
        done < <(module -t avail 2>&1 | awk '/(^|\/)Python\/3\.(1[2-9]|[2-9][0-9])/{sub(/\(.*/, "", $1); print $1}')
    fi

    for module_name in "${module_candidates[@]}"; do
        module load "${module_name}" >/dev/null 2>&1 || continue
        break
    done
fi

is_good_python='import sys, sysconfig; raise SystemExit(0 if sys.version_info >= (3, 12) and not sysconfig.get_config_var("Py_GIL_DISABLED") else 1)'
# ^ rejects anything below 3.12, and rejects free-threaded builds (e.g. a
# stray pyenv "3.14.0t") - most of our compiled dependencies (polars,
# numpy...) have no prebuilt wheels for those yet, so uv/pip would try to
# build them from source and fail.

# pyenv shims only resolve to whichever version is currently active (its
# global/local pointer) - if that happens to be an unrelated version with no
# "python3.12" of its own, the shim fails even though 3.12 is installed. Ask
# pyenv directly for its installed versions first, bypassing the shims.
if command -v pyenv >/dev/null 2>&1; then
    while IFS= read -r pyenv_version; do
        candidate="$(pyenv root)/versions/${pyenv_version}/bin/python3"
        if [[ -x "${candidate}" ]] && "${candidate}" -c "${is_good_python}"; then
            python_bin="${candidate}"
            break
        fi
    done < <(pyenv versions --bare 2>/dev/null | sort -rV)
fi

if [[ -z "${python_bin}" ]]; then
    for candidate in python3.12 python3.13 python3.14 python3 python; do
        if command -v "${candidate}" >/dev/null 2>&1 && "${candidate}" -c "${is_good_python}"; then
            python_bin="${candidate}"
            break
        fi
    done
fi

if [[ -z "${python_bin}" ]]; then
    echo "[error] python 3.12 or newer is required." >&2
    exit 1
fi

if [[ -e "${venv}" && ! -x "${venv}/bin/python" ]]; then
    echo "[error] ${venv} exists but is not a valid virtual environment." >&2
    exit 1
fi

if [[ ! -x "${venv}/bin/python" ]]; then
    echo "[tracking] creating ${venv} with ${python_bin}..."
    "${python_bin}" -m venv "${venv}"
fi

if ! "${venv}/bin/python" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)'; then
    echo "[error] existing ${venv} does not use python 3.12 or newer." >&2
    echo "[error] rename it before rerunning: mv ${venv} ${venv}.old" >&2
    exit 1
fi

if command -v uv >/dev/null 2>&1; then
    echo "[tracking] installing locked dependencies with uv..."
    uv sync --project "${repo_root}" --active
else
    echo "[tracking] uv not found; installing with pip..."
    "${venv}/bin/python" -m pip install --progress-bar on -e .
fi

"${venv}/bin/python" -c 'import drugmr; print("[done] drugmr environment ready")'

# Nextflow needs a matching version too - self-manage it here rather than
# relying on whatever "module load Nextflow" happens to offer on a given
# cluster. nextflow.config's manifest is the single source of truth for
# which version.
if ! command -v java >/dev/null 2>&1 && type module >/dev/null 2>&1; then
    # Nextflow needs Java 17+ - take the lowest version that still clears
    # that bar, not just whatever module lists first (often a stale Java 11)
    # and not the cluster's newest default either (often too new to trust).
    # Lmod's terse listing renders an aliased module as e.g.
    # "Java/17(@Java/17.0.15)" on one line - strip the "(@...)" part so we
    # pass module load a real name, not that whole alias annotation.
    java_module="$(module -t avail 2>&1 | awk -F/ '/(^|\/)Java\// && $2+0 >= 17 {sub(/\(.*/, ""); print; exit}')"
    if [[ -n "${java_module}" ]]; then
        module load "${java_module}" >/dev/null 2>&1 || true
    fi
fi

if ! command -v java >/dev/null 2>&1; then
    echo "[error] java is required to run nextflow." >&2
    exit 1
fi

nxf_version="$(grep 'nextflowVersion' "${repo_root}/nextflow.config" 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
if [[ -z "${nxf_version}" ]]; then
    echo "[error] could not find a nextflowVersion in ${repo_root}/nextflow.config." >&2
    exit 1
fi
nxf_dir="${repo_root}/env/bin"

if [[ ! -x "${nxf_dir}/nextflow" ]]; then
    echo "[tracking] installing the nextflow launcher..."
    mkdir -p "${nxf_dir}"
    (cd "${nxf_dir}" && curl -fsSL https://get.nextflow.io | bash) >/dev/null
fi

echo "[tracking] fetching nextflow ${nxf_version}..."
NXF_VER="${nxf_version}" "${nxf_dir}/nextflow" -version >/dev/null
echo "[done] nextflow ${nxf_version} ready"

if [[ "${this_script}" != "${0}" ]]; then
    # Preserve the historical behaviour when this script is sourced.
    # shellcheck disable=SC1091
    source "${venv}/bin/activate"
    export DRUGMR_ROOT="${repo_root}"
    echo "[done] activated ${venv}"
else
    echo "[next] run: source env/activate.sh"
fi
