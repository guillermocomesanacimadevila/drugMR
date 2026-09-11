#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv="${repo_root}/.venv"

cd "${repo_root}"

# Make the module command available in non-interactive shells when the cluster
# provides the standard Environment Modules initialisation files.
if ! type module >/dev/null 2>&1; then
    for modules_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r "${modules_init}" ]]; then
            # shellcheck disable=SC1090
            source "${modules_init}"
            break
        fi
    done
fi

if [[ ! -x "${venv}/bin/python" ]]; then
    python_bin=""

    # On module-based HPC systems, Python 3.12 may exist but not be loaded.
    # Respect an explicit module name first, then discover a suitable module.
    if type module >/dev/null 2>&1; then
        module_candidates=()
        if [[ -n "${DRUGMR_PYTHON_MODULE:-}" ]]; then
            module_candidates+=("${DRUGMR_PYTHON_MODULE}")
        else
            while IFS= read -r module_name; do
                [[ -n "${module_name}" ]] && module_candidates+=("${module_name}")
            done < <(module -t avail 2>&1 | awk '/(^|\/)Python\/3\.(1[2-9]|[2-9][0-9])/{print $1}')
        fi

        for module_name in "${module_candidates[@]}"; do
            module load "${module_name}" >/dev/null 2>&1 || continue
            break
        done
    fi

    for candidate in python3.12 python3.13 python3.14 python3 python; do
        if command -v "${candidate}" >/dev/null 2>&1 &&
           "${candidate}" -c \
             'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)'
        then
            python_bin="${candidate}"
            break
        fi
    done

    if [[ -z "${python_bin}" ]]; then
        echo "[error] python 3.12 or newer is required." >&2
        exit 1
    fi

    echo "[tracking] creating ${venv} with ${python_bin}..."
    "${python_bin}" -m venv "${venv}"
fi

if ! "${venv}/bin/python" -c \
    'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)'
then
    echo "[error] existing ${venv} does not use python 3.12 or newer." >&2
    echo "[error] rename it before rerunning: mv ${venv} ${venv}.old" >&2
    exit 1
fi

echo "[tracking] installing drugmr into ${venv}..."
"${venv}/bin/python" -m pip install -e .

echo "[done] drugmr environment ready."
echo "activate with: source ${venv}/bin/activate"
