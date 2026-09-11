#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv="${repo_root}/.venv"

cd "${repo_root}"

if [[ ! -x "${venv}/bin/python" ]]; then
    python_bin=""

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
