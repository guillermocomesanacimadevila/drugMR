#!/usr/bin/env bash
set -euo pipefail

# -> this script is to build .venv/ after git clone ...
# to avoid use of conda when re-using pipeline

# check whether .venv/ is there already
if [ ! -d .venv ]; then
    echo "[launch.sh] creating virtual environment (.venv/)..."
    python3.12 -m venv .venv
fi

source .venv/bin/activate

echo "[launch.sh] installing drugmr + notebook deps..."
pip install -q --upgrade pip 
pip install -q -e ".[notebook]"

echo "[launch.sh] launching Jupyter Lab..."
exec jupyter lab notebooks/00_drugmr.ipynb