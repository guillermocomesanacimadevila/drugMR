#!/usr/bin/env python3
import synapseclient
import polars as pl 
import os
import argparse
import subprocess
import requests

# syn51364943 for UKB-PPP
synapse_id = "syn51364943"
syn = synapseclient.Synapse()
syn.login() # pull from base/syn...
ancestry = "European (discovery)" # arg
slurm_id = ""

def download_into_hpc():

    """
    Steps:
    1. Log into SLURM HPC (make sure synapse client file in base dir & synapseclient package update within GRCH Docker container)
    2. Ensure .sif container within ./drugMR == ok
    3. Re-make the loop above but with the ukbppp_dl package
    4. Run a re-adapted version of scripts/ukb_ppp/ukb_ppp.py for 1 single file and that rm -rf the previous one 
    5. Map rsIDs
    """

    def login(slurm_id: str):
        login_cmd = f"""
set -euo pipefail
ssh {slurm_id}.falconlogin.cf.ac.uk
        """

        # first login 
        subprocess.run(login_cmd, shell=True, check=True, executable="/bin/bash")

        # check and load up .sif
        remote = f"/shared/home1/{slurm_id}/drugMR"
        sif = f"{remote}/env/drugmr.sif"

        # cd "$HOME/drugMR"
        # chmod +x bin/bootstrap_hpc.sh
        # bash bin/bootstrap_hpc.sh

        cmd_singularity = """
set -euo pipefail
cd "$HOME/drugMR"
chmod +x bin/bootstrap_hpc.sh
bash bin/bootstrap_hpc.sh # this loads all the singularity stuff
"""

    # Establish out_dir
    out_dir = "dat/pQTL/ukb_ppp/" # all .parquet files here 
    os.makedirs(out_dir, exist_ok=True)

    # download files
    for item in syn.getChildren(synapse_id):
        if item["name"] == "UKB-PPP pGWAS summary statistics":
            print(f"Entering: {item['name']}\n")
            for folder in syn.getChildren(item["id"]):
                if folder["name"] == ancestry:
                    print(f"Entering: {folder['name']}\n")
                    rows = []
                    for file in syn.getChildren(folder["id"]):
                        downloaded = syn.get(file["id"], downloadLocation=out_dir)
                        downloaded_path = downloaded.path
                        print(f"Downloaded: {downloaded_path}")

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        #  placeholder for CI/CD testing
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        txt_out = os.path.join(out_dir, file["name"].replace(".tar", "_first_line.txt"))
                        with open(downloaded_path, "rb") as f:
                            first_line = f.readline()
                        with open(txt_out, "wb") as f:
                            f.write(first_line) 

                        print(f"Saved first line: {txt_out}")

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        os.remove(downloaded_path)
                        print(f"Deleted: {downloaded_path}")
                       # continue 
                    break
            break

# Rebuild .sif file 
if __name__ == "__main__":
    download_into_hpc()