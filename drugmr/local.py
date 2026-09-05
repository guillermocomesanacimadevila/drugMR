#!/usr/bin/env python3
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from drugmr import paths, registry
from drugmr.config import Config

# LOCAL RESULTS / DASHBOARD STUFF
# Need to create local running functions including a pulling docker from container function
# So then still the QC+MR runs in micromamba Docker env
# PostgreSQL db pulling and dashboard == jupyter (with .toml in ./)
# Goal == have a flagging variable (local/hpc)

def cmd_base(cmd: list):
    """
    Baseline function for parsing and running CLI-based .py scripts.
    cmd is an argv list (no shell) - interpolated values (image names,
    dataset/config values from params/*.yaml) are passed as literal argv
    elements rather than parsed as shell syntax.
    """
    return subprocess.run(cmd, check=True)

def check_output(path: Path, step: str, overwrite: bool = False):
    # run step if overwrite == True
    if overwrite:
        print(f"[TRACKING] Overwrite enabled - rerunning {step}...")
        return False

    # run step if output does not exist
    if not path.exists():
        print(f"[TRACKING] No existing {step} output found - running step...")
        return False

    # run step if output exists but is empty
    if path.stat().st_size == 0:
        print(f"[CONCERN] {step} output exists but is empty - rerunning step...")
        return False

    print(f"[TRACKING] {step} already completed: {path}")
    print(f"[TRACKING] Skipping {step}...")
    return True

def check_cis_regions(cis_dir: Path, overwrite: bool = False):
    # run step if overwrite == True
    if overwrite:
        print("[TRACKING] Overwrite enabled - rerunning cis-region preparation...")
        return False

    # run step if cis-region directory does not exist
    if not cis_dir.exists():
        print("[TRACKING] No existing cis-region directory found - running step...")
        return False

    # check whether complete pQTL + GWAS cis-regions actually exist
    pqtl_loci = {file.parent.name for file in cis_dir.glob("*/pqtl.parquet")}
    gwas_loci = {file.parent.name for file in cis_dir.glob("*/gwas.parquet")}
    complete_loci = pqtl_loci.intersection(gwas_loci)
    incomplete_loci = pqtl_loci.symmetric_difference(gwas_loci)

    if len(complete_loci) == 0:
        print("[CONCERN] cis-region directory exists but no complete pQTL/GWAS loci were found - rerunning step...")
        return False

    if len(incomplete_loci) > 0:
        print(f"[CONCERN] Found {len(incomplete_loci)} incomplete cis-region loci - rerunning step...")
        print(f"[CONCERN] Example incomplete loci: {sorted(incomplete_loci)[:10]}")
        return False

    print(f"[TRACKING] cis-regions already completed: {len(complete_loci)} complete loci found")
    print("[TRACKING] Skipping cis-region preparation...")
    return True

def require_output(path: Path, step: str, required_for: str):
    # do not run downstream step where required upstream output does not exist
    if not path.exists():
        raise FileNotFoundError(
            f"{required_for} cannot run because {step} output was not found: {path}"
        )

    # do not run downstream step where required upstream output is empty
    if path.stat().st_size == 0:
        raise RuntimeError(
            f"{required_for} cannot run because {step} output is empty: {path}"
        )

def results(
    config: str,
    db_id: str = "drugmr",
    dashboard_script: str = "dashboard/mr_app.py",
    db_script: str = "bin/load_db_into_postgres.py",
    write_secrets_script: str = "bin/write_streamlit_secrets.py",
    port_number: int = 5433,
):
    project_root = Path(__file__).resolve().parents[1]
    cfg = Config(project_root / config)
    pqtl_dataset = cfg.pqtl_dataset
    pheno_id = cfg.pheno_id

    run_id = registry.get_latest_run_id(pheno_id, pqtl_dataset, root=str(project_root / "runs"))
    if run_id is None:
        raise FileNotFoundError(
            f"No recorded run found for pheno_id={pheno_id!r}, pqtl_dataset={pqtl_dataset!r}. "
            "Run dm.local(config=...) (or dm.hpc(...)) first."
        )
    out_dir = str(paths.run_results_dir(run_id, root=str(project_root / "runs")))

    mr_res = project_root / paths.mr_out(pqtl_dataset, pheno_id, out_dir)
    coloc_res = project_root / paths.coloc_out(pqtl_dataset, pheno_id, out_dir)
    db_script = project_root / db_script
    dashboard_script = project_root / dashboard_script
    write_secrets_script = project_root / write_secrets_script

    require_output(mr_res, "cis-MR", "PostgreSQL loading")
    require_output(coloc_res, "COLOC", "PostgreSQL loading")

    print("[TRACKING] Starting PostgreSQL (docker compose)...")
    subprocess.run(
        ["docker", "compose", "up", "-d"],
        cwd=str(project_root),
        check=True,
    )

    print("[TRACKING] Writing Streamlit secrets...")
    subprocess.run(
        [sys.executable, str(write_secrets_script)],
        check=True,
    )

    print("[TRACKING] Loading MR results into PostgreSQL...")

    subprocess.run(
        [
            sys.executable,
            str(db_script),
            "--results_file",
            str(mr_res),
            "--db_id",
            db_id,
            "--run_id",
            run_id,
            "--pqtl_dataset",
            pqtl_dataset,
            "--table",
            "cis_mr_results",
        ],
        check=True,
    )

    print("[TRACKING] Loading COLOC results into PostgreSQL...")

    subprocess.run(
        [
            sys.executable,
            str(db_script),
            "--results_file",
            str(coloc_res),
            "--db_id",
            db_id,
            "--run_id",
            run_id,
            "--pqtl_dataset",
            pqtl_dataset,
            "--table",
            "coloc_results",
        ],
        check=True,
    )

    print("[TRACKING] Launching Streamlit dashboard...")

    # Streamlit resolves .streamlit/config.toml (the custom theme - primaryColor,
    # background/text colors, font) relative to its OWN process working
    # directory, not this script's location. Without an explicit cwd here, this
    # subprocess inherits whatever cwd the CALLER (e.g. a notebook in notebooks/,
    # which has no config.toml of its own) happened to have, so the theme
    # silently fell back to Streamlit defaults instead of the one actually
    # defined at <project_root>/.streamlit/config.toml. Pinning cwd=project_root
    # fixes that, and also makes the relative dashboard_script path resolve
    # correctly regardless of the caller's own cwd.
    subprocess.run(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(project_root / dashboard_script),
            "--",
            "--db_name",
            db_id,
            "--port_number",
            str(port_number),
            "--phenotype",
            pheno_id,
            "--pqtl_dataset",
            pqtl_dataset
        ],
        cwd=str(project_root),
        check=True,
    )

# secondment functions
# if run with local -> load up docker container
# i.e. check whether the container exists within current local machine 
# if so -> load docker container -> and run within that
# produce the same output as in the cloud -> which then run local scripts to load into postgres db and then dashboard
# to check docker container - create a function within local()

def local(config: str, run_id: str = None):
    # config has no default on purpose - there's no single correct params file
    # anymore now that each (pheno_id, pqtl_dataset) pair has its own under
    # params/ (e.g. params/AD.wingo_brain.yaml) - pass one explicitly.
    #
    # run_id defaults to None, which keeps the deterministic
    # (pheno_id, pqtl_dataset, day, git commit) behaviour below. Pass an
    # existing runs/<run_id> value explicitly to resume/retry into that same
    # run dir instead - e.g. after fixing a bug in one step, to pick up
    # where a previous day's run left off rather than starting a fresh
    # runs/<run_id> dir (and therefore rerunning every step from scratch).
    project_root = Path(__file__).resolve().parents[1]
    cfg = Config(project_root / config)
    pheno_id = cfg.pheno_id
    sumstats = cfg.sumstats
    n_cases = cfg.n_cases
    n_controls = cfg.n_controls
    pqtl_dataset = cfg.pqtl_dataset
    pqtl_dir = cfg.pqtl_dir
    ref_bfile = cfg.ref_bfile
    snp_col = cfg.snp_col
    a1_col = cfg.a1_col
    a2_col = cfg.a2_col
    beta_col = cfg.beta_col
    se_col = cfg.se_col
    p_col = cfg.p_col
    pos_col = cfg.pos_col
    chr_col = cfg.chr_col
    af_col = cfg.af_col
    genome_build = cfg.genome_build
    target_build = cfg.target_build
    maf = getattr(cfg, "maf", 0.01)
    info_threshold = getattr(cfg, "info_threshold", None)
    info_col = getattr(cfg, "info_col", None)
    remove_mhc = getattr(cfg, "remove_mhc", True)
    remove_apoe = getattr(cfg, "remove_apoe", False)
    overwrite = getattr(cfg, "overwrite", False)
    image_uri = getattr(cfg, "image_uri", "ghcr.io/guillermocomesanacimadevila/drugmr:latest")
    image_name = getattr(cfg, "image_name", "ghcr.io/guillermocomesanacimadevila/drugmr:latest")
    run_smr = getattr(cfg, "run_smr", True)
    bulk_eqtl_datasets = getattr(cfg, "bulk_eqtl_datasets", [])
    sc_eqtl_dataset = getattr(cfg, "sc_eqtl_dataset", "")

    # cis-MR / coloc gate thresholds - see params/schema.json's gates block;
    # defaults match what bin/coloc_targets.py used to hardcode
    wald_fdr_q = cfg.gate("cis_mr", "wald_fdr_q", 0.05)
    ivw_fdr_q = cfg.gate("cis_mr", "ivw_fdr_q", 0.05)
    cochran_q_pval = cfg.gate("cis_mr", "cochran_q_pval", 0.05)
    egger_intercept_pval_min = cfg.gate("cis_mr", "egger_intercept_pval_min", 0)
    min_instruments_for_ivw = cfg.gate("cis_mr", "min_instruments_for_ivw", 3)
    pp4_threshold = cfg.gate("coloc", "pp4_threshold", 0.7)


    # set projectDir()
    project_root = Path(__file__).resolve().parents[1] # i.e. "Users/.../drugMR"

    # run_id is deterministic for a given (pheno_id, pqtl_dataset, day, git commit) -
    # rerunning today on the same commit reuses the same runs/<run_id>/ dir (and its
    # check_output()/require_output() skip behavior); a new day or a new commit starts
    # a fresh, separately-tracked run rather than silently overwriting the old one -
    # unless run_id is passed explicitly, in which case that existing run dir is
    # reused as-is (its own check_output()/require_output() gates then decide what
    # still needs to run)
    git_sha7 = registry.current_git_sha7(cwd=project_root)
    date_str = datetime.now().strftime("%Y%m%d")
    if run_id is None:
        run_id = paths.make_run_id(pheno_id, pqtl_dataset, date_str, git_sha7)
    elif not run_id.startswith(f"{pheno_id}_{pqtl_dataset}_"):
        raise ValueError(
            f"run_id {run_id!r} does not match config's pheno_id={pheno_id!r}, "
            f"pqtl_dataset={pqtl_dataset!r} - refusing to write into a run dir "
            "for a different (pheno_id, pqtl_dataset) pair."
        )
    out_dir = str(paths.run_results_dir(run_id))
    print(f"[TRACKING] run_id: {run_id}")

    # define all outputs first so pipeline knows what has already been ran
    qc_out = project_root / paths.qc_out(pheno_id)
    cis_dir = project_root / "dat" / "cis_regions" / pqtl_dataset
    mr_out = project_root / paths.mr_out(pqtl_dataset, pheno_id, out_dir)
    mr_instruments_out = project_root / paths.mr_instruments_out(pqtl_dataset, pheno_id, out_dir)
    coloc_out = project_root / paths.coloc_out(pqtl_dataset, pheno_id, out_dir)

    # PWCoCo (conditional coloc) - runs alongside coloc_out above, not instead of it
    # (see project_pwcoco_wiring memory); a failure here is logged as a [CONCERN] and
    # does not halt the run, since standard COLOC is the required path and PWCoCo is
    # a complementary annotation on top of it
    pwcoco_out = project_root / paths.pwcoco_out(pqtl_dataset, pheno_id, out_dir)

    # PWCoCo (eQTL-informed) - eQTL-pQTL / eQTL-GWAS on SMR-passing targets, gated
    # on the same complementary-not-required basis as pwcoco_out above
    pwcoco_qtl_out = project_root / paths.pwcoco_eqtl_pqtl_out(pqtl_dataset, pheno_id, out_dir)

    # SMR (bulk and/or single-cell) - promising target output per eQTL mode
    # bulk_eqtl_datasets is a list (MetaBrain / GTEx_v10 etc. are pre-computed
    # separately under results/SMR/bulk/{dataset}/) so its per-dataset outputs are built
    # inside the SMR step below rather than up front here
    smr_sc_out = project_root / paths.smr_sc_out(pqtl_dataset, pheno_id, sc_eqtl_dataset, out_dir)

    # final harmonised target stats
    target_stats_out = project_root / paths.target_stats_out(pqtl_dataset, pheno_id, out_dir)

    # phewas - ukb_ppp_AD_PheWAS-FinnGen.tsv
    phewas_out = project_root / paths.phewas_out(pqtl_dataset, pheno_id, out_dir)

    # phewas ukbb_out
    phewas_ukbb_out = project_root / paths.phewas_ukbb_out(pqtl_dataset, pheno_id, out_dir)

    # FinnGen PheWAS coverage manifest - UKB PheWAS reads this to run only on the
    # fallback set (targets with zero retained instruments in FinnGen)
    phewas_finngen_coverage_out = project_root / paths.phewas_finngen_coverage_out(pqtl_dataset, pheno_id, out_dir)

    def check_docker():
        # genuinely multi-statement (if/else, command availability checks) -
        # written directly in Python instead of a shell script string, so
        # image_name/image_uri (from params/*.yaml) never pass through a shell
        print("[TRACKING] Checking Docker...")
        if shutil.which("docker") is None:
            print("[ERROR] Mate, install Docker before you run this locally.")
            sys.exit(1)

        if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
            print("[ERROR] Mate, Docker is installed but it is not running.")
            print("[HINT] Open Docker Desktop and try again.")
            sys.exit(1)

        print("[TRACKING] Docker is installed and running.")

        inspect = subprocess.run(["docker", "image", "inspect", image_name], capture_output=True)
        if inspect.returncode == 0:
            print("[TRACKING] drugMR Docker image already exists locally.")
        else:
            print("[TRACKING] drugMR Docker image not found locally.")
            print("[TRACKING] Pulling drugMR image from GHCR...")
            subprocess.run(["docker", "pull", image_uri], check=True)

    # call check docker function
    check_docker()

    info_args = []
    if info_col is not None:
        info_args += ["--info-col", str(info_col)]
    if info_threshold is not None:
        info_args += ["--info-threshold", str(info_threshold)]

    flag_args = []
    if remove_mhc:
        flag_args.append("--remove_mhc")
    if remove_apoe:
        flag_args.append("--remove_apoe")

    # running individual modules
    cmd_qc = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        image_name,
        "python", "bin/qc_gwas.py",
        "--pheno-id", pheno_id,
        "--sumstats", str(sumstats),
        "--out-dir", str(paths.qc_out(pheno_id).parent),
        "--maf", str(maf),
        "--snp-col", snp_col,
        "--a1-col", a1_col,
        "--a2-col", a2_col,
        "--beta-col", beta_col,
        "--se-col", se_col,
        "--p-col", p_col,
        "--pos-col", pos_col,
        "--chr-col", chr_col,
        "--af_col", af_col,
        "--genome_build", genome_build,
        "--target_build", target_build,
        "--n_cases", str(n_cases),
        "--n_controls", str(n_controls),
        "--falcon-user", "local",
        *info_args,
        *flag_args,
    ]

    if not check_output(qc_out, "GWAS QC", overwrite):
        print("[TRACKING] Running GWAS QC locally via Docker...")
        cmd_base(cmd_qc)

    require_output(qc_out, "GWAS QC", "cis-region preparation")

    # cis-region module
    cmd_cis = [
        "docker", "run", "--rm",
        "--platform", "linux/amd64",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/prep_cis_regions.py",
        "--pqtl_dataset", pqtl_dataset,
        "--pheno_id", pheno_id,
        "--pqtl_dir", str(pqtl_dir),
    ]

    if not check_cis_regions(cis_dir, overwrite):
        print("[TRACKING] Preparing cis-regions locally...")
        cmd_base(cmd_cis)

    print(f"[TRACKING] Checking cis-region output: {cis_dir}")

    pqtl_loci = {file.parent.name for file in cis_dir.glob("*/pqtl.parquet")}
    gwas_loci = {file.parent.name for file in cis_dir.glob("*/gwas.parquet")}
    complete_loci = pqtl_loci.intersection(gwas_loci)
    incomplete_loci = pqtl_loci.symmetric_difference(gwas_loci)
    print(f"[TRACKING] Complete cis-region loci generated: {len(complete_loci)}")

    if len(complete_loci) == 0:
        raise RuntimeError("No complete cis-region files generated. Check pqtl_dir path.")

    if len(incomplete_loci) > 0:
        raise RuntimeError(
            f"{len(incomplete_loci)} incomplete cis-region loci found. "
            f"Example loci: {sorted(incomplete_loci)[:10]}"
        )

    # cis-MR module
    cmd_mr = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        image_name,
        "Rscript", "bin/cis_mr.R",
        pqtl_dataset,
        f"dat/cis_regions/{pqtl_dataset}",
        pheno_id,
        str(paths.qc_out(pheno_id)),
        str(ref_bfile),
        out_dir,
    ]

    if not check_output(mr_out, "cis-MR", overwrite):
        print("[TRACKING] Running cis-MR locally via Docker...")
        cmd_base(cmd_mr)

    require_output(mr_out, "cis-MR", "COLOC")

    # coloc target module
    cmd_coloc = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/coloc_targets.py",
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", out_dir,
        "--pqtl_dir", f"dat/cis_regions/{pqtl_dataset}",
        "--pheno_id", pheno_id,
        "--n_cases", str(n_cases),
        "--n_controls", str(n_controls),
        "--wald_fdr_q", str(wald_fdr_q),
        "--ivw_fdr_q", str(ivw_fdr_q),
        "--cochran_q_pval", str(cochran_q_pval),
        "--egger_intercept_pval_min", str(egger_intercept_pval_min),
        "--min_instruments_for_ivw", str(min_instruments_for_ivw),
    ]

    if not check_output(coloc_out, "COLOC", overwrite):
        print("[TRACKING] Running COLOC locally...")
        cmd_base(cmd_coloc)

    require_output(coloc_out, "COLOC", "Top cis-hit compilation")

    # PWCoCo (conditional coloc) - complementary to standard COLOC above, not a
    # replacement (see project_pwcoco_wiring memory): runs on the same cis-MR-passing
    # targets and its results are joined against coloc_out downstream (dashboard
    # coloc_support annotation), not used to gate anything in this orchestration.
    cmd_pwcoco = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/pwcoco_wrapper.py",
        "--pqtl_dataset", pqtl_dataset,
        "--pheno_id", pheno_id,
        "--ref_bfile", str(ref_bfile),
        "--n_cases", str(n_cases),
        "--n_controls", str(n_controls),
        "--local_results_dir", out_dir,
        "--cochran_q_pval", str(cochran_q_pval),
        "--wald_fdr_q", str(wald_fdr_q),
    ]

    if not check_output(pwcoco_out, "PWCoCo", overwrite):
        print("[TRACKING] Running PWCoCo locally...")
        try:
            cmd_base(cmd_pwcoco)
        except subprocess.CalledProcessError as error:
            print(f"[CONCERN] PWCoCo run failed - continuing without it: {error}")

    # compile final hits
    cmd_compile_top_hits = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=/work",
        image_name,
        "python", "bin/compile_cis_hit_info.py",
        "--pheno_id", pheno_id,
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", out_dir,
    ]

    if not check_output(target_stats_out, "Top cis-hit compilation", overwrite):
        print("[TRACKING] Compiling harmonised top cis-hit table...")
        cmd_base(cmd_compile_top_hits)

    require_output(
        target_stats_out,
        "Top cis-hit compilation",
        "pipeline completion"
    )

    # SMR module (bulk and/or single-cell eQTL, run right after coloc + top-cis-hit compilation)
    # -> targets which survive cis-MR + COLOC are checked against SMR + HEIDI in the
    #    configured eQTL dataset(s), alleles aligned to the AD risk allele
    if run_smr:
        if bulk_eqtl_datasets:
            # bulk eQTL SMR (MetaBrain / GTEx_v10) is pre-computed elsewhere -
            # bin/sort_smr.py ingests results/SMR/bulk/{dataset}/ rather than re-running SMR
            for bulk_dataset in bulk_eqtl_datasets:
                smr_bulk_out = project_root / paths.smr_bulk_out(pqtl_dataset, pheno_id, bulk_dataset, out_dir)

                cmd_smr_bulk = [
                    "docker", "run", "--rm",
                    "-v", f"{project_root}:/work",
                    "-w", "/work",
                    "-e", "PYTHONPATH=.",
                    image_name,
                    "python", "bin/sort_smr.py",
                    "--pheno_id", pheno_id,
                    "--sumstats", str(paths.qc_out(pheno_id)),
                    "--pqtl_dataset", pqtl_dataset,
                    "--eqtl_dataset", bulk_dataset,
                    "--eqtl_mode", "bulk",
                    "--ref_bfile", str(ref_bfile),
                    "--maf", str(maf),
                    "--local_results_dir", out_dir,
                ]

                if not check_output(smr_bulk_out, f"Bulk SMR ({bulk_dataset})", overwrite):
                    print(f"[TRACKING] Ingesting pre-computed bulk eQTL SMR for {bulk_dataset} via Docker...")
                    cmd_base(cmd_smr_bulk)
        else:
            print("[TRACKING] No bulk_eqtl_datasets specified, skipping bulk SMR.")

        if sc_eqtl_dataset:
            cmd_smr_sc = [
                "docker", "run", "--rm",
                "-v", f"{project_root}:/work",
                "-w", "/work",
                "-e", "PYTHONPATH=.",
                image_name,
                "python", "bin/sort_smr.py",
                "--pheno_id", pheno_id,
                "--sumstats", str(paths.qc_out(pheno_id)),
                "--pqtl_dataset", pqtl_dataset,
                "--eqtl_dataset", sc_eqtl_dataset,
                "--eqtl_mode", "single_cell",
                "--ref_bfile", str(ref_bfile),
                "--maf", str(maf),
                "--local_results_dir", out_dir,
            ]

            if not check_output(smr_sc_out, "Single-cell SMR", overwrite):
                print("[TRACKING] Running single-cell eQTL SMR locally via Docker...")
                cmd_base(cmd_smr_sc)
        else:
            print("[TRACKING] No sc_eqtl_dataset specified, skipping single-cell SMR.")
    else:
        print("[TRACKING] run_smr is False, skipping SMR entirely.")

    # PWCoCo (eQTL-informed) - eQTL-pQTL / eQTL-GWAS PWCoCo on every SMR-passing
    # target, then compared for shared colocalising SNPs against the pQTL-GWAS
    # PWCoCo above (see project_pwcoco_wiring memory) - runs only when SMR did,
    # since it depends on smr_final_targets_out; non-fatal like PWCoCo above.
    if run_smr:
        cmd_pwcoco_qtl = [
            "docker", "run", "--rm",
            "-v", f"{project_root}:/work",
            "-w", "/work",
            "-e", "PYTHONPATH=.",
            image_name,
            "python", "bin/pwcoco_qtl_wrapper.py",
            "--pqtl_dataset", pqtl_dataset,
            "--pheno_id", pheno_id,
            "--ref_bfile", str(ref_bfile),
            "--n_cases", str(n_cases),
            "--n_controls", str(n_controls),
            "--local_results_dir", out_dir,
        ]

        if not check_output(pwcoco_qtl_out, "PWCoCo (eQTL)", overwrite):
            print("[TRACKING] Running PWCoCo (eQTL) locally...")
            try:
                cmd_base(cmd_pwcoco_qtl)
            except subprocess.CalledProcessError as error:
                print(f"[CONCERN] PWCoCo (eQTL) run failed - continuing without it: {error}")

    # HyPrColoc module (bulk and/or single-cell eQTL) - run right after SMR so
    # the combined final multi-omics target table (bulk + single-cell) is complete.
    # For every target x cell-type/tissue hit supported by a configured eQTL dataset,
    # runs a 3-trait (pQTL / GWAS / eQTL) HyPrColoc restricted to that target's
    # cis-region, matched on shared SNPs (see drugmr.extract_common_snps). Each
    # dataset is run (and gated) independently so bulk and single-cell compose.
    hyprcoloc_out = project_root / paths.hyprcoloc_out(pqtl_dataset, pheno_id, out_dir)

    hyprcoloc_eqtl_datasets = list(bulk_eqtl_datasets) + ([sc_eqtl_dataset] if sc_eqtl_dataset else [])

    if run_smr and hyprcoloc_eqtl_datasets:
        for hc_dataset in hyprcoloc_eqtl_datasets:
            hc_dataset_out = project_root / paths.hyprcoloc_dataset_out(pqtl_dataset, hc_dataset, pheno_id, out_dir)

            cmd_hyprcoloc = [
                "docker", "run", "--rm",
                "-v", f"{project_root}:/work",
                "-w", "/work",
                "-e", "PYTHONPATH=.",
                image_name,
                "python", "bin/hyprcoloc_targets.py",
                "--pqtl_dataset", pqtl_dataset,
                "--pheno_id", pheno_id,
                "--eqtl_dataset", hc_dataset,
                "--local_results_dir", out_dir,
            ]

            if not check_output(hc_dataset_out, f"HyPrColoc ({hc_dataset})", overwrite):
                print(f"[TRACKING] Running HyPrColoc for {hc_dataset} locally via Docker...")
                cmd_base(cmd_hyprcoloc)
    else:
        print("[TRACKING] No bulk_eqtl_datasets or sc_eqtl_dataset specified (or run_smr is False), skipping HyPrColoc.")

    # PheWAS stuff (for FinnGen)
    cmd_phewas = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/phewas_cis_pqtls.py",
        "--pheno_id", pheno_id,
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", out_dir,
    ]

    # PheWAS depends on the pairwise COLOC results + cis-MR instruments
    require_output(coloc_out, "COLOC", "FinnGen PheWAS")
    require_output(mr_instruments_out, "cis-MR instruments", "FinnGen PheWAS")

    if not check_output(phewas_out, "PheWAS safety analysis on FinnGen", overwrite):
        print("[TRACKING] Running PheWAS (FinnGen) safety analysis locally...")
        cmd_base(cmd_phewas)

    require_output(phewas_out, "PheWAS safety analysis for FinnGen", "pipeline completion")
    print(f"[TRACKING] FinnGen PheWAS safety results found: {phewas_out}")


    # PheWAS (for UKBB)
    cmd_phewas_ukbb = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/work",
        "-w", "/work",
        "-e", "PYTHONPATH=.",
        image_name,
        "python", "bin/ukb_phewas.py",
        "--pheno_id", pheno_id,
        "--pqtl_dataset", pqtl_dataset,
        "--local_results_dir", out_dir,
    ]

    require_output(coloc_out, "COLOC", "UKBB PheWAS")
    require_output(mr_instruments_out, "cis-MR instruments", "UKBB PheWAS")
    # UKB PheWAS is a fallback - only run for targets with zero FinnGen instrument
    # coverage - and reads this manifest internally to build that fallback set
    require_output(phewas_finngen_coverage_out, "FinnGen PheWAS coverage manifest", "UKBB PheWAS")

    if not check_output(phewas_ukbb_out, "PheWAS safety analysis on UKBB", overwrite):
        print("[TRACKING] Running PheWAS (UKBB) safety analysis locally...")
        cmd_base(cmd_phewas_ukbb)
    require_output(phewas_ukbb_out, "PheWAS safety analysis for UKBB", "pipeline completion")
    print(f"[TRACKING] UKBB PheWAS safety results found: {phewas_ukbb_out}")

    print(f"[TRACKING] cis-MR instruments found: {mr_instruments_out}")
    print(f"[TRACKING] Final harmonised target stats found: {target_stats_out}")

    # Reached only if every check_output()/require_output() gate above passed - a
    # docker/subprocess failure earlier raises and this is never reached, so the
    # registry can never point at a partial/failed run.
    registry.write_manifest(
        run_id,
        {
            "pheno_id": pheno_id,
            "pqtl_dataset": pqtl_dataset,
            "git_sha7": git_sha7,
            "date": date_str,
            "created_at": datetime.now().isoformat(),
            "mode": "local",
            "image_name": image_name,
            "overwrite": overwrite,
        },
        root=str(project_root / "runs"),
    )
    registry.record_successful_run(pheno_id, pqtl_dataset, run_id, root=str(project_root / "runs"))
    print(f"[TRACKING] Recorded successful run in registry: {run_id}")
    print("[DONE] Local Docker run completed.")