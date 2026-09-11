import argparse
import json
from datetime import datetime
import yaml

from drugmr import paths, registry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--run_id", required=True)
    p.add_argument("--root", default="runs")
    p.add_argument("--git_sha7", default="unknown")
    p.add_argument("--image_name", default="unknown")
    p.add_argument("--host", default="unknown")
    p.add_argument("--params_json", required=True)
    args = p.parse_args()

    run_params = json.loads(args.params_json)
    params_lock = paths.run_params_lock_path(args.run_id, root=args.root)
    params_lock.parent.mkdir(parents=True, exist_ok=True)
    params_lock.write_text(yaml.safe_dump(run_params, sort_keys=False))

    # load a Nextflow-produced run too.
    registry.write_manifest(
        args.run_id,
        {
            "pheno_id": args.pheno_id,
            "pqtl_dataset": args.pqtl_dataset,
            "git_sha7": args.git_sha7,
            "date": datetime.now().strftime("%Y%m%d"),
            "created_at": datetime.now().isoformat(),
            "mode": "nextflow",
            "image_name": args.image_name,
            "host": args.host,
            "overwrite": False,
            "params_lock": params_lock.name,
        },
        root=args.root,
    )
    registry.record_successful_run(args.pheno_id, args.pqtl_dataset, args.run_id, root=args.root)
    print(f"[DONE] Recorded successful run in registry: {args.run_id}")


if __name__ == "__main__":
    main()
