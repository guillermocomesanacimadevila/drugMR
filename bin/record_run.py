import argparse
from datetime import datetime

from drugmr import registry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--run_id", required=True)
    p.add_argument("--root", default="runs")
    p.add_argument("--git_sha7", default="unknown")
    p.add_argument("--image_name", default="unknown")
    p.add_argument("--host", default="unknown")
    args = p.parse_args()

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
        },
        root=args.root,
    )
    registry.record_successful_run(args.pheno_id, args.pqtl_dataset, args.run_id, root=args.root)
    print(f"[DONE] Recorded successful run in registry: {args.run_id}")


if __name__ == "__main__":
    main()
