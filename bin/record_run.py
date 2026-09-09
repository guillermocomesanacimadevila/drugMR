import argparse

from drugmr import registry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--run_id", required=True)
    p.add_argument("--root", default="runs")
    args = p.parse_args()
    registry.record_successful_run(args.pheno_id, args.pqtl_dataset, args.run_id, root=args.root)
    print(f"[DONE] Recorded successful run in registry: {args.run_id}")


if __name__ == "__main__":
    main()
