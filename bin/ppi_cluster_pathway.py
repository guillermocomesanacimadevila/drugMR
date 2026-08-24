import argparse
import polars as pl
from drugmr.string_ppi import StringPPI
from drugmr.unsupervised_algorithms import MarkovClustering
from drugmr.enrichr import EnrichR
from drugmr.utils import filter_mr_targets, filter_coloc_targets, filter_phewas, strip_protein_id
from drugmr.paths import (
    mr_out,
    coloc_out,
    phewas_out,
    pathway_candidate_targets_out,
    pathway_string_ppi_out,
    pathway_mcl_clusters_out,
    pathway_mcl_plot_out,
    pathway_enrichr_out,
)


# TO DO's
# Save into datasets and update postgres DB accodingly with new tables
# Add score gait into params/ for STRING
# Include changes into local.py and hpc.py wrappers within drugmr
# Adapt dashboard code (new page)
# Re-build docker image

def extract_candidate_targets(pqtl_dataset: str, pheno_id: str, out_dir: str = "results"):

    """
    Reqs:
    - cis-MR beta P_FDR (wald/IVW) < 0.05
    - Coloc (default priors) PP.H4 > tresh
    - Pass FinnGen safety -> UKBB (not needed)
    """

    mr_df = pl.read_csv(mr_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=out_dir), separator="\t")
    coloc_df = pl.read_csv(coloc_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=out_dir), separator="\t")
    phewas_df = pl.read_csv(phewas_out(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=out_dir), separator="\t")
    mr_targets = set(strip_protein_id(filter_mr_targets(mr_df)))
    coloc_targets = set(strip_protein_id(filter_coloc_targets(coloc_df)))
    phewas_targets = set(strip_protein_id(filter_phewas(phewas_df)))
    final_targets = sorted(mr_targets & coloc_targets & phewas_targets)
    out_file = pathway_candidate_targets_out(pqtl_dataset, pheno_id, out_dir)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame({"protein": final_targets}).write_csv(str(out_file), separator="\t")
    print(f"[DONE] Saved candidate targets: {out_file}")
    return final_targets


def query_and_filter_ppi(targets: list[str], string_score: float, pqtl_dataset: str, pheno_id: str, out_dir: str = "results"):

    """
    Query STRING PPI for the given targets and filter by score.
    """

    ppi_network = StringPPI().query_ppi(targets) # just a df for STRING results
    # add score gait into params/ for STRING ##################
    ppi_network = ppi_network.filter(pl.col("score").cast(pl.Float64) > string_score)
    for row in ppi_network.iter_rows(named=True):
        gene1 = row["preferredName_A"]
        gene2 = row["preferredName_B"]
        sc = row["score"]
        print(f"[TRACKING] STRING score between {gene1} and {gene2}: {sc}")

    out_file = pathway_string_ppi_out(pqtl_dataset, pheno_id, out_dir)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    ppi_network.write_csv(str(out_file), separator="\t")
    print(f"[DONE] Saved filtered STRING PPI: {out_file}")
    return ppi_network


def build_ppi_network(pqtl_dataset: str, pheno_id: str, string_score: float, out_dir: str = "results"):

    """
    Grab candidate targets from previous function -> build + filter STRING PPI.
    """

    targets = extract_candidate_targets(pqtl_dataset=pqtl_dataset, pheno_id=pheno_id, out_dir=out_dir)
    return query_and_filter_ppi(
        targets=targets,
        string_score=string_score,
        pqtl_dataset=pqtl_dataset,
        pheno_id=pheno_id,
        out_dir=out_dir,
    )


def perform_mcl_clustering(string_ppi: pl.DataFrame, pqtl_dataset: str, pheno_id: str, out_dir: str = "results"):

    """
    MCL clustering
    """

    mcl = MarkovClustering()

    matrix = mcl.df_to_matrix(
        gene_1_col=string_ppi["preferredName_A"],
        gene_2_col=string_ppi["preferredName_B"],
        score_col=string_ppi["score"]
    )

    result, index_clusters, gene_clusters, Q = mcl.markov_clustering(matrix=matrix)

    for cluster in gene_clusters:
        print(f"[TRACKING] MCL cluster ({len(cluster)} genes): {cluster}")

    plot_file = pathway_mcl_plot_out(pqtl_dataset, pheno_id, out_dir)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    mcl.draw_clusters(matrix, index_clusters, output_path=str(plot_file))
    print(f"[DONE] Saved MCL cluster plot: {plot_file}")

    clusters_file = pathway_mcl_clusters_out(pqtl_dataset, pheno_id, out_dir)
    clusters_df = pl.DataFrame({
        "cluster_id": [i for i, cluster in enumerate(gene_clusters) for _ in cluster],
        "gene": [gene for cluster in gene_clusters for gene in cluster],
    })
    clusters_df.write_csv(str(clusters_file), separator="\t")
    print(f"[DONE] Saved MCL clusters ({len(gene_clusters)} clusters, modularity Q={Q}): {clusters_file}")
    return gene_clusters


def run_enrichr(cluster_list: list[list], pqtl_dataset: str, pheno_id: str, out_dir: str = "results"):

    """
    For each cluster -> run EnrichR
    """

    enrichr = EnrichR()
    results = []
    for cluster_id, cluster in enumerate(cluster_list):
        go_df, kegg_df = enrichr.run_enrichr(gene_list=cluster)
        print(f"[TRACKING] EnrichR done for cluster ({len(cluster)} genes): {cluster}")
        go_file = pathway_enrichr_out(pqtl_dataset, pheno_id, cluster_id, "GO", out_dir)
        go_file.parent.mkdir(parents=True, exist_ok=True)
        go_df.write_csv(str(go_file), separator="\t")
        kegg_file = pathway_enrichr_out(pqtl_dataset, pheno_id, cluster_id, "KEGG", out_dir)
        kegg_df.write_csv(str(kegg_file), separator="\t")
        print(f"[DONE] Saved EnrichR results for cluster {cluster_id}: {go_file}, {kegg_file}")
        results.append((cluster, go_df, kegg_df))
    return results


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--string_score", type=float, default=0.4)
    p.add_argument("--out_dir", default="results")
    args = p.parse_args()

    ppi_network = build_ppi_network(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        string_score=args.string_score,
        out_dir=args.out_dir,
    )

    gene_clusters = perform_mcl_clustering(
        ppi_network,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        out_dir=args.out_dir,
    )

    results = run_enrichr(
        gene_clusters,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        out_dir=args.out_dir,
    )
    
    for cluster, go_df, kegg_df in results:
        print(f"[TRACKING] cluster: {cluster}")
        print(go_df)
        print(kegg_df)


if __name__ == "__main__":
    main()