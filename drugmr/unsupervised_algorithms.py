import markov_clustering as mc
import networkx as nx
import polars as pl
import numpy as np 
import scipy.sparse as sp
import matplotlib.pyplot as plt 


class MarkovClustering:

    """
    MCL (network-based) clustering for STRING-derived PPIs
    """

    def __init__(self):
        self.network = nx.Graph()
        self.node_order = []


    def df_to_matrix(self, gene_1_col: list[str], gene_2_col: list[str], score_col: list[float]):

        """
        For each row in df (i.e STRING results) - asumming they're filtered by conf,
        grab each pair of proteins (i.e. network walk) and  create a df ->
        pl.DataFrame({
            "nameA": ["GENE1", "GENE2"],
            "nameB": ["GENEX", "GENEY"],
            "score": ["1", "2"]
        })
        """

        network_df = pl.DataFrame({
            "preferredName_A": gene_1_col,
            "preferredName_B": gene_2_col,
            "score": score_col
        })

        for row in network_df.iter_rows(named=True):
            self.network.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

        self.node_order = list(self.network.nodes())
        matrix = sp.csr_matrix(nx.to_scipy_sparse_array(self.network, nodelist=self.node_order, weight="weight", format="csr"))
        return matrix
    

    # def markov_clustering(self, matrix, inflation: float = 2.0):
    def markov_clustering(self, matrix, inflation: float = 2.0): # change to -> looping over inflation coefficients

        """
        Run MCL and convert the resulting matrix indices back into gene names.
        """

        result = mc.run_mcl(matrix, inflation=inflation)
        index_clusters = mc.get_clusters(result)
        gene_clusters = [[self.node_order[index] for index in cluster] for cluster in index_clusters]
        Q = mc.modularity(matrix=result, clusters=index_clusters)
        print("inflation:", inflation, "modularity:", Q)
        return result, index_clusters, gene_clusters, Q


    def draw_clusters(self, matrix, index_clusters, output_path: str = "string_mcl_clusters.png"):

        """
        Draw and save the clustered STRING network.
        """

        gene_positions = nx.spring_layout(self.network, seed=42)
        index_positions = {index: gene_positions[gene] for index, gene in enumerate(self.node_order)}
        labels = {index: gene for index, gene in enumerate(self.node_order)}
        mc.draw_graph(matrix, index_clusters, pos=index_positions, labels=labels, node_size=800, font_size=8, with_labels=True, edge_color="silver")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()