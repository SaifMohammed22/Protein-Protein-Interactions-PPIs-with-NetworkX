from pathlib import Path
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


# Parse graph data from PathLinker file
def parse_data(data_path: str, sample_frac: float | None = None, max_edges: int | None = None, random_state: int | None = None):
    """Parse graph data from a PathLinker file.

    Optional subsampling can be applied by fraction (`sample_frac`) or by
    maximum number of edges (`max_edges`). If both are provided, `max_edges`
    takes precedence. `random_state` ensures reproducible sampling.
    """
    tails = []
    heads = []
    edge_weights = []
    edge_types = []
    with open(data_path, "r") as file_handle:
        for line in file_handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            tail, head, weight, edge_type = line.split("\t")
            tails.append(tail)
            heads.append(head)
            edge_weights.append(float(weight))
            edge_types.append(edge_type)

    n = len(tails)
    if n == 0:
        return tails, heads, edge_weights, edge_types

    # Determine number of samples
    if max_edges is not None:
        k = min(max_edges, n)
    elif sample_frac is not None:
        k = max(1, int(np.floor(sample_frac * n)))
    else:
        k = n

    if k < n:
        rng = np.random.RandomState(random_state)
        indices = rng.choice(n, size=k, replace=False)
        indices.sort()
        tails = [tails[i] for i in indices]
        heads = [heads[i] for i in indices]
        edge_weights = [edge_weights[i] for i in indices]
        edge_types = [edge_types[i] for i in indices]

    return tails, heads, edge_weights, edge_types


# Build the protein-protein graph using NetworkX
def build_graph(tails, heads, edge_weights):
    graph = nx.DiGraph()
    graph.add_weighted_edges_from(zip(tails, heads, edge_weights))
    return graph


def export_adjacency_matrix(graph: nx.DiGraph, output_path: str):
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    adjacency_matrix = nx.to_numpy_array(graph, weight=None, dtype=int)
    np.savetxt(output_file, adjacency_matrix, delimiter=",", fmt="%d")
    return adjacency_matrix


def visualize_graph(graph: nx.DiGraph, output_path: str | None = None):
    plt.figure(figsize=(14, 10))
    positions = nx.circular_layout(graph)

    nx.draw_networkx_nodes(
        graph, positions, node_size=30, node_color="#1a1d20")
    nx.draw_networkx_edges(graph, positions, arrows=True, alpha=0.2, width=0.5)
    plt.axis("off")
    plt.tight_layout()

    if output_path is not None:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    data_path = Path(__file__).resolve().parent.parent / \
        "data" / "PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    # Allow quick subsampling for large datasets via environment variables.
    sample_frac = 0.1
    max_edges = 5000
    random_state = 42
    sample_frac_val = float(sample_frac) if sample_frac is not None else None
    max_edges_val = int(max_edges) if max_edges is not None else None
    random_state_val = int(random_state) if random_state is not None else None

    tails, heads, edge_weights, edge_types = parse_data(
        data_path,
        sample_frac=sample_frac_val,
        max_edges=max_edges_val,
        random_state=random_state_val,
    )
    graph = build_graph(tails, heads, edge_weights)
    export_adjacency_matrix(
        graph,
        Path(__file__).resolve().parent.parent /
        "results" / "4_adjacency_matrix.csv",
    )
    visualize_graph(
        graph,
        Path(__file__).resolve().parent.parent /
        "images" / "image.png",
    )
