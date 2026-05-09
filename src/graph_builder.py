from pathlib import Path
import math

import networkx as nx
import numpy as np


# Parse graph data from PathLinker file
def parse_data(data_path: str, sample_frac: float | None = None, max_edges: int | None = None, random_state: int | None = None):
    """
    Parse graph data from a PathLinker file.
    Returns: tails, heads, edge_weights, edge_distances, edge_types
    """
    tails = []
    heads = []
    edge_weights = []
    edge_distances = []
    edge_types = []
    with open(data_path, "r") as file_handle:
        for line in file_handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            tail, head, weight, edge_type = line.split("\t")
            probability = float(weight)
            distance_val = - \
                math.log(probability) if probability > 0 else float('inf')
            tails.append(tail)
            heads.append(head)
            edge_weights.append(probability)
            edge_distances.append(distance_val)
            edge_types.append(edge_type)

    n = len(tails)
    if n == 0:
        return tails, heads, edge_weights, edge_distances, edge_types

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
        edge_distances = [edge_distances[i] for i in indices]
        edge_types = [edge_types[i] for i in indices]

    return tails, heads, edge_weights, edge_distances, edge_types


# Build the protein-protein graph using NetworkX
def build_graph(tails, heads, edge_weights, edge_distances=None):
    graph = nx.DiGraph()
    if edge_distances is None:
        edge_distances = [None] * len(tails)
    for tail, head, weight, distance in zip(tails, heads, edge_weights, edge_distances):
        graph.add_edge(tail, head, weight=weight, distance=distance)
    return graph


def build_demo_graph():
    graph = nx.DiGraph()
    demo_edges = [
        ("P1", "P2", 0.8),
        ("P1", "P3", 0.9),
        ("P2", "P4", 0.7),
        ("P3", "P4", 0.85),
        ("P4", "P5", 0.6),
    ]
    for tail, head, prob in demo_edges:
        distance = -math.log(prob) if prob > 0 else float('inf')
        graph.add_edge(tail, head, weight=prob, distance=distance)
    return graph


def export_adjacency_matrix(graph: nx.DiGraph, output_path: str):
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    adjacency_matrix = nx.to_numpy_array(graph, weight=None, dtype=int)
    np.savetxt(output_file, adjacency_matrix, delimiter=",", fmt="%d")
    return adjacency_matrix


if __name__ == "__main__":
    # Load and process the PathLinker dataset
    data_path = Path(__file__).resolve().parent.parent / \
        "data" / "PathLinker_2018_human-ppi-weighted-cap0_75.txt"

    # Parse a small sample for testing
    tails, heads, edge_weights, edge_distances, edge_types = parse_data(
        data_path,
        sample_frac=None,
        max_edges=100,
        random_state=42,
    )

    # Build the graph and export adjacency matrix
    graph = build_graph(tails, heads, edge_weights, edge_distances)
    export_adjacency_matrix(
        graph,
        Path(__file__).resolve().parent.parent /
        "results" / "4_adjacency_matrix.csv",
    )
