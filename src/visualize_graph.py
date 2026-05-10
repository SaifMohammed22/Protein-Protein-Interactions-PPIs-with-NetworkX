import sys
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# Import from graph_builder
from graph_builder import parse_data, build_graph

def visualize_and_save(graph: nx.DiGraph, output_path: Path, title: str, show_labels: bool = False):
    """
    Visualizes the networkx graph and saves it to the specified output path.
    """
    print(f"Generating visualization for: {title} ({graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges)")
    
    if graph.number_of_nodes() == 0:
        print(f"Warning: Graph for '{title}' is empty. Skipping.")
        return

    plt.figure(figsize=(14, 14))
    plt.title(title, fontsize=20, pad=20)
    
    # Optimization: Use spring_layout with fewer iterations for large graphs
    if graph.number_of_nodes() < 200:
        pos = nx.kamada_kawai_layout(graph)
    else:
        print("  Computing layout... (this may take a moment)")
        # k controls the optimal distance between nodes
        pos = nx.spring_layout(graph, seed=42, k=1/np.sqrt(graph.number_of_nodes()), iterations=30)
    
    # Calculate node degrees for sizing
    degrees = dict(graph.degree())
    
    if show_labels:
        node_sizes = [v * 100 + 50 for v in degrees.values()]
    else:
        # Smaller nodes for large graphs
        node_sizes = [15 for _ in degrees]
    
    # Draw nodes
    nx.draw_networkx_nodes(
        graph, pos, 
        node_size=node_sizes, 
        node_color=list(degrees.values()), 
        cmap=plt.cm.Blues, 
        alpha=0.8
    )
    
    # Draw edges with transparency
    alpha = max(0.05, min(0.4, 1000 / max(1, graph.number_of_edges())))
    nx.draw_networkx_edges(
        graph, pos, 
        arrowstyle='->', 
        arrowsize=5 if graph.number_of_nodes() < 1000 else 1, 
        edge_color='gray', 
        alpha=alpha,
        width=0.5
    )
    
    if show_labels and graph.number_of_nodes() < 100:
        nx.draw_networkx_labels(graph, pos, font_size=8, font_family='sans-serif')
    
    plt.axis('off')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight') # Reduced DPI for speed
    plt.close()
    print(f"  Saved to {output_path}")

def get_largest_component(graph: nx.DiGraph):
    """Returns the largest weakly connected component as a subgraph."""
    if graph.number_of_nodes() == 0:
        return graph
    components = sorted(nx.weakly_connected_components(graph), key=len, reverse=True)
    return graph.subgraph(components[0]).copy()

def main():
    base_dir = Path(__file__).resolve().parent.parent
    data_path = base_dir / "data" / "PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    results_dir = base_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        print(f"Error: Data file not found at {data_path}")
        return

    # 1. Visualize a Sub-sample (Connected Component)
    print("--- Phase 1: Sub-sample Visualization ---")
    # Parse 1000 edges to find a decent sized component
    tails, heads, weights, dists, _ = parse_data(data_path, max_edges=1000, random_state=42)
    full_sub_graph = build_graph(tails, heads, weights, dists)
    
    sub_graph = get_largest_component(full_sub_graph)
    sub_sample_out = results_dir / "sub_sample_graph.png"
    visualize_and_save(sub_graph, sub_sample_out, "PPI Sub-sample (Largest Component)", show_labels=True)

    # 2. Visualize a Representative Large Sample
    print("\n--- Phase 2: Representative Sample Visualization ---")
    # Reduced to 3000 edges for faster rendering while still showing complexity
    tails_l, heads_l, weights_l, dists_l, _ = parse_data(data_path, max_edges=3000, random_state=42)
    large_graph = build_graph(tails_l, heads_l, weights_l, dists_l)
    
    large_graph_out = results_dir / "large_graph_sample.png"
    visualize_and_save(large_graph, large_graph_out, "PPI Representative Sample (3000 edges)", show_labels=False)

if __name__ == "__main__":
    main()
