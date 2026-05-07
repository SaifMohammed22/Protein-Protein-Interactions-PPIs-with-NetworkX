import os
import sys
from pathlib import Path

# Ensure the current directory is in the path so we can import modules directly
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from graph_builder import parse_data, build_graph, export_adjacency_matrix
from connectivity_analysis import run_connectivity_analysis, print_degree_statistics
from id_converter import get_gene_names

# --- Phase 5: Test Variables ---
# Recommended test proteins from the implementation plan
PROTEIN_A = "P15056" 
PROTEIN_B = "P08069" 
PROTEIN_LIST = ["P15056", "P08069", "O15111"]

def main():
    # 1. Setup paths
    base_dir = Path(__file__).resolve().parent.parent
    data_path = base_dir / "data" / "PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    results_dir = base_dir / "results"
    
    # Ensure results directory exists
    os.makedirs(results_dir, exist_ok=True)

    print("=== PPI Network Analysis Pipeline ===")
    
    # 2. Phase 1: Foundation & Data Parsing
    print("\n--- Phase 1: Building Graph ---")
    # Subsampling for performance during testing (optional)
    tails, heads, weights, _ = parse_data(str(data_path), max_edges=10000)
    graph = build_graph(tails, heads, weights)
    
    # Export adjacency matrix
    adj_matrix_path = results_dir / "4_adjacency_matrix.csv"
    export_adjacency_matrix(graph, str(adj_matrix_path))
    print(f"Graph built with {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges.")
    print(f"Adjacency matrix saved to {adj_matrix_path}")

    # 3. Phase 2: Connectivity & Statistics
    print("\n--- Phase 2: Connectivity Analysis ---")
    run_connectivity_analysis(graph, str(results_dir))
    print_degree_statistics(graph)

    # 4. Phase 3: Pathfinding (Member 3 - Placeholder)
    print("\n--- Phase 3: Pathfinding ---")
    try:
        from path_analysis import run_path_analysis
        run_path_analysis(graph, PROTEIN_A, PROTEIN_B, str(results_dir))
    except (ImportError, AttributeError):
        print("[Notice] Phase 3 (path_analysis) is not yet implemented or lacks 'run_path_analysis'. Skipping.")

    # 5. Phase 4: Biological Context (ID Mapping)
    print("\n--- Phase 4: Biological Context ---")
    gene_map = get_gene_names(PROTEIN_LIST)
    print("\nProtein ID to Gene Name Mapping:")
    for uid, gene in gene_map.items():
        print(f"  {uid} -> {gene}")

    print("\n=== Pipeline Execution Completed Successfully ===")

if __name__ == "__main__":
    main()
