"""
Connectivity Analysis Module (Member 2: The Statistician)
Handles node-level statistics, degree calculations, protein ranking, and visualization.
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def calculate_single_protein_degrees(graph: nx.Graph) -> dict:
    """
    Calculate the degree of each protein (node) in the graph.
    
    Args:
        graph: NetworkX graph object
        
    Returns:
        Dictionary with protein names as keys and their degrees as values
    """
    degrees = dict(graph.degree())
    return degrees


def rank_multiprotein_by_degree(graph: nx.Graph) -> list:
    """
    Rank all proteins by their degree (connectivity) in descending order.
    
    Args:
        graph: NetworkX graph object
        
    Returns:
        List of tuples (protein_name, degree) sorted by degree in descending order
    """
    degrees = dict(graph.degree())
    ranked_proteins = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    return ranked_proteins


def save_single_protein_degrees(graph: nx.Graph, output_path: str) -> None:
    """
    Calculate and save the degree of each protein to a text file.
    
    Args:
        graph: NetworkX graph object
        output_path: Path where results/2_single_protein_degree.txt will be saved
    """
    degrees = calculate_single_protein_degrees(graph)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write degrees to file
    with open(output_path, 'w') as f:
        f.write("Protein\tDegree\n")
        f.write("=" * 40 + "\n")
        for protein, degree in sorted(degrees.items()):
            f.write(f"{protein}\t{degree}\n")


def save_ranked_proteins(graph: nx.Graph, output_path: str) -> None:
    """
    Rank proteins by degree and save the ranked list to a text file.
    
    Args:
        graph: NetworkX graph object
        output_path: Path where results/3_multiprotein_ranked.txt will be saved
    """
    ranked = rank_multiprotein_by_degree(graph)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write ranked proteins to file
    with open(output_path, 'w') as f:
        f.write("Rank\tProtein\tDegree\n")
        f.write("=" * 50 + "\n")
        for rank, (protein, degree) in enumerate(ranked, 1):
            f.write(f"{rank}\t{protein}\t{degree}\n")


def plot_degree_distribution(graph: nx.Graph, output_path: str) -> None:
    """
    Calculate degree distribution and create a histogram visualization.
    
    Args:
        graph: NetworkX graph object
        output_path: Path where results/3_multiprotein_histogram.png will be saved
    """
    degrees = [degree for node, degree in graph.degree()]
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create histogram
    plt.figure(figsize=(10, 6))
    plt.hist(degrees, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    plt.xlabel('Degree (Number of Interactions)', fontsize=12)
    plt.ylabel('Frequency (Number of Proteins)', fontsize=12)
    plt.title('Protein Degree Distribution', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add statistics text
    mean_degree = np.mean(degrees)
    median_degree = np.median(degrees)
    max_degree = np.max(degrees)
    
    stats_text = f'Mean: {mean_degree:.2f}\nMedian: {median_degree:.2f}\nMax: {max_degree}'
    plt.text(0.98, 0.97, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def run_connectivity_analysis(graph: nx.Graph, results_dir: str = "results") -> None:
    """
    Run all connectivity analysis functions and generate all output files.
    
    Args:
        graph: NetworkX graph object
        results_dir: Directory where results will be saved (default: "results")
    """
    # Define output paths
    degree_output = os.path.join(results_dir, "2_single_protein_degree.txt")
    ranked_output = os.path.join(results_dir, "3_multiprotein_ranked.txt")
    histogram_output = os.path.join(results_dir, "3_multiprotein_histogram.png")
    
    # Run all analyses
    print("Calculating single protein degrees...")
    save_single_protein_degrees(graph, degree_output)
    print(f"✓ Saved to {degree_output}")
    
    print("Ranking proteins by connectivity...")
    save_ranked_proteins(graph, ranked_output)
    print(f"✓ Saved to {ranked_output}")
    
    print("Generating degree distribution histogram...")
    plot_degree_distribution(graph, histogram_output)
    print(f"✓ Saved to {histogram_output}")
    
    print("\nConnectivity analysis completed successfully!")


# Statistics helper functions
def get_degree_statistics(graph: nx.Graph) -> dict:
    """
    Get comprehensive degree statistics for the graph.
    
    Args:
        graph: NetworkX graph object
        
    Returns:
        Dictionary containing various degree statistics
    """
    degrees = [degree for node, degree in graph.degree()]
    
    stats = {
        'mean_degree': np.mean(degrees),
        'median_degree': np.median(degrees),
        'std_degree': np.std(degrees),
        'min_degree': np.min(degrees),
        'max_degree': np.max(degrees),
        'num_nodes': graph.number_of_nodes(),
        'num_edges': graph.number_of_edges()
    }
    
    return stats


def print_degree_statistics(graph: nx.Graph) -> None:
    """
    Print comprehensive degree statistics to console.
    
    Args:
        graph: NetworkX graph object
    """
    stats = get_degree_statistics(graph)
    
    print("\n" + "=" * 50)
    print("CONNECTIVITY STATISTICS")
    print("=" * 50)
    print(f"Number of nodes (proteins):    {stats['num_nodes']}")
    print(f"Number of edges (interactions): {stats['num_edges']}")
    print(f"Mean degree:                   {stats['mean_degree']:.2f}")
    print(f"Median degree:                 {stats['median_degree']:.2f}")
    print(f"Std deviation:                 {stats['std_degree']:.2f}")
    print(f"Min degree:                    {stats['min_degree']}")
    print(f"Max degree:                    {stats['max_degree']}")
    print("=" * 50 + "\n")


def get_protein_neighbors(graph: nx.DiGraph, protein_id: str, output_path: str) -> dict:
    """
    Get all directly connected proteins (neighbors) of a given protein.
    
    Returns both incoming and outgoing neighbors with their interaction weights.
    For a directed graph, this includes:
    - Predecessors: proteins that interact WITH this protein (incoming edges)
    - Successors: proteins that this protein interacts WITH (outgoing edges)
    
    Args:
        graph: NetworkX directed graph object
        protein_id: UniProt ID of the query protein
        output_path: Path where results will be saved (e.g., "results/0_protein_neighbors.txt")
        
    Returns:
        Dictionary with 'degree', 'predecessors', and 'successors' info
        
    Raises:
        nx.NodeNotFound: If protein_id not in graph
    """
    if protein_id not in graph:
        raise nx.NodeNotFound(f"Protein '{protein_id}' not found in graph.")
    
    # Get incoming and outgoing neighbors
    predecessors = dict(graph.pred[protein_id])  # Proteins pointing TO this one
    successors = dict(graph.succ[protein_id])    # Proteins this one points TO
    
    # Calculate total degree
    total_degree = len(predecessors) + len(successors)
    in_degree = len(predecessors)
    out_degree = len(successors)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("PROTEIN NEIGHBORS ANALYSIS\n")
        f.write(f"Query Protein: {protein_id}\n")
        f.write(f"Total Degree: {total_degree}\n")
        f.write(f"  - In-Degree (incoming interactions): {in_degree}\n")
        f.write(f"  - Out-Degree (outgoing interactions): {out_degree}\n")
        f.write("=" * 60 + "\n\n")
        
        # Write incoming neighbors (predecessors)
        if predecessors:
            f.write("--- INCOMING INTERACTIONS (Predecessors) ---\n")
            f.write("Protein ID\t\tInteraction Weight\n")
            f.write("-" * 40 + "\n")
            for neighbor_id, edge_data in sorted(predecessors.items()):
                weight = edge_data.get('weight', 'N/A')
                f.write(f"{neighbor_id}\t\t{weight}\n")
            f.write("\n")
        
        # Write outgoing neighbors (successors)
        if successors:
            f.write("--- OUTGOING INTERACTIONS (Successors) ---\n")
            f.write("Protein ID\t\tInteraction Weight\n")
            f.write("-" * 40 + "\n")
            for neighbor_id, edge_data in sorted(successors.items()):
                weight = edge_data.get('weight', 'N/A')
                f.write(f"{neighbor_id}\t\t{weight}\n")
            f.write("\n")
    
    return {
        'protein_id': protein_id,
        'total_degree': total_degree,
        'in_degree': in_degree,
        'out_degree': out_degree,
        'predecessors': predecessors,
        'successors': successors
    }
