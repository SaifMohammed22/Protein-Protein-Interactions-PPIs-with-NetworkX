"""
Test script for connectivity_analysis.py (Member 2)
Tests all functions with a fake 5-node graph
"""

import networkx as nx
import sys
import os
from pathlib import Path

# Add src to path so we can import connectivity_analysis
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from connectivity_analysis import (
    calculate_single_protein_degrees,
    rank_multiprotein_by_degree,
    save_single_protein_degrees,
    save_ranked_proteins,
    plot_degree_distribution,
    get_degree_statistics,
    print_degree_statistics,
    run_connectivity_analysis
)


def create_fake_graph():
    """Create a simple 5-node fake graph for testing.
    
    Nodes: ProteinA, ProteinB, ProteinC, ProteinD, ProteinE
    Edges: A-B, B-C, C-D, D-E, A-C, B-D
    
    Degrees: A=2, B=3, C=3, D=3, E=1
    """
    G = nx.Graph()
    edges = [
        ('ProteinA', 'ProteinB'),
        ('ProteinB', 'ProteinC'),
        ('ProteinC', 'ProteinD'),
        ('ProteinD', 'ProteinE'),
        ('ProteinA', 'ProteinC'),
        ('ProteinB', 'ProteinD')
    ]
    G.add_edges_from(edges)
    return G


def test_individual_functions():
    """Test each function individually."""
    print("\n" + "="*60)
    print("TEST 1: Individual Function Tests")
    print("="*60)
    
    # Create fake graph
    G = create_fake_graph()
    print(f"\n[OK] Created fake 5-node graph")
    print(f"  - Nodes: {list(G.nodes())}")
    print(f"  - Edges: {list(G.edges())}")
    
    # Test degree calculation
    print("\n--- Test: calculate_single_protein_degrees() ---")
    degrees = calculate_single_protein_degrees(G)
    print(f"Degrees calculated: {degrees}")
    assert len(degrees) == 5, "Should have 5 proteins"
    print("[OK] PASSED")
    
    # Test protein ranking
    print("\n--- Test: rank_multiprotein_by_degree() ---")
    ranked = rank_multiprotein_by_degree(G)
    print(f"Ranked proteins (top 3): {ranked[:3]}")
    assert len(ranked) == 5, "Should have 5 proteins"
    assert ranked[0][1] >= ranked[-1][1], "Should be sorted descending"
    print("[OK] PASSED")
    
    # Test statistics
    print("\n--- Test: get_degree_statistics() ---")
    stats = get_degree_statistics(G)
    print(f"Mean degree: {stats['mean_degree']:.2f}")
    print(f"Median degree: {stats['median_degree']:.2f}")
    print(f"Max degree: {stats['max_degree']}")
    assert stats['num_nodes'] == 5, "Should have 5 nodes"
    assert stats['num_edges'] == 6, "Should have 6 edges"
    print("[OK] PASSED")
    
    print_degree_statistics(G)


def test_output_files():
    """Test that output files are created correctly."""
    print("\n" + "="*60)
    print("TEST 2: Output File Generation")
    print("="*60)
    
    # Create results directory
    results_dir = Path(__file__).parent / "__test_tmp__"
    results_dir.mkdir(exist_ok=True)
    
    # Create fake graph
    G = create_fake_graph()
    
    # Test file outputs
    print("\n--- Test: save_single_protein_degrees() ---")
    degree_file = results_dir / "2_single_protein_degree.txt"
    save_single_protein_degrees(G, str(degree_file))
    assert degree_file.exists(), f"File {degree_file} should exist"
    with open(degree_file) as f:
        content = f.read()
        assert "Protein" in content and "Degree" in content
    print(f"[OK] Created: {degree_file.name}")
    print(f"  File size: {degree_file.stat().st_size} bytes")
    
    print("\n--- Test: save_ranked_proteins() ---")
    ranked_file = results_dir / "3_multiprotein_ranked.txt"
    save_ranked_proteins(G, str(ranked_file))
    assert ranked_file.exists(), f"File {ranked_file} should exist"
    with open(ranked_file) as f:
        content = f.read()
        assert "Rank" in content and "Degree" in content
    print(f"[OK] Created: {ranked_file.name}")
    print(f"  File size: {ranked_file.stat().st_size} bytes")
    
    print("\n--- Test: plot_degree_distribution() ---")
    histogram_file = results_dir / "3_multiprotein_histogram.png"
    plot_degree_distribution(G, str(histogram_file))
    assert histogram_file.exists(), f"File {histogram_file} should exist"
    print(f"[OK] Created: {histogram_file.name}")
    print(f"  File size: {histogram_file.stat().st_size} bytes")


def test_full_pipeline():
    """Test the full pipeline run_connectivity_analysis()."""
    print("\n" + "="*60)
    print("TEST 3: Full Pipeline (run_connectivity_analysis)")
    print("="*60)
    
    # Create results directory
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    # Create fake graph
    G = create_fake_graph()
    
    print("\nRunning full connectivity analysis pipeline...")
    run_connectivity_analysis(G, str(results_dir))
    
    # Verify all output files exist
    expected_files = [
        "2_single_protein_degree.txt",
        "3_multiprotein_ranked.txt",
        "3_multiprotein_histogram.png"
    ]
    
    print("\nVerifying output files:")
    for filename in expected_files:
        filepath = results_dir / filename
        assert filepath.exists(), f"Expected file {filename} not created!"
        size = filepath.stat().st_size
        print(f"  [OK] {filename} ({size} bytes)")
    
    print("\n[OK] PASSED - All output files created successfully!")


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("CONNECTIVITY ANALYSIS MODULE - COMPREHENSIVE TEST SUITE")
    print("="*60)
    
    try:
        test_individual_functions()
        test_output_files()
        test_full_pipeline()
        
        print("\n" + "="*60)
        print("[OK][OK][OK] ALL TESTS PASSED! [OK][OK][OK]")
        print("="*60)
        print("\nYour connectivity_analysis.py is working correctly!")
        print("Output files are in the 'results/' directory")
        
    except AssertionError as e:
        print(f"\n[FAIL] TEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
