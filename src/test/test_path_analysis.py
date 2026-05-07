import math
import os
import sys
from pathlib import Path

# Ensure the project root is on sys.path so `src.` imports work
# when running this file directly (e.g. python src/test/test_path_analysis.py)
_project_root = str(Path(__file__).resolve().parent.parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import matplotlib.pyplot as plt
import networkx as nx
from src.path_analysis import find_all_shortest_paths, analyze_shortest_paths, run_path_analysis


def generate_mock_graph() -> nx.DiGraph:
    """Create a small directed PPI graph with weighted edges for testing.

    Returns a ``nx.DiGraph`` with 6 nodes named after UniProt-style
    identifiers and directed, weighted edges that guarantee:

    * The graph is **fully connected** (every node is reachable from
      every other node when edge direction is ignored).
    * There are **at least two distinct acyclic paths** from the
      designated source ``"P10001"`` to the designated target
      ``"P10006"``, enabling shortest-path algorithm testing.

    Graph topology
    --------------
    The edges are arranged so that two clear routes exist between
    P10001 → P10006:

        Path A  (upper route):  P10001 → P10002 → P10004 → P10006
        Path B  (lower route):  P10001 → P10003 → P10005 → P10006

    Additional cross-links (P10002 → P10005, P10003 → P10004) tie the
    two branches together and keep the graph fully connected.

    ASCII sketch::

                  0.85        0.60
        P10001 -------→ P10002 -------→ P10004
          |               |    ↗            |
          |  0.90         | 0.40   0.75     | 0.95
          ↓               ↓                 ↓
        P10003 -------→ P10005          P10006
                  0.50          ←-----------↑
                                     (via P10005, wt 0.70)

    Edge weights represent interaction confidence scores (0.01 – 1.0).

    Testing notes
    -------------
    * **Source node for path tests:** ``"P10001"``
    * **Target node for path tests:** ``"P10006"``
    * Path A cumulative weight: 0.85 + 0.60 + 0.95 = 2.40
    * Path B cumulative weight: 0.90 + 0.50 + 0.70 = 2.10  ← lighter
    """

    # ── Initialise a directed graph ──────────────────────────────────
    graph = nx.DiGraph()

    # ── Define the 6 UniProt-style protein nodes ─────────────────────
    nodes = ["P10001", "P10002", "P10003", "P10004", "P10005", "P10006"]
    graph.add_nodes_from(nodes)

    # ── Define directed, weighted edges ──────────────────────────────
    # Each tuple is (tail, head, weight).
    # "tail → head" mirrors the biological convention where the tail
    # protein activates / signals to the head protein.
    edges = [
        # --- Upper route: P10001 → P10002 → P10004 → P10006 ---------
        ("P10001", "P10002", 0.85),   # high-confidence direct link
        ("P10002", "P10004", 0.60),   # moderate confidence
        ("P10004", "P10006", 0.95),   # very high confidence

        # --- Lower route: P10001 → P10003 → P10005 → P10006 ---------
        ("P10001", "P10003", 0.90),   # high confidence
        ("P10003", "P10005", 0.50),   # moderate confidence
        ("P10005", "P10006", 0.70),   # good confidence

        # --- Cross-links (keep graph fully connected) -----------------
        ("P10002", "P10005", 0.40),   # weak cross-branch link
        ("P10003", "P10004", 0.75),   # strong cross-branch link
    ]

    # Add all edges with their weight attribute in one call
    graph.add_weighted_edges_from(edges)

    # ── Compute 'distance' = -log(weight) for each edge ──────────────
    # The pathfinding algorithm (Dijkstra) minimises sum-of-distances,
    # which is mathematically equivalent to maximising the product of
    # raw probabilities.  See find_all_shortest_paths() for full
    # explanation of the negative-log transform.
    for u, v, data in graph.edges(data=True):
        data["distance"] = -math.log(data["weight"])

    return graph


if __name__ == "__main__":
    import time
    import random
    import traceback

    passed = 0
    failed = 0

    def report(name, ok, detail=""):
        global passed, failed
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}] {name}" + (f"  -- {detail}" if detail else ""))
        if ok:
            passed += 1
        else:
            failed += 1

    print("=" * 60)
    print("  PATH_ANALYSIS.PY  --  STRESS TESTS & DEBUGGING")
    print("=" * 60)

    # ------------------------------------------------------------------
    # TEST 1: Basic happy-path with the mock graph
    # ------------------------------------------------------------------
    print("\n--- TEST 1: Mock graph basic run ---")
    G = generate_mock_graph()
    results = find_all_shortest_paths(G, "P10001", "P10006")
    report("At least one path found", len(results) >= 1, f"{len(results)} path(s)")

    # Every result dict must have the three expected keys
    keys_ok = all(
        {"path", "edge_weights", "path_score"} <= set(r.keys())
        for r in results
    )
    report("Result dicts have correct keys", keys_ok)

    # Verify score = product of edge weights (within floating-point tolerance)
    for i, r in enumerate(results):
        expected = math.prod(r["edge_weights"])
        close = math.isclose(r["path_score"], expected, rel_tol=1e-9)
        report(f"Path {i+1} score matches product", close,
               f"score={r['path_score']:.8f}  expected={expected:.8f}")

    # Results should be sorted descending by path_score
    scores = [r["path_score"] for r in results]
    report("Paths sorted descending by score", scores == sorted(scores, reverse=True))

    # ------------------------------------------------------------------
    # TEST 2: No path exists (reversed direction)
    # ------------------------------------------------------------------
    print("\n--- TEST 2: No path exists (target -> source) ---")
    try:
        find_all_shortest_paths(G, "P10006", "P10001")
        report("NetworkXNoPath raised", False, "Exception was NOT raised")
    except nx.NetworkXNoPath:
        report("NetworkXNoPath raised", True)
    except Exception as e:
        report("NetworkXNoPath raised", False, f"Wrong exception: {type(e).__name__}")

    # ------------------------------------------------------------------
    # TEST 3: Missing node (NodeNotFound)
    # ------------------------------------------------------------------
    print("\n--- TEST 3: Missing source / target node ---")
    try:
        find_all_shortest_paths(G, "FAKE_NODE", "P10006")
        report("NodeNotFound for bad source", False)
    except nx.NodeNotFound:
        report("NodeNotFound for bad source", True)

    try:
        find_all_shortest_paths(G, "P10001", "FAKE_NODE")
        report("NodeNotFound for bad target", False)
    except nx.NodeNotFound:
        report("NodeNotFound for bad target", True)

    # ------------------------------------------------------------------
    # TEST 4: Trivial single-edge graph
    # ------------------------------------------------------------------
    print("\n--- TEST 4: Single-edge graph ---")
    G_tiny = nx.DiGraph()
    G_tiny.add_edge("A", "B", weight=0.42, distance=-math.log(0.42))
    res = find_all_shortest_paths(G_tiny, "A", "B")
    report("Single-edge path found", len(res) == 1)
    report("Score equals edge weight",
           math.isclose(res[0]["path_score"], 0.42, rel_tol=1e-9))

    # ------------------------------------------------------------------
    # TEST 5: Source == Target (zero-length path)
    # ------------------------------------------------------------------
    print("\n--- TEST 5: Source equals target ---")
    try:
        res = find_all_shortest_paths(G, "P10001", "P10001")
        # NetworkX returns the trivial path [source] with length 0
        report("Trivial self-path returned", len(res) == 1 and res[0]["path"] == ["P10001"])
        report("Self-path score is 1.0 (empty product)",
               math.isclose(res[0]["path_score"], 1.0, rel_tol=1e-9))
    except Exception as e:
        report("Self-path handling", False, f"{type(e).__name__}: {e}")

    # ------------------------------------------------------------------
    # TEST 6: Disconnected components (two isolated clusters)
    # ------------------------------------------------------------------
    print("\n--- TEST 6: Disconnected components ---")
    G_disc = nx.DiGraph()
    G_disc.add_edge("X1", "X2", weight=0.5, distance=-math.log(0.5))
    G_disc.add_edge("Y1", "Y2", weight=0.8, distance=-math.log(0.8))
    try:
        find_all_shortest_paths(G_disc, "X1", "Y2")
        report("NetworkXNoPath for disconnected", False)
    except nx.NetworkXNoPath:
        report("NetworkXNoPath for disconnected", True)

    # ------------------------------------------------------------------
    # TEST 7: Edge-weight boundary values (very small / very large)
    # ------------------------------------------------------------------
    print("\n--- TEST 7: Boundary edge weights ---")
    G_bnd = nx.DiGraph()
    G_bnd.add_edge("S", "M", weight=0.01, distance=-math.log(0.01))   # minimum allowed weight
    G_bnd.add_edge("M", "T", weight=1.0, distance=-math.log(1.0))    # maximum allowed weight
    res = find_all_shortest_paths(G_bnd, "S", "T")
    report("Boundary weights: path found", len(res) == 1)
    report("Boundary score correct",
           math.isclose(res[0]["path_score"], 0.01 * 1.0, rel_tol=1e-9))

    # ------------------------------------------------------------------
    # TEST 8: Multiple equal-cost shortest paths
    # ------------------------------------------------------------------
    print("\n--- TEST 8: Multiple equal-cost shortest paths ---")
    G_eq = nx.DiGraph()
    #  S --0.5--> A --0.5--> T   total distance = -log(0.5)*2
    #  S --0.5--> B --0.5--> T   total distance = -log(0.5)*2  (equal)
    G_eq.add_edge("S", "A", weight=0.5, distance=-math.log(0.5))
    G_eq.add_edge("A", "T", weight=0.5, distance=-math.log(0.5))
    G_eq.add_edge("S", "B", weight=0.5, distance=-math.log(0.5))
    G_eq.add_edge("B", "T", weight=0.5, distance=-math.log(0.5))
    res = find_all_shortest_paths(G_eq, "S", "T")
    report("Two equal-cost paths found", len(res) == 2,
           f"got {len(res)}")
    # Both scores should be 0.25 (0.5 * 0.5)
    all_eq = all(math.isclose(r["path_score"], 0.25) for r in res)
    report("Both scores equal 0.25", all_eq)

    # ------------------------------------------------------------------
    # TEST 9: Long chain (10-hop linear path)
    # ------------------------------------------------------------------
    print("\n--- TEST 9: Long linear chain (10 hops) ---")
    G_chain = nx.DiGraph()
    chain_weight = 0.9
    for i in range(10):
        G_chain.add_edge(f"N{i}", f"N{i+1}",
                         weight=chain_weight,
                         distance=-math.log(chain_weight))
    res = find_all_shortest_paths(G_chain, "N0", "N10")
    expected_score = chain_weight ** 10
    report("10-hop path found", len(res) == 1)
    report("10-hop score correct",
           math.isclose(res[0]["path_score"], expected_score, rel_tol=1e-9),
           f"score={res[0]['path_score']:.10f}  expected={expected_score:.10f}")

    # ------------------------------------------------------------------
    # TEST 10: Stress test -- large random graph
    # ------------------------------------------------------------------
    print("\n--- TEST 10: Stress test (large random graph) ---")
    random.seed(42)
    N_NODES = 500
    N_EDGES = 3000

    G_big = nx.DiGraph()
    nodes = [f"P{i:05d}" for i in range(N_NODES)]
    G_big.add_nodes_from(nodes)

    # Add random directed edges with random weights
    edges_added = 0
    while edges_added < N_EDGES:
        u = random.choice(nodes)
        v = random.choice(nodes)
        if u != v and not G_big.has_edge(u, v):
            w = round(random.uniform(0.01, 1.0), 4)
            G_big.add_edge(u, v, weight=w, distance=-math.log(w))
            edges_added += 1

    print(f"  Graph: {G_big.number_of_nodes()} nodes, {G_big.number_of_edges()} edges")

    # Pick two random nodes and try to find paths
    src, tgt = "P00000", "P00499"
    t0 = time.perf_counter()
    try:
        res = find_all_shortest_paths(G_big, src, tgt)
        elapsed = time.perf_counter() - t0
        report(f"Large graph: found {len(res)} path(s)", True,
               f"in {elapsed:.4f}s")
        report("Large graph: completed under 5s", elapsed < 5.0,
               f"{elapsed:.4f}s")
    except nx.NetworkXNoPath:
        elapsed = time.perf_counter() - t0
        report("Large graph: no path (expected for sparse random)", True,
               f"in {elapsed:.4f}s")
    except Exception as e:
        elapsed = time.perf_counter() - t0
        report("Large graph: unexpected error", False,
               f"{type(e).__name__}: {e} ({elapsed:.4f}s)")

    # ------------------------------------------------------------------
    # TEST 11: File output verification
    # ------------------------------------------------------------------
    print("\n--- TEST 11: File output sanity check ---")
    import tempfile
    import shutil

    tmp_dir = os.path.join(os.path.dirname(__file__), "__test_tmp__")
    G_mock = generate_mock_graph()
    analyze_shortest_paths(G_mock, "P10001", "P10006", tmp_dir)

    txt_file = os.path.join(tmp_dir, "1_shortest_paths.txt")
    png_file = os.path.join(tmp_dir, "1_shortest_paths_subnetwork.png")

    report("Text report created", os.path.isfile(txt_file))
    report("PNG plot created", os.path.isfile(png_file))

    # Check text file is non-empty and contains expected header
    if os.path.isfile(txt_file):
        with open(txt_file, "r", encoding="utf-8") as f:
            content = f.read()
        report("Report contains source ID",   "P10001" in content)
        report("Report contains target ID",   "P10006" in content)
        report("Report contains 'Path Score'", "Path Score" in content)

    # Check PNG is a valid file (starts with PNG magic bytes)
    if os.path.isfile(png_file):
        with open(png_file, "rb") as f:
            header = f.read(8)
        is_png = header[:4] == b"\x89PNG"
        report("PNG has valid header", is_png)

    # ------------------------------------------------------------------
    # TEST 12: run_path_analysis graceful error handling
    # ------------------------------------------------------------------
    print("\n--- TEST 12: run_path_analysis wrapper ---")
    result = run_path_analysis(G, "DOES_NOT_EXIST", "P10006", tmp_dir)
    report("Returns None for missing node", result is None)

    result = run_path_analysis(G, "P10006", "P10001", tmp_dir)
    report("Returns None for no-path case", result is None)

    # ------------------------------------------------------------------
    # TEST 13: Debugging helper -- print full graph structure
    # ------------------------------------------------------------------
    print("\n--- DEBUG: Mock graph structure ---")
    G_dbg = generate_mock_graph()
    print(f"  Nodes ({G_dbg.number_of_nodes()}): {sorted(G_dbg.nodes)}")
    print(f"  Edges ({G_dbg.number_of_edges()}):")
    for u, v, d in G_dbg.edges(data=True):
        print(f"    {u} -> {v}  weight={d['weight']}  distance={d['distance']:.4f}")
    print(f"  Connected (undirected): {nx.is_connected(G_dbg.to_undirected())}")
    all_paths = list(nx.all_simple_paths(G_dbg, "P10001", "P10006"))
    print(f"  All simple paths P10001->P10006: {len(all_paths)}")
    for p in all_paths:
        ws = [G_dbg[p[i]][p[i+1]]["weight"] for i in range(len(p)-1)]
        print(f"    {' -> '.join(p)}  weights={ws}  product={math.prod(ws):.6f}")

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    total = passed + failed
    print(f"  RESULTS: {passed}/{total} passed, {failed} failed")
    if failed == 0:
        print("  ALL TESTS PASSED")
    else:
        print("  ** SOME TESTS FAILED -- review output above **")
    print("=" * 60)
