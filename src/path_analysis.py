import math
import os

import matplotlib.pyplot as plt
import networkx as nx


# ═══════════════════════════════════════════════════════════════════════
#  PHASE 1 & 2: PATHFINDING AND SCORING LOGIC
# ═══════════════════════════════════════════════════════════════════════
def find_all_shortest_paths(
    graph: nx.DiGraph,
    source: str,
    target: str,
) -> list[dict]:
    """Find every shortest simple path between two proteins and score them.

    We use `networkx.all_shortest_paths` with `weight="weight"`
    rather than `all_simple_paths` because we specifically want the
    shortest (minimum-weight) set, not the exhaustive set of all
    possible acyclic walks. Internally NetworkX delegates to Dijkstra's
    algorithm, which is well-suited for non-negative edge weights.

    Each edge weight is a probability (confidence that the interaction
    truly exists). For a complete signalling cascade (path) to occur,
    every interaction along the route must be present. Assuming
    statistical independence, the joint probability is the product of
    the individual edge probabilities: `Path Score = w₁ * w₂ * ... * wₙ`

    A higher path score indicates a more reliable signalling route.

    Args:
        graph (nx.DiGraph): The PPI directed graph. Every edge must carry a weight attribute (float, 0 < w <= 1).
        source (str): UniProt-style ID of the source (signal-origin) protein.
        target (str): UniProt-style ID of the target (signal-destination) protein.

    Returns:
        list[dict]: Each dict contains:
            - "path": list[str] -- ordered node sequence
            - "edge_weights": list[float] -- weight of each edge in the path
            - "path_score": float -- product of all edge weights

    Raises:
        nx.NodeNotFound: If source or target is not present in the graph.
        nx.NetworkXNoPath: If no directed path exists from source to target.
    """

    # ── Guard: make sure both endpoints exist in the graph ───────────
    for node_id, label in ((source, "Source"), (target, "Target")):
        if node_id not in graph:
            raise nx.NodeNotFound(
                f"{label} protein '{node_id}' is not present in the graph."
            )

    # ── Find all shortest paths using Dijkstra's algorithm ───────────
    # networkx.all_shortest_paths returns a generator of node-lists.
    # The `weight` parameter tells Dijkstra which edge attribute to minimise.
    # We minimise 'distance' (= -log(probability))
    # If no path exists, NetworkXNoPath is raised.
    try:
        shortest_paths_gen = nx.all_shortest_paths(
            graph,
            source=source,
            target=target,
            weight="distance",        # Minimise sum of -log(prob)
        )
        # Materialise the generator into a list so we can iterate twice
        shortest_paths = list(shortest_paths_gen)
    except nx.NetworkXNoPath:
        raise nx.NetworkXNoPath(
            f"No directed path exists from '{source}' to '{target}'."
        )

    # ── Score each path ──────────────────────────────────────────────
    scored_paths: list[dict] = []

    for path_nodes in shortest_paths:
        # Extract the weight for every consecutive edge in the path.
        # For a path [A, B, C] we need edges (A→B) and (B→C).
        edge_weights: list[float] = []
        for i in range(len(path_nodes) - 1):
            tail = path_nodes[i]
            head = path_nodes[i + 1]
            # Access the edge-data dict and pull the 'weight' attribute.
            weight = graph[tail][head]["weight"]
            edge_weights.append(weight)

        # Multiplicative score: joint probability of the interaction cascade.
        path_score = math.prod(edge_weights)

        scored_paths.append({
            "path":         path_nodes,
            "edge_weights": edge_weights,
            "path_score":   path_score,
        })

    # Sort descending by path_score so the most reliable route is first
    scored_paths.sort(key=lambda entry: entry["path_score"], reverse=True)

    return scored_paths


# ═══════════════════════════════════════════════════════════════════════
#  PHASE 3: FILE OUTPUT
# ═══════════════════════════════════════════════════════════════════════
def save_paths_to_file(
    scored_paths: list[dict],
    source: str,
    target: str,
    output_path: str,
) -> None:
    """Write scored shortest-path results to a human-readable text file.

    The output format is:

        ========================================
        SHORTEST PATH ANALYSIS
        Source Protein: P10001
        Target Protein: P10006
        Number of shortest paths found: 2
        ========================================

        --- Path 1 (Best) ---
        Node sequence : P10001 → P10002 → P10006
        Edge weights  : [0.85, 0.60]
        Path Score    : 0.510000
        ...

    Args:
        scored_paths (list[dict]): Output of `find_all_shortest_paths`.
        source (str): Source protein identifier.
        target (str): Target protein identifier.
        output_path (str): Path for the output text file.
    """

    # Ensure the parent directory tree exists before writing
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Use utf-8 encoding so directional arrows render correctly on
    # all platforms (Windows default cp1252 cannot encode '\u2192').

    with open(output_path, "w", encoding="utf-8") as fh:
        # ── Header block ─────────────────────────────────────────────
        fh.write("=" * 55 + "\n")
        fh.write("SHORTEST PATH ANALYSIS\n")
        fh.write(f"Source Protein : {source}\n")
        fh.write(f"Target Protein : {target}\n")
        fh.write(f"Number of shortest paths found : {len(scored_paths)}\n")
        fh.write("=" * 55 + "\n\n")

        # ── One block per path ───────────────────────────────────────
        for idx, entry in enumerate(scored_paths, start=1):
            # Label the top-scoring path explicitly
            rank_label = f"Path {idx}"
            if idx == 1:
                rank_label += " (Best)"

            fh.write(f"--- {rank_label} ---\n")

            # Format the node sequence with directional arrows
            node_seq_str = " → ".join(entry["path"])
            fh.write(f"  Node sequence : {node_seq_str}\n")

            # List individual edge weights for transparency
            weights_str = ", ".join(f"{w:.4f}" for w in entry["edge_weights"])
            fh.write(f"  Edge weights  : [{weights_str}]\n")

            # Final multiplicative path score
            fh.write(f"  Path Score    : {entry['path_score']:.6f}\n")
            fh.write("\n")


# ═══════════════════════════════════════════════════════════════════════
#  PHASE 4: SUB-NETWORK VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════
def visualize_path_subnetwork(
    graph: nx.DiGraph,
    scored_paths: list[dict],
    source: str,
    target: str,
    output_path: str,
) -> None:
    """Visualise only the edges/nodes that belong to shortest paths.

    A sub-graph is extracted from the full PPI network containing
    the nodes and directed edges that appear in the shortest paths.

    We use `nx.kamada_kawai_layout` which produces clean,
    force-directed positions that tend to spread nodes evenly and
    make the directional flow from source → target clearly visible.

    Args:
      graph (nx.DiGraph): The original PPI directed graph (used to pull edge attributes).
      scored_paths (list[dict]): Output of `find_all_shortest_paths`.
      source (str): Source protein ID (highlighted green).
      target (str): Target protein ID (highlighted red).
      output_path (str): File path for the saved PNG image.
    """

    # ── Build the sub-graph from path edges ──────────────────────────
    subgraph = nx.DiGraph()

    for entry in scored_paths:
        path_nodes = entry["path"]
        for i in range(len(path_nodes) - 1):
            tail = path_nodes[i]
            head = path_nodes[i + 1]
            # Copy edge data (weight) from the original graph
            edge_data = graph[tail][head]
            subgraph.add_edge(tail, head, **edge_data)

    # ── Determine node colours ───────────────────────────────────────
    # Source = green, target = red/salmon, intermediates = sky-blue.
    node_colors = []
    for node in subgraph.nodes():
        if node == source:
            node_colors.append("#2ecc71")   # emerald green
        elif node == target:
            node_colors.append("#e74c3c")   # alizarin red
        else:
            node_colors.append("#3498db")   # peter-river blue

    # ── Layout: Kamada-Kawai produces clean, force-directed positions ─
    pos = nx.kamada_kawai_layout(subgraph)

    # ── Edge labels: show the weight on each arrow ───────────────────
    edge_labels = {
        (u, v): f"{d['weight']:.2f}"
        for u, v, d in subgraph.edges(data=True)
    }

    # ── Draw ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_title(
        f"Shortest-Path Sub-Network\n{source}  \u2192  {target}",
        fontsize=14,
        fontweight="bold",
    )

    # Nodes
    nx.draw_networkx_nodes(
        subgraph, pos, ax=ax,
        node_size=1500,
        node_color=node_colors,
        edgecolors="#2c3e50",       # dark border for contrast
        linewidths=1.5,
    )

    # Node labels (protein IDs)
    nx.draw_networkx_labels(
        subgraph, pos, ax=ax,
        font_size=9,
        font_weight="bold",
        font_color="white",
    )

    # Directed edges with visible arrowheads
    nx.draw_networkx_edges(
        subgraph, pos, ax=ax,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=18,
        edge_color="#7f8c8d",
        width=2.0,
        connectionstyle="arc3,rad=0.1",   # slight curve avoids overlap
    )

    # Edge-weight labels
    nx.draw_networkx_edge_labels(
        subgraph, pos, ax=ax,
        edge_labels=edge_labels,
        font_size=8,
        font_color="#c0392b",
        label_pos=0.4,
    )

    ax.axis("off")
    fig.tight_layout()

    # ── Save to disk ─────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
#  PHASE 5: INTEGRATION — PUBLIC ENTRY POINTS
# ═══════════════════════════════════════════════════════════════════════
def analyze_shortest_paths(
    graph: nx.DiGraph,
    source_node: str,
    target_node: str,
    results_dir: str = "results",
) -> list[dict]:
    """End-to-end shortest-path analysis pipeline (Phases 1 → 4).

    Performs the complete analysis workflow:
      1. Find all shortest acyclic paths between `source_node`
         and `target_node` using Dijkstra's algorithm.
      2. Score each path by multiplying edge-weight probabilities.
      3. Save the results to a formatted text report.
      4. Visualize the sub-network of discovered paths.

    Args:
      graph (nx.DiGraph): The PPI directed graph with `weight` edge attributes.
      source_node (str): UniProt-style ID of the source protein.
      target_node (str): UniProt-style ID of the target protein.
      results_dir (str): Output directory (default `results`).

    Returns:
      list[dict]: The scored paths (same structure as `find_all_shortest_paths`).
    """

    # ── Phase 1 & 2: Pathfinding + Scoring ───────────────────────────
    print(f"Finding shortest paths: {source_node} -> {target_node} ...")
    scored_paths = find_all_shortest_paths(graph, source_node, target_node)
    print(f"  Found {len(scored_paths)} shortest path(s).")

    # ── Phase 3: Text report ─────────────────────────────────────────
    txt_path = os.path.join(results_dir, "1_shortest_paths.txt")
    save_paths_to_file(scored_paths, source_node, target_node, txt_path)
    print(f"  [OK] Report saved to {txt_path}")

    # ── Phase 4: Sub-network plot ────────────────────────────────────
    png_path = os.path.join(results_dir, "1_shortest_paths_subnetwork.png")
    visualize_path_subnetwork(
        graph, scored_paths, source_node, target_node, png_path,
    )
    print(f"  [OK] Sub-network plot saved to {png_path}")

    return scored_paths


def run_path_analysis(
    graph: nx.DiGraph,
    source_node: str,
    target_node: str,
    results_dir: str = "results",
) -> list[dict] | None:
    """Thin wrapper around `analyze_shortest_paths` for main.py.

    This function catches `NetworkXNoPath` / `NodeNotFound`
    and prints a user-friendly notice instead of propagating the exception.

    Args:
      graph (nx.DiGraph): The PPI directed graph.
      source_node (str): Source protein UniProt ID.
      target_node (str): Target protein UniProt ID.
      results_dir (str): Output directory (default `results`).

    Returns:
      list[dict] | None: Scored paths on success, or `None` if no path was found.
    """
    try:
        return analyze_shortest_paths(
            graph, source_node, target_node, results_dir,
        )
    except (nx.NetworkXNoPath, nx.NodeNotFound) as exc:
        print(f"  [Warning] Path analysis skipped -- {exc}")
        return None
