# Protein-Protein Interactions with NetworkX

This project analyzes a directed protein-protein interaction (PPI) network using Python and NetworkX. The dataset is the PathLinker 2018 human interactome (`PathLinker_2018_human-ppi-weighted-cap0_75.txt`). Each line in the file describes one directed interaction with a source protein, target protein, confidence score, and annotation method.

## What the project does

- Builds a directed PPI graph from the PathLinker interactome.
- Converts the graph to an unweighted adjacency matrix and saves it.
- Computes node degree statistics and produces a ranked protein list and degree histogram.
- Finds shortest acyclic paths between selected proteins, scores them, and saves a text report plus subnetwork figure.
- Maps UniProt protein IDs to gene names using UniProt via `bioservices`.

## Repository structure

```text
├── data
│   └── PathLinker_2018_human-ppi-weighted-cap0_75.txt
├── README.md
├── requirements.txt
├── results
│   ├── 2_single_protein_degree.txt
│   ├── 3_multiprotein_histogram.png
│   ├── 3_multiprotein_ranked.txt
│   └── 4_adjacency_matrix.csv
├── src
│   ├── __init__.py
│   ├── __pycache__
│   │   ├── __init__.cpython-310.pyc
│   │   ├── graph_builder.cpython-310.pyc
│   │   └── graph_builder.cpython-312.pyc
│   ├── connectivity_analysis.py
│   ├── graph_builder.py
│   ├── id_converter.py
│   ├── main.py
│   ├── path_analysis.py
│   └── test
│       ├── __init__.py
│       └── test_path_analysis.py
└── test_connectivity.py
```

## Requirements

Install the Python dependencies listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

## How to run

Run the full pipeline from the project root:

```bash
python3 src/main.py
```

This will:

1. Parse a sample of the PathLinker network.
2. Build a directed NetworkX graph `nx.DiGraph()`.
3. Save the adjacency matrix to `results/4_adjacency_matrix.csv`.
4. Generate connectivity outputs in `results/`.
5. Run shortest-path analysis if the chosen source and target proteins exist in the graph.
6. Query UniProt for gene-name mapping.

## Output files

- `results/4_adjacency_matrix.csv`
- `results/2_single_protein_degree.txt`
- `results/3_multiprotein_ranked.txt`
- `results/3_multiprotein_histogram.png`
- `results/1_shortest_paths.txt`
- `results/1_shortest_paths_subnetwork.png`

## Notes

- The PathLinker file is directed: each row defines an interaction from tail to head.
- Edge confidence values are stored as `weight` and converted to `distance = -log(weight)` for shortest-path calculations.
- UniProt lookups require network access.

## Team Members

- Saif Mohammed
- Mohammed Ashraf
- Mostafa Ashraf
- Mazen Atef
