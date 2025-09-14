
# A Network Perspective on Protein Complex Predictions

This repository contains the code and resources used in the 2025 paper **"A Network Perspective on Protein Complex Predictions"**. 

It includes implementations for generating and analyzing prediction networks of protein complexes, evaluation scripts, and data processing tools. The repository is intended to support reproducibility and further exploration of network-based protein complex prediction methods.

The project includes a main script `main.py` that takes two parameters. The first parameter is a folder containing prediction and reference files: predictions must end with `clusters.txt` and references must end with `complexes.txt`. The second parameter is the overlap score threshold used for analysis.

## Contents

- `main.py` – Main script.  
- `PCC.tsv` – GO terms exported as TSV from QuickGO.  
- `SGD_features.tab` – GO protein annotations to GO terms.



## Description of node and edge attributes in the GDF output

| **Attribute**       | **Description**                                                   |
|----------------------|-------------------------------------------------------------------|
| **Node attributes**  |                                                                   |
| name                | Unique identifier of the node                                     |
| label               | Proteins in the prediction or reference separated by space        |
| size                | Number of proteins                                                |
| match_reference     | Indicator if the node matches at least one reference (0/1)        |
| match_references    | Number of reference nodes matched                                 |
| type                | Node type ((matched) prediction / reference)                      |
| max_OS              | Maximum overlap score with any neighbor                           |
| max_OS_ref          | Maximum overlap score with any reference                          |
| matching_ref_ids    | IDs of matched reference nodes                                    |
| best_go_term        | GO term with the highest F1 score                                 |
| best_go_f1          | F1 score with the best-matching GO term                           |
| **Edge attributes** |                                                                   |
| node1               | Source node of the edge                                           |
| node2               | Target node of the edge                                           |
| weight              | Edge weight (e.g., overlap score)                                 |
| common              | Number of common proteins between connected nodes                 |
| wang_similarity     | Wang semantic similarity between GO terms of connected nodes      |
