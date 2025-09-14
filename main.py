import sys
import os
import glob
import pandas as pd
from pygosemsim import similarity, graph

def precision(pred: set, ref: set) -> float:
    return len(pred & ref) / len(pred)

def recall(pred: set, ref: set) -> float:
    return len(pred & ref) / len(ref)

def f1_score(precision: float, recall: float) -> float:
    if precision == 0 or recall == 0:
        return 0
    return (2 * precision * recall) / (precision + recall)


def calculate_go_f1_scores(cluster_genes, go_terms_dict):
    """Calculate F1 scores between a cluster and all GO terms"""
    go_f1_scores = {}
    
    for go_term, go_genes in go_terms_dict.items():
        # Calculate precision and recall
        prec = precision(cluster_genes, go_genes)
        rec = recall(cluster_genes, go_genes)
        f1 = f1_score(prec, rec)
        
        if f1 > 0:  # Only include GO terms with non-zero F1 score
            go_f1_scores[go_term] = f1
    
    return go_f1_scores


def get_best_go_term(cluster_genes, go_terms_dict):
    """Get the GO term with the highest F1 score for a cluster"""
    go_f1_scores = calculate_go_f1_scores(cluster_genes, go_terms_dict)
    
    if not go_f1_scores:
        return None, 0.0
    
    # Find the GO term with highest F1 score
    best_go_term = max(go_f1_scores.items(), key=lambda x: x[1])
    return best_go_term[0], best_go_term[1]


G = graph.from_resource("go-basic")
def calculate_go_similarity(go_term1, go_term2):
    """Calculate GO term similarity between two GO terms using Wang similarity"""
    if not go_term1 or not go_term2:
        return 0.0
    try:
        return similarity.wang(G, go_term1, go_term2)
    except:
        return 0.0

QUICK_GO = "PCC.tsv"

quick_go = pd.read_csv(QUICK_GO, sep="\t")
quick_go.columns = ["GENE PRODUCT DB", "GENE PRODUCT ID", "SYMBOL", "QUALIFIER", "GO TERM", "GO NAME", "ECO ID", "GO EVIDENCE CODE", "REFERENCE", "WITH/FROM", "TAXON ID", "ASSIGNED BY", "ANNOTATION EXTENSION", "GO ASPECT"]

sgd_features = pd.read_csv("SGD_features.tab", sep="\t")
sgd_features.columns = ["SGDID", "Feature_type", "Qualifier", "Systematic_name", "Gene_name", "Alias", "Parent_feature", "Secondary_SGDID", "Chromosome", "Start", "End", "Strand", "Genetic_position", "Coordinate_version", "Sequence_version", "Description"]

go_terms = quick_go.groupby("GO TERM")
go_terms_dict = {}
for key, group_df in go_terms:
    gene_names = group_df["SYMBOL"].unique()
    orf_names = set()
    for gene_name in gene_names:
        if gene_name in sgd_features["Gene_name"].values:
            orf_names.add(sgd_features[sgd_features["Gene_name"] == gene_name]["Systematic_name"].values[0])
    if len(orf_names) > 0:
        go_terms_dict[key] = orf_names

import glob
import sys

NETWORK_NAME = sys.argv[1]
THRESHOLD = 0.25

cluster_files = glob.glob(f"{NETWORK_NAME}/*clusters.txt")
# Dictionary to store cluster data
clusters_dict = {}
# Read each cluster file
for cluster_file in cluster_files:
    pred_clusters = set()
    with open(cluster_file) as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
            # Split line into cluster ID and gene list
            pred_cluster = line.strip().split()
            # Convert gene list to set of genes (as strings, not frozenset)
            pred_cluster = set(pred_cluster)
            pred_clusters.add(frozenset(pred_cluster))
        clusters_dict[cluster_file] = pred_clusters


# Find all cluster files and the first reference file

cluster_files = sorted(glob.glob(f"{NETWORK_NAME}/*_clusters.txt"))
reference_files = sorted(glob.glob(f"{NETWORK_NAME}/*_complexes.txt"))

if not cluster_files:
    print("No cluster files found.")
    sys.exit(1)
print(f"Found {len(cluster_files)} prediction files")
if not reference_files:
    print("No reference file found.")
    sys.exit(1)
REFERENCE_FILE = reference_files[0]
print(REFERENCE_FILE)
references = set()
reference_indices = []
with open(REFERENCE_FILE) as f:
    for idx, line in enumerate(f):
        ref = frozenset(line.strip().split())
        if len(ref) >= 2 and ref not in references:
            references.add(ref)
            reference_indices.append(idx)

references = list(references)

# Convert reference indices to string format once, outside the loop
reference_indices = [f"ref_{idx}" for idx in reference_indices]

def matching_score(set1, set2):
    """Calculates the matching score between two sets (e.g., a cluster and a complex)
    using the approach of Bader et al, 2001"""
    return (len(set1.intersection(set2))**2 / (float(len(set1)) * len(set2)))


for NETWORK_FILE in cluster_files:
    print(f"Processing {NETWORK_FILE}")
    predictions = set()

    # Read clusters from file and keep track of original line indices

    clusters = []
    cluster_indices = []
    low = []
    with open(NETWORK_FILE) as f:
        for idx, l in enumerate(f.read().splitlines()):
            pred = frozenset(l.strip().split())
            if len(pred) >= 3 and pred not in clusters:
                clusters.append(pred)
                cluster_indices.append(idx)
            elif len(pred) <= 2 or pred in clusters:
                low.append(pred)
                

    # Use references list and their original line indices
    references_list = references

    n = len(clusters)
    print(f"Number of clusters: {n}")
    print(f"Number of clusters <= 2: {len(low)}")
    print(f"Number of references: {len(references_list)}")
    
    
    # Calculate matching scores between clusters
    matching_matrix = [[None for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            else:
                matching_matrix[i][j] = matching_score(clusters[i], clusters[j])

    edgelist = []
    for i in range(n):
        for j in range(i + 1, n):  # Only consider upper triangle (i < j)
            if matching_matrix[i][j] is not None and matching_matrix[i][j] >= THRESHOLD:
                cluster_i = " ".join(sorted(clusters[i]))
                cluster_j = " ".join(sorted(clusters[j]))
                idx_i = cluster_indices[i]
                idx_j = cluster_indices[j]
                edgelist.append((idx_i, idx_j, matching_matrix[i][j], None))

    # Add edges between clusters and references
    for i, cluster in enumerate(clusters):
        for j, ref in enumerate(references_list):
            score = matching_score(cluster, ref)
            if score >= THRESHOLD:
                idx_i = cluster_indices[i]
                idx_j = reference_indices[j]
                # Find all matching references for this cluster
                matching_refs = []
                for k, other_ref in enumerate(references_list):
                    other_score = matching_score(cluster, other_ref)
                    if other_score >= THRESHOLD:
                        ref_proteins = " ".join(sorted(other_ref))
                        matching_refs.append((ref_proteins, other_score))
                edgelist.append((idx_i, idx_j, score, matching_refs))

    # Add edges between references
    for i in range(len(references_list)):
        for j in range(i + 1, len(references_list)):
            score = matching_score(references_list[i], references_list[j])
            if score >= THRESHOLD:
                idx_i = reference_indices[i]
                idx_j = reference_indices[j]
                edgelist.append((idx_i, idx_j, score, None))

    # GDF
    network_base = os.path.basename(NETWORK_FILE)
    reference_base = os.path.basename(REFERENCE_FILE)
    output_gdf = f"{network_base}-{reference_base}-network.gdf"
    
    # Map from file index to cluster list index
    fileidx_to_clusteridx = {file_idx: i for i, file_idx in enumerate(cluster_indices)}
    # Map from reference index to reference list index
    refidx_to_refidx = {ref_idx: i for i, ref_idx in enumerate(reference_indices)}

    # Prepare per-node CC/BP/MF values when using CSV
    lineindex_to_vals = {}

    with open(output_gdf, "w") as f:
        f.write("nodedef>name VARCHAR,label VARCHAR,size INT,match_reference INT,match_references INT,type VARCHAR, max_OS DOUBLE, max_OS_ref DOUBLE, matching_ref_ids VARCHAR, best_go_term VARCHAR, best_go_f1 DOUBLE, CC_Value DOUBLE, BP_Value DOUBLE, MF_Value DOUBLE\n")
        
        # Write cluster nodes (including isolated ones)
        for idx, cluster in zip(cluster_indices, clusters):
            label = " ".join(sorted(cluster))
            size = len(cluster)
            match_reference = False
            matching_ref_ids = []
            for i, ref in enumerate(references_list):
                if matching_score(cluster, ref) >= THRESHOLD:
                    match_reference = True
                    if reference_indices[i] not in matching_ref_ids:
                        matching_ref_ids.append(reference_indices[i])
            match_references = len(matching_ref_ids)
            
            # Calculate max OS from neighbors (other clusters that have edges)
            max_OS = 0
            cluster_idx = cluster_indices.index(idx)
            has_connected_neighbors = False
            for j, other_cluster in enumerate(clusters):
                if j != cluster_idx:  # Don't compare with itself
                    os_score = matching_score(cluster, other_cluster)
                    if os_score >= THRESHOLD:
                        has_connected_neighbors = True
                        if os_score > max_OS:
                            max_OS = os_score
            
            # If no connected neighbors, max_OS remains 0
            
            # Get best GO term for this cluster
            best_go_term, best_go_f1 = get_best_go_term(cluster, go_terms_dict)
            best_go_term_str = best_go_term if best_go_term else ""
            

            if match_references > 0:
                # Calculate max OS from connected references only
                max_OS_ref = 0
                for ref in references_list:
                    os_score = matching_score(cluster, ref)
                    if os_score >= THRESHOLD and os_score > max_OS_ref:
                        max_OS_ref = os_score
                matching_refs_str = ";".join(matching_ref_ids)
                f.write(f"{idx},'{label}',{size},{int(match_reference)},{match_references},'matched_prediction',{max_OS:.3f},{max_OS_ref:.3f},'{matching_refs_str}','{best_go_term_str}',{best_go_f1:.3f}\n")
            else:
                max_OS_ref = 0
                matching_refs_str = ";".join(matching_ref_ids)
                f.write(f"{idx},'{label}',{size},{int(match_reference)},{match_references},'prediction',{max_OS:.3f},{max_OS_ref:.3f},'{matching_refs_str}','{best_go_term_str}',{best_go_f1:.3f}\n")

        # Write reference nodes (including isolated ones)
        for idx, ref in zip(reference_indices, references_list):
            label = " ".join(sorted(ref))
            size = len(ref)
            match_reference = True  # References always match themselves
            match_references = 1    # References match themselves
            
            # Check if this reference matches any clusters
            matches_any_cluster = False
            for cluster in clusters:
                if matching_score(cluster, ref) >= THRESHOLD:
                    matches_any_cluster = True
                    break
            
            # Determine reference type based on whether it matches any clusters
            if matches_any_cluster:
                ref_type = 'matched_reference'
                max_OS_ref = 1.0
            else:
                ref_type = 'reference'
                max_OS_ref = 1.0
            
            # Get best GO term for this reference
            best_go_term, best_go_f1 = get_best_go_term(ref, go_terms_dict)
            best_go_term_str = best_go_term if best_go_term else ""

            # Reference nodes: set one-hot flags from reference aspect
            f.write(f"{idx},'{label}',{size},{int(match_reference)},{match_references},'{ref_type}',1,{max_OS_ref:.3f},'','{best_go_term_str}',{best_go_f1:.3f},{0},{0},{0}\n")

        f.write("edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE,common INT,matching_refs VARCHAR,wang_similarity DOUBLE\n")
        for edge in edgelist:
            i, j, weight, ref_info = edge
            
            # Initialize GO similarity variables
            wang_similarity = 0.0
            ref_info_str = ""
            # Handle cluster-cluster edges
            if i in fileidx_to_clusteridx and j in fileidx_to_clusteridx:
                ci = fileidx_to_clusteridx[i]
                cj = fileidx_to_clusteridx[j]
                common = len(clusters[ci].intersection(clusters[cj]))
                # Get best GO terms for both clusters
                best_go1, _ = get_best_go_term(clusters[ci], go_terms_dict)
                best_go2, _ = get_best_go_term(clusters[cj], go_terms_dict)
                # Calculate Wang similarity between best GO terms
                wang_similarity = calculate_go_similarity(best_go1, best_go2)
            
            # Handle cluster-reference edges
            elif i in fileidx_to_clusteridx and j in refidx_to_refidx:
                ci = fileidx_to_clusteridx[i]
                rj = refidx_to_refidx[j]
                common = len(clusters[ci].intersection(references_list[rj]))
                best_go1, _ = get_best_go_term(clusters[ci], go_terms_dict)
                best_go2, _ = get_best_go_term(references_list[rj], go_terms_dict)
                wang_similarity = calculate_go_similarity(best_go1, best_go2)
                # Format matching references info for cluster-reference edges
                if ref_info is not None:
                    ref_info_str = ";".join([f"'{ref_proteins}'({score:.3f})" for ref_proteins, score in ref_info])
                else:
                    ref_info_str = ""
            elif i in refidx_to_refidx and j in fileidx_to_clusteridx:
                ri = refidx_to_refidx[i]
                cj = fileidx_to_clusteridx[j]
                common = len(references_list[ri].intersection(clusters[cj]))
                ref_info_str = ""
                best_go1, _ = get_best_go_term(references_list[ri], go_terms_dict)
                best_go2, _ = get_best_go_term(clusters[cj], go_terms_dict)
                wang_similarity = calculate_go_similarity(best_go1, best_go2)
            # Handle reference-reference edges
            elif i in refidx_to_refidx and j in refidx_to_refidx:
                ri = refidx_to_refidx[i]
                rj = refidx_to_refidx[j]
                common = len(references_list[ri].intersection(references_list[rj]))
                ref_info_str = ""
                best_go1, _ = get_best_go_term(references_list[ri], go_terms_dict)
                best_go2, _ = get_best_go_term(references_list[rj], go_terms_dict)
                wang_similarity = calculate_go_similarity(best_go1, best_go2)
            else:
                common = 0
                ref_info_str = ""
                
            f.write(f"{i},{j},{weight},{common},'{ref_info_str}',{wang_similarity:.3f}\n")
