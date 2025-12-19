import pandas as pd
import networkx as nx
import pickle
import numpy as np
from pathlib import Path
import time

# --- Configuration ---
INPUT_FILE = 'drugbank_final_features_with_nlp.csv'
GDA_FILE = 'disgenet_gda_subset.csv'
GRAPH_INPUT = 'drug_repurposing_graph.pkl'
GRAPH_OUTPUT = 'drug_repurposing_graph_final.pkl'

# Thresholds
TANIMOTO_THRESHOLD = 0.7
COSINE_THRESHOLD = 0.85

def build_kg_edges_optimized():
    start_time = time.time()
    
    if not Path(GRAPH_INPUT).exists():
        print(f"Error: {GRAPH_INPUT} not found. Run 3.2 first.")
        return

    # Load data
    with open(GRAPH_INPUT, 'rb') as f:
        G = pickle.load(f)
    
    print("Loading dataframes...")
    df = pd.read_csv(INPUT_FILE)
    df_gda = pd.read_csv(GDA_FILE)
    drug_ids = df['DrugBank ID'].tolist()

    # 1. Known Biological Edges (Fast)
    print("1/4: Adding Drug-Target Edges...")
    dt_edges = []
    for _, row in df.iterrows():
        d_id = row['DrugBank ID']
        targets = str(row['UniProt_IDs']).replace(' ', '').split('|')
        for t in targets:
            if t in G:
                dt_edges.append((d_id, t, {'label': 'targets'}))
    G.add_edges_from(dt_edges)

    print("2/4: Adding Target-Disease Edges...")
    td_edges = []
    for _, row in df_gda.iterrows():
        if row['UniProt_ID'] in G and row['Disease_Name'] in G:
            td_edges.append((row['UniProt_ID'], row['Disease_Name'], 
                             {'label': 'associated_with', 'score': row['Score']}))
    G.add_edges_from(td_edges)

    # 2. Chemical Similarity (Optimized)
    print("3/4: Calculating Chemical Similarity (Tanimoto)...")
    fp_cols = [c for c in df.columns if c.startswith('FP_')]
    fp_matrix = df[fp_cols].values.astype(np.float32)
    
    dot_product = np.dot(fp_matrix, fp_matrix.T)
    sums = fp_matrix.sum(axis=1)
    tanimoto = dot_product / (sums[:, None] + sums[None, :] - dot_product + 1e-10)
    
    # Free up memory
    del fp_matrix
    del dot_product
    
    print("   - Filtering and adding chemical edges...")
    rows, cols = np.where(tanimoto > TANIMOTO_THRESHOLD)
    chem_edges = []
    for i, j in zip(rows, cols):
        if i < j:
            chem_edges.append((drug_ids[i], drug_ids[j], 
                               {'label': 'chem_similarity', 'weight': float(tanimoto[i, j])}))
    
    G.add_edges_from(chem_edges)
    del tanimoto # Clear matrix
    print(f"   - Added {len(chem_edges)} chemical similarity edges.")

    # 3. Mechanism Similarity (Optimized)
    print("4/4: Calculating Mechanism Similarity (Cosine)...")
    nlp_cols = [c for c in df.columns if c.startswith('NLP_Vec_')]
    nlp_matrix = df[nlp_cols].values.astype(np.float32)
    
    # Normalized dot product is cosine similarity
    norms = np.linalg.norm(nlp_matrix, axis=1)
    nlp_norm = nlp_matrix / (norms[:, None] + 1e-10)
    cosine = np.dot(nlp_norm, nlp_norm.T)
    
    print("   - Filtering and adding mechanism edges...")
    rows, cols = np.where(cosine > COSINE_THRESHOLD)
    mech_edges = []
    for i, j in zip(rows, cols):
        if i < j:
            mech_edges.append((drug_ids[i], drug_ids[j], 
                               {'label': 'mech_similarity', 'weight': float(cosine[i, j])}))
    
    G.add_edges_from(mech_edges)
    print(f"   - Added {len(mech_edges)} mechanism similarity edges.")

    # Save
    print(f"\n✅ Knowledge Graph Complete!")
    print(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")
    print(f"Total time: {time.time() - start_time:.2f} seconds")
    
    with open(GRAPH_OUTPUT, 'wb') as f:
        pickle.dump(G, f)

if __name__ == "__main__":
    build_kg_edges_optimized()