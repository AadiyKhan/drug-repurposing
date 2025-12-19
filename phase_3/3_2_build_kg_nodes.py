import pandas as pd
import networkx as nx
import pickle
from pathlib import Path

# --- Configuration ---
INPUT_FILE = 'drugbank_final_features_with_nlp.csv'
ENRICHMENT_FILE = 'disease_target_enrichment.csv'
GRAPH_OUTPUT = 'drug_repurposing_graph.pkl'

def build_graph_nodes():
    if not Path(INPUT_FILE).exists():
        print(f"Error: {INPUT_FILE} not found.")
        return

    print(f"Loading drug data from {INPUT_FILE}...")
    df_drugs = pd.read_csv(INPUT_FILE)
    
    print(f"Loading enrichment data from {ENRICHMENT_FILE}...")
    df_diseases = pd.read_csv(ENRICHMENT_FILE)
    
    # Initialize the Graph
    G = nx.Graph()
    
    # 1. Add Drug Nodes
    print("Adding drug nodes...")
    for _, row in df_drugs.iterrows():
        # We store the Name as an attribute for visualization/referencing later
        G.add_node(row['DrugBank ID'], 
                   type='drug', 
                   name=row['Name'])
    
    # 2. Add Target Nodes (Proteins)
    print("Adding target nodes...")
    all_targets = set()
    # Split the '|' separated UniProt IDs and add to a set
    df_drugs['UniProt_IDs'].str.replace(' ', '').str.split('|').dropna().apply(all_targets.update)
    # Remove the placeholder if it exists
    all_targets.discard('NO_TARGETS_LISTED')
    
    for target in all_targets:
        G.add_node(target, type='target')
        
    # 3. Add Disease Nodes
    print("Adding disease nodes...")
    for _, row in df_diseases.iterrows():
        G.add_node(row['Disease_Name'], type='disease')
        
    # Summary
    print("\n✅ Knowledge Graph Nodes Initialized.")
    print(f"Total Nodes: {G.number_of_nodes()}")
    print(f"- Drug Nodes: {len(df_drugs)}")
    print(f"- Target Nodes: {len(all_targets)}")
    print(f"- Disease Nodes: {len(df_diseases)}")
    
    # Save the graph object to disk
    with open(GRAPH_OUTPUT, 'wb') as f:
        pickle.dump(G, f)
    print(f"\nGraph structure saved to {GRAPH_OUTPUT}")

if __name__ == "__main__":
    build_graph_nodes()