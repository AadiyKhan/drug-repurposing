import pandas as pd
import networkx as nx
import pickle
from pathlib import Path

# --- Configuration ---
GRAPH_INPUT = 'drug_repurposing_graph_final.pkl'
DRUG_DATA = 'drugbank_final_features_with_nlp.csv'
OUTPUT_REPORT = 'repurposing_candidates_final.csv'

def predict_repurposing_candidates():
    if not Path(GRAPH_INPUT).exists():
        print(f"Error: {GRAPH_INPUT} not found.")
        return

    # Load the graph and drug data
    print("Loading Knowledge Graph...")
    with open(GRAPH_INPUT, 'rb') as f:
        G = pickle.load(f)
    
    df_drugs = pd.read_csv(DRUG_DATA)
    drug_names = dict(zip(df_drugs['DrugBank ID'], df_drugs['Name']))
    
    # Identify node types
    disease_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'disease']
    drug_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'drug']

    results = []

    print(f"Scanning for candidates across {len(disease_nodes)} diseases...")

    for disease in disease_nodes:
        print(f"  > Processing: {disease}")
        
        # 1. Find 'Direct' Drugs (already linked via a target)
        # Path: Drug -> Target -> Disease
        direct_drugs = set()
        for target in G.neighbors(disease):
            for drug in G.neighbors(target):
                if G.nodes[drug].get('type') == 'drug':
                    direct_drugs.add(drug)

        # 2. Find 'Repurposing' Candidates
        # Path: New Drug --(Similarity)--> Direct Drug --(Target)--> Disease
        candidates = {}
        
        for d_drug in direct_drugs:
            for neighbor in G.neighbors(d_drug):
                # Check if the neighbor is a drug and NOT already a direct drug
                if G.nodes[neighbor].get('type') == 'drug' and neighbor not in direct_drugs:
                    
                    # Calculate a simple score based on edge weights
                    edge_data = G.get_edge_data(d_drug, neighbor)
                    weight = edge_data.get('weight', 0.5) # Default weight if not similarity edge
                    
                    if neighbor not in candidates:
                        candidates[neighbor] = 0
                    candidates[neighbor] += weight

        # 3. Sort and save top 5 for this disease
        sorted_candidates = sorted(candidates.items(), key=lambda x: x[1], reverse=True)[:5]
        
        for cand_id, score in sorted_candidates:
            results.append({
                'Disease': disease,
                'Candidate_DrugBank_ID': cand_id,
                'Candidate_Name': drug_names.get(cand_id, "Unknown"),
                'Repurposing_Score': round(score, 4),
                'Num_Known_Similar_Drugs': len([d for d in direct_drugs if G.has_edge(cand_id, d)])
            })

    # Save to CSV
    final_df = pd.DataFrame(results)
    final_df.to_csv(OUTPUT_REPORT, index=False)
    
    print(f"\n✅ Prediction Engine Complete!")
    print(f"Top candidates saved to: {OUTPUT_REPORT}")
    
    # Display results
    print("\n--- SAMPLE REPURPOSING PREDICTIONS ---")
    if not final_df.empty:
        print(final_df.sort_values(by='Repurposing_Score', ascending=False).head(15).to_markdown(index=False))
    else:
        print("No candidates found with current thresholds.")

if __name__ == "__main__":
    predict_repurposing_candidates()