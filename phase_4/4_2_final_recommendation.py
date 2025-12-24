import pandas as pd
import networkx as nx
import pickle

# --- Configuration ---
GRAPH_INPUT = 'drug_repurposing_graph_final.pkl'
CANDIDATE_FILE = 'repurposing_candidates_final.csv' # or 'repurposing_candidates_filtered.csv'

def generate_evidence_report():
    print("🔬 Extracting Biological Evidence for Top Candidates...")
    
    # Load Data
    with open(GRAPH_INPUT, 'rb') as f:
        G = pickle.load(f)
    df_cand = pd.read_csv(CANDIDATE_FILE)
    
    # Get top 3 unique candidates based on Repurposing Score
    top_hits = df_cand.sort_values(by='Repurposing_Score', ascending=False).head(3)
    
    for _, row in top_hits.iterrows():
        drug = row['Candidate_DrugBank_ID']
        drug_name = row['Candidate_Name']
        disease = row['Disease']
        
        print(f"\n" + "="*50)
        print(f"REPORT FOR: {drug_name} ({drug})")
        print(f"PREDICTED USE: {disease}")
        print("-" * 50)
        
        # 1. Find similar drugs that are already used for this disease
        direct_drugs = []
        for target in G.neighbors(disease):
            for neighbor in G.neighbors(target):
                if G.nodes[neighbor].get('type') == 'drug':
                    direct_drugs.append(neighbor)
        
        similar_known_drugs = [n for n in G.neighbors(drug) if n in direct_drugs]
        
        print(f"💡 AI Logic: This drug is similar to {len(similar_known_drugs)} drugs already linked to {disease}.")
        print(f"Top Similarity Peers: {', '.join([G.nodes[d]['name'] for d in similar_known_drugs[:5]])}")
        
        # 2. Extract Shared Biological Targets
        # Path: Candidate -> Peer -> Target -> Disease
        shared_targets = set()
        for peer in similar_known_drugs:
            peer_targets = [t for t in G.neighbors(peer) if G.nodes[t].get('type') == 'target']
            disease_targets = [t for t in G.neighbors(disease) if G.nodes[t].get('type') == 'target']
            
            # Common targets between the peer drug and the disease
            overlap = set(peer_targets).intersection(set(disease_targets))
            shared_targets.update(overlap)
            
        print(f"\n🧬 Key Biological Bridges (Targets to investigate):")
        if shared_targets:
            for t in list(shared_targets)[:10]:
                print(f" - {t}")
        else:
            print(" - No direct protein bridge found; similarity is likely purely chemical/NLP-based.")
        
        print("="*50)

if __name__ == "__main__":
    generate_evidence_report()