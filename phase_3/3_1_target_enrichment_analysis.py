import pandas as pd
from pathlib import Path

# --- Configuration ---
INPUT_DRUG_FILE = 'drugbank_final_features_with_nlp.csv'
OUTPUT_ENRICHMENT_FILE = 'disease_target_enrichment.csv'
DISGENET_MOCK_FILE = 'disgenet_gda_subset.csv'

# --- 1. Create Mock DisGeNET GDA Dataset ---
# In a real scenario, this would be a downloaded and pre-processed file from DisGeNET
def create_mock_disgenet():
    """Creates a small mock file to simulate the Gene-Disease Association data."""
    
    # We use common UniProt IDs and associated diseases for a realistic mock-up
    data = {
        'UniProt_ID': [
            'P00533', 'P00533', 'P04626', 'P04626', 'P04626', 'P08630', 'P08630', 
            'P10600', 'P10600', 'Q13065', 'Q13065', 'Q13065', 'P14416', 'P14416', 
            'P04150', 'P04150', 'P04150', 'P08630'
        ],
        'Disease_Name': [
            'Breast Cancer', 'Non-Small Cell Lung Cancer', 'Alzheimer Disease', 
            'Depressive Disorder', 'Schizophrenia', 'Rheumatoid Arthritis', 
            'Systemic Lupus Erythematosus', 'Type 2 Diabetes Mellitus', 
            'Obesity', 'Hypertension', 'Stroke', 'Myocardial Infarction',
            'HIV Infections', 'AIDS', 'Primary Immunodeficiency', 'Transplant Rejection',
            'Severe Combined Immunodeficiency', 'Rheumatoid Arthritis' # Duplicated to increase count
        ],
        'Score': [0.95, 0.88, 0.99, 0.75, 0.70, 0.90, 0.85, 0.92, 0.80, 0.93, 0.87, 0.82, 0.95, 0.90, 0.99, 0.98, 0.97, 0.91]
    }
    
    df_gda = pd.DataFrame(data)
    df_gda.to_csv(DISGENET_MOCK_FILE, index=False)
    print(f"Created mock DisGeNET GDA file: {DISGENET_MOCK_FILE}")
    return df_gda

# --- 2. Target Enrichment Analysis ---

def run_target_enrichment_analysis():
    
    # Check for drug input file
    if not Path(INPUT_DRUG_FILE).exists():
        print(f"Error: Required drug features file '{INPUT_DRUG_FILE}' not found.")
        return

    # Create and load the mock GDA file
    df_gda = create_mock_disgenet()
    
    print(f"Loading final drug feature data from {INPUT_DRUG_FILE}...")
    df_drugs = pd.read_csv(INPUT_DRUG_FILE)
    
    # 3. Process Drug Targets
    # Explode the pipe-separated UniProt_IDs column into separate rows
    df_targets = df_drugs[['DrugBank ID', 'UniProt_IDs']].copy()
    
    # Clean up the ID column before splitting (remove extra spaces)
    df_targets['UniProt_IDs'] = df_targets['UniProt_IDs'].str.replace(' ', '')
    
    # Split the '|' separated IDs and stack them
    df_targets_expanded = (df_targets.assign(UniProt_ID=df_targets['UniProt_IDs'].str.split('|'))
                          .explode('UniProt_ID')
                          .drop(columns=['UniProt_IDs'])
                          .dropna(subset=['UniProt_ID'])
                          .reset_index(drop=True)
                         )
    
    # Remove rows where UniProt_ID is a placeholder like 'NO_TARGETS_LISTED'
    df_targets_expanded = df_targets_expanded[~df_targets_expanded['UniProt_ID'].str.contains('NO_TARGETS', na=False)]
    
    print(f"Total unique drug-target associations found: {len(df_targets_expanded)}")
    
    # 4. Merge Drug Targets with DisGeNET Diseases
    # Join the expanded drug targets (DrugBank ID and UniProt_ID) with the DisGeNET GDA data
    df_merged = pd.merge(
        df_targets_expanded,
        df_gda,
        on='UniProt_ID',
        how='inner' # Only keep drugs whose targets are associated with a disease in DisGeNET
    )
    
    print(f"Total drug-disease associations after merging with DisGeNET: {len(df_merged)}")

    # 5. Enrichment Analysis: Count Unique Drugs per Disease
    # We count how many UNIQUE DrugBank IDs are associated with each Disease_Name
    # This gives us the "enrichment score" or the number of known drugs that could target
    # the disease's underlying genes.
    df_enrichment = (df_merged
                     .groupby('Disease_Name')['DrugBank ID']
                     .nunique()
                     .sort_values(ascending=False)
                     .reset_index()
                     .rename(columns={'DrugBank ID': 'Num_Associated_Drugs'})
                    )
    
    # 6. Save and Report Top 10
    top_10_diseases = df_enrichment.head(10)
    top_10_diseases.to_csv(OUTPUT_ENRICHMENT_FILE, index=False)

    print(f"\n✅ Target Enrichment Analysis Complete. Top 10 diseases saved to {OUTPUT_ENRICHMENT_FILE}")
    
    print("\n--- Top 10 Most Enriched Diseases (High-Priority Candidates) ---")
    print(top_10_diseases.to_markdown(index=False, numalign="left", stralign="left"))
    
    print("\nThis analysis suggests the best diseases to focus on for repurposing.")


if __name__ == "__main__":
    # Ensure you have 'drugbank_final_features_with_nlp.csv' in the directory!
    run_target_enrichment_analysis()