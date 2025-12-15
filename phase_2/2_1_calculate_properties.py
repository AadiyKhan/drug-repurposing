import pandas as pd
from pathlib import Path

# Try to import RDKit modules
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("🚨 ERROR: RDKit not found.")
    print("Please run: pip install rdkit")
    exit()

# --- Configuration ---
INPUT_FILE = 'drugbank_with_smiles.csv'
OUTPUT_FILE = 'drugbank_small_molecules_features.csv'

# --- RDKit Calculation Function ---

def calculate_rdkit_properties(smiles):
    """Calculates five key molecular properties for a given SMILES string."""
    # Attempt to create an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        # If RDKit fails to parse the SMILES, return None for all properties
        return pd.Series([None, None, None, None, None], index=['MW', 'AlogP', 'HBD', 'HBA', 'TPSA'])
    
    # Calculate the desired descriptors
    try:
        mw = Descriptors.MolWt(mol)
        alogp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        
        return pd.Series([mw, alogp, hbd, hba, tpsa], index=['MW', 'AlogP', 'HBD', 'HBA', 'TPSA'])
    except:
        # Catch unexpected errors during calculation (e.g. invalid atom valence)
        return pd.Series([None, None, None, None, None], index=['MW', 'AlogP', 'HBD', 'HBA', 'TPSA'])


# --- Main Processing Script ---

def run_feature_calculation():
    
    if not Path(INPUT_FILE).exists():
        print(f"Error: Required input file '{INPUT_FILE}' not found.")
        return

    print(f"Loading data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE)
    
    # 1. Initial Filtering: Remove rows where SMILES_String is null (removes primary biologics)
    initial_count = len(df)
    df_filtered = df.dropna(subset=['SMILES_String']).copy()
    
    # 2. Apply the RDKit property calculation function
    print(f"Starting RDKit calculation for {len(df_filtered)} drugs...")
    
    # The .apply method is fast for this type of operation
    new_features = df_filtered['SMILES_String'].apply(calculate_rdkit_properties)
    
    # Merge the new feature columns into the main DataFrame
    df_with_features = pd.concat([df_filtered.reset_index(drop=True), new_features], axis=1)

    # 3. Secondary Filtering: Remove rows where RDKit failed to calculate properties
    # This filters out the highly complex molecules that caused RDKit parse errors.
    final_df = df_with_features.dropna(subset=['MW']).copy() 
    
    filtered_count = len(final_df)

    print(f"Initial row count: {initial_count}")
    print(f"Filtered to small molecules (valid SMILES): {len(df_filtered)}")
    print(f"Final dataset size (RDKit calculable): {filtered_count}")
    print(f"Total rows removed: {initial_count - filtered_count}")
    
    # 4. Save the new, enriched DataFrame
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    print(f"\n✅ Data with features saved to {OUTPUT_FILE}")
    
    # Display summary of new data structure
    print("\nSample of the new Feature-Rich Dataset (Confirming features are numerical):")
    
    # Display the key columns (Identifiers, Targets, SMILES, and new Features)
    display_cols = ['DrugBank ID', 'UniProt_IDs', 'SMILES_String', 'MW', 'AlogP', 'HBD', 'HBA', 'TPSA']
    print(final_df[display_cols].head().to_markdown(index=False))

    print("\nData Structure after Calculation and Final Filtering:")
    print(final_df.info(verbose=False, memory_usage="deep"))

if __name__ == "__main__":
    run_feature_calculation()