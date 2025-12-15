import pandas as pd
import requests
import time
from pathlib import Path
from urllib.parse import quote

# --- Configuration ---
INPUT_FILE = 'drugbank_normalized_targets.csv'
OUTPUT_FILE = 'drugbank_with_smiles.csv'
DRUGBANK_VOCABULARY_FILE = 'drugbank_vocabulary.csv' # Required local file
CACTUS_DELAY_SECONDS = 0.5 # Delay to respect the public server's load
REQUEST_TIMEOUT = 30 
CACTUS_BASE_URL = "https://cactus.nci.nih.gov/chemical/structure/"

# --- Local Data Loading ---

def load_inchikey_map(filename):
    """
    Loads DrugBank ID -> Standard InChI Key mapping from the local vocabulary file.
    """
    print(f"Loading official DrugBank vocabulary file: {filename} to get InChI Keys...")
    file_path = Path(filename)
    
    if not file_path.exists():
        print(f"🚨 ERROR: Required file '{filename}' not found.")
        print("Please ensure 'drugbank_vocabulary.csv' is in your project directory.")
        return None

    try:
        # Load the vocabulary file, using the two identified columns
        df_map = pd.read_csv(
            file_path, 
            usecols=['DrugBank ID', 'Standard InChI Key'], 
            dtype={'DrugBank ID': str, 'Standard InChI Key': str}
        )
        
        # Drop rows where InChI Key is missing (biologics/non-chemicals)
        df_map.dropna(subset=['Standard InChI Key'], inplace=True)
        
        # Create a dictionary map for fast lookup
        inchikey_map = df_map.set_index('DrugBank ID')['Standard InChI Key'].to_dict()
        
        print(f"  Loaded {len(inchikey_map)} DrugBank ID to InChI Key mappings.")
        return inchikey_map
        
    except KeyError:
        print("🚨 ERROR: The downloaded file is missing required columns ('DrugBank ID' or 'Standard InChI Key').")
        return None
    except Exception as e:
        print(f"🚨 ERROR loading {filename}: {e}")
        return None

# --- Cactus API Function (Highly Stable InChIKey Lookup) ---

def fetch_smiles_from_inchikey_cactus(inchikey):
    """
    Fetches Canonical SMILES from the Cactus NCI Resolver using InChIKey.
    """
    encoded_inchikey = quote(inchikey)
    
    # API structure: BASE_URL / [identifier] / [representation]
    smiles_url = f"{CACTUS_BASE_URL}{encoded_inchikey}/smiles"
    
    try:
        # Request the SMILES string
        response_smiles = requests.get(smiles_url, timeout=REQUEST_TIMEOUT)
        
        # 404 means the identifier was not found (usually because it's a biologic or a very new/rare compound)
        if response_smiles.status_code == 404:
            return None
            
        # Raise an error for other HTTP failures (e.g., 500 server error)
        response_smiles.raise_for_status()
        
        # The result is the raw SMILES string (plus potentially a newline/whitespace)
        smiles = response_smiles.text.strip()
        
        # The API sometimes returns a message like "None" if lookup fails but returns 200
        if smiles and "cannot be resolved" not in smiles:
            return smiles
        
        return None
        
    except requests.exceptions.RequestException:
        # Catch network timeout or other non-HTTP errors
        return None
    except Exception:
        return None


# --- Main Processing Script ---

def run_smiles_fetching():
    
    # 0. Load the local DrugBank ID -> InChI Key map
    inchikey_map = load_inchikey_map(DRUGBANK_VOCABULARY_FILE)
    if inchikey_map is None:
        return

    if not Path(INPUT_FILE).exists():
        print(f"\nError: Required input file '{INPUT_FILE}' not found.")
        return

    print(f"\nLoading data from {INPUT_FILE} for SMILES retrieval...")
    df = pd.read_csv(INPUT_FILE)
    
    if 'SMILES_String' in df.columns and df['SMILES_String'].first_valid_index() is not None:
        print("SMILES_String column already exists and contains data. Skipping retrieval.")
        df.to_csv(OUTPUT_FILE, index=False)
        return

    df['SMILES_String'] = None
    total_drugs = len(df)
    success_count = 0
    
    print(f"Starting hybrid lookup for {total_drugs} drugs (Local InChIKey -> CACTUS SMILES)...")

    for index, row in df.iterrows():
        db_id = row['DrugBank ID']
        
        # 1. Local Lookup: Get InChIKey
        inchikey = inchikey_map.get(db_id)

        if inchikey:
            # 2. API Lookup: Get SMILES from CACTUS using InChIKey
            smiles = fetch_smiles_from_inchikey_cactus(inchikey)
            
            if smiles:
                df.at[index, 'SMILES_String'] = smiles
                success_count += 1
            
            # Pause to respect the public API's rate limits
            time.sleep(CACTUS_DELAY_SECONDS)
        
        if (index + 1) % 100 == 0:
            print(f"Processing drug {index + 1}/{total_drugs}. Found {success_count} SMILES strings so far.")
        # Only print updates for the first ~500 drugs
        elif (index + 1) % 10 == 0 and index < 500:
            print(f"  -- Current Drug: {index + 1}/{total_drugs} ({db_id}) --")
            
    print(f"\nTotal SMILES strings successfully fetched: {success_count} / {total_drugs}")
    
    df.to_csv(OUTPUT_FILE, index=False)
    
    print(f"✅ Data with SMILES saved to {OUTPUT_FILE}")
    
    print("\nSample of new 'SMILES_String' column (First 5 drugs that have SMILES):")
    
    # Filter for the first 5 drugs that actually have SMILES for a good sample
    smiles_sample = df[df['SMILES_String'].notna()].head(5)
    
    if smiles_sample.empty:
         # Show the first 5 rows even if they are biologics
         print(df[['DrugBank ID', 'SMILES_String']].head(5).to_markdown(index=False))
         print("\nNote: The first drugs in the file are biologics. If this is a final run, the SMILES count should be high.")
    else:
        print(smiles_sample[['DrugBank ID', 'SMILES_String']].to_markdown(index=False))

    
if __name__ == "__main__":
    run_smiles_fetching()