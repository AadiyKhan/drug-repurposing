import pandas as pd
import numpy as np
import requests
import time
import pickle
import re
from pathlib import Path

# --- Configuration ---
INPUT_FILE = 'drugbank_cleaned_step1.csv'
OUTPUT_FILE = 'drugbank_normalized_targets.csv'
MAPPING_DICT_FILE = 'uniprot_target_mapping.pkl'

# NEW ENDPOINT
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# Rate limiting control (mandatory for repeated API calls)
SEARCH_DELAY_SECONDS = 0.5 # Delay between individual search queries

# --- API Interaction Function (Revised for Search Endpoint) ---

import re # Ensure this is at the top of your file!

# ... rest of the imports ...

# --- API Interaction Function (Revised for Search Endpoint) ---

def search_uniprot_for_id(target_name):
    """
    Searches UniProtKB for a target name and returns the top UniProt Accession ID.
    We clean the target name to prevent 400 Bad Request errors.
    """
    
    # 1. Clean the target name to remove characters that break the UniProt search query language
    # Removing parentheses as they are the immediate cause of the 400 error.
    cleaned_name = re.sub(r'[()]', '', target_name).strip()
    
    # If the name is empty after cleaning, skip it
    if not cleaned_name:
        return None

    # 2. Build the search query string
    # We use the cleaned name directly, without surrounding parentheses, 
    # and limit the query to human, reviewed proteins.
    query_str = f'{cleaned_name} AND (reviewed:true) AND (organism_id:9606)'
    
    # 3. Build the API request parameters
    params = {
        'query': query_str,
        'format': 'tsv',
        'fields': 'accession,id,organism_id', # We only need the accession ID
        'size': 1 # Only fetch the top match
    }
    
    try:
        response = requests.get(
            UNIPROT_SEARCH_URL,
            params=params,
            headers={"User-Agent": "DrugRepurposingProject"},
            timeout=30 # Set a timeout for the request
        )
        response.raise_for_status() # Raise exception for bad status codes

        # Read the TSV response (which includes a header row)
        data = response.text.strip().split('\n')
        
        if len(data) > 1:
            tsv_parts = data[1].split('\t')
            return tsv_parts[0] # Return the UniProt Accession (e.g., P00734)
        else:
            return None # No result found for this query

    except requests.exceptions.RequestException as e:
        # Log the failure, but return None to continue the process
        # We don't log the 400 error if we filtered the name correctly, but keep general logging
        if response.status_code != 400:
            print(f"  [ERROR] API request failed for '{target_name}': {e}")
        return None
# --- Main Processing Script ---

def run_target_normalization():
    """Loads, processes, and normalizes targets using the UniProt search API."""
    
    print(f"Loading data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE)
    
    # 1. Extract Unique Targets (Same as before)
    df_targets = df[df['Targets_Clean'] != 'NO_TARGETS_LISTED'].copy()
    all_targets = set()
    for targets_str in df_targets['Targets_Clean'].dropna():
        names = re.split(r'[;,\s|]', targets_str)
        for name in names:
            name = name.strip()
            if name and len(name) > 2:
                all_targets.add(name)

    unique_targets = sorted(list(all_targets))
    print(f"Extracted {len(unique_targets)} unique target names for mapping.")
    
    # 2. Sequential Search and Mapping
    target_map = {}
    successful_maps = 0
    total_targets = len(unique_targets)
    
    print(f"Starting sequential UniProt search for {total_targets} targets...")

    for i, name in enumerate(unique_targets):
        
        # Display progress update every 100 targets
        if (i + 1) % 100 == 0 or i == total_targets - 1:
            print(f"Processing target {i + 1}/{total_targets} ({name}). Found {successful_maps} IDs so far.")
            
        uniprot_id = search_uniprot_for_id(name)
        
        if uniprot_id:
            target_map[name] = uniprot_id
            successful_maps += 1
            
        # Enforce rate limit delay
        time.sleep(SEARCH_DELAY_SECONDS)

    print(f"\nTotal unique target names mapped to UniProt IDs: {successful_maps} / {total_targets}")
    
    # 3. Save the mapping dictionary (CRUCIAL artifact for Phase 3/4)
    with open(MAPPING_DICT_FILE, 'wb') as f:
        pickle.dump(target_map, f)
    print(f"✅ Target mapping dictionary saved to {MAPPING_DICT_FILE}")

    # 4. Apply Mapping to DataFrame (Same as before)
    
    def apply_uniprot_mapping(targets_str):
        """Converts raw target names in a cell to a pipe-separated list of UniProt IDs."""
        if targets_str == 'NO_TARGETS_LISTED':
            return targets_str
            
        names = re.split(r'[;,\s|]', targets_str)
        uniprot_ids = []
        for name in names:
            name = name.strip()
            if name in target_map:
                uniprot_ids.append(target_map[name])
                
        return " | ".join(sorted(list(set(uniprot_ids))))

    df['UniProt_IDs'] = df['Targets_Clean'].apply(apply_uniprot_mapping)
    df['UniProt_IDs'] = df['UniProt_IDs'].replace('', 'NO_IDS_FOUND')

    # 5. Final Save
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nFinal data saved to {OUTPUT_FILE}. Ready for EDA.")
    
    print("\nSample Data with new UniProt_IDs column:")
    print(df[['DrugBank ID', 'Name', 'Targets_Clean', 'UniProt_IDs']].head())

if __name__ == "__main__":
    run_target_normalization()