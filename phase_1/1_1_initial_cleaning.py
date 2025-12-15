import pandas as pd
import numpy as np
import re
from pathlib import Path

# --- Configuration ---
RAW_DATA_FILE = 'drugbank_data.csv'
CLEANED_DATA_FILE = 'drugbank_cleaned_step1.csv'
TARGET_COLUMNS = ['DrugBank ID', 'Name', 'Description', 'Targets']

# --- Cleaning Functions ---

def clean_text_for_nlp(text):
    """Applies basic text cleaning: lowercasing, removing non-essential characters."""
    if pd.isna(text):
        return text
    text = str(text).lower()
    # Remove citation brackets (e.g., [L41539] or [A246609]) which are common in DrugBank descriptions
    text = re.sub(r'\[[a-z0-9]+\]', '', text)
    # Remove common HTML/XML remnants if present (though less likely in a clean CSV)
    text = re.sub(r'<[^>]*>', '', text)
    # Normalize excessive whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    return text

def standardize_id_and_targets(item):
    """Ensures targets/IDs are treated as strings or NaN."""
    if pd.isna(item):
        return np.nan
    return str(item).strip()


## --- Main Processing Script ---

def run_initial_cleaning():
    """Loads, cleans, and saves the initial DrugBank data."""
    
    if not Path(RAW_DATA_FILE).exists():
        print(f"Error: Raw data file '{RAW_DATA_FILE}' not found.")
        print("Please ensure 'drugbank_data.csv' is in the current directory.")
        return

    print(f"Loading raw data from {RAW_DATA_FILE}...")
    df = pd.read_csv(RAW_DATA_FILE, usecols=TARGET_COLUMNS)

    print(f"Initial shape: {df.shape}")

    # 1. Standardize ID and Name columns
    for col in ['DrugBank ID', 'Name']:
        df[col] = df[col].apply(standardize_id_and_targets)

    # 2. Clean the Description column (NLP Feature)
    df['Description_Clean'] = df['Description'].apply(clean_text_for_nlp)

    # 3. Handle Missing Values in Core Columns
    
    # A. Drop records with missing DrugBank ID or Name (essential identifiers)
    df.dropna(subset=['DrugBank ID', 'Name'], inplace=True)
    print(f"Shape after dropping null IDs/Names: {df.shape}")

    # B. Impute missing Targets and Descriptions (important for later feature generation)
    
    # Missing Targets are crucial: treat them as 'NO_TARGETS_LISTED'
    df['Targets_Clean'] = df['Targets'].fillna('NO_TARGETS_LISTED').apply(standardize_id_and_targets)
    
    # Missing Descriptions: treat them as 'NO_DESCRIPTION' (clean text version)
    df['Description_Clean'] = df['Description_Clean'].fillna('no description')

    # 4. Final Data Check: Remove duplicate DrugBank IDs (keeping the first entry)
    df.drop_duplicates(subset=['DrugBank ID'], keep='first', inplace=True)
    print(f"Final shape after deduplication: {df.shape}")

    # Select final columns for the next phase
    final_cols = ['DrugBank ID', 'Name', 'Description_Clean', 'Targets_Clean']
    df_cleaned = df[final_cols]
    
    # Save the cleaned intermediate file
    df_cleaned.to_csv(CLEANED_DATA_FILE, index=False)
    print(f"\n✅ Initial cleaning complete. Data saved to {CLEANED_DATA_FILE}")

    # Display a sample of the cleaned data
    print("\nSample Cleaned Data:")
    print(df_cleaned.head(3))

if __name__ == "__main__":
    # Create the environment for running the script
    if not Path('./.venv').exists():
        print("Please create and activate your virtual environment first (python -m venv .venv).")
    else:
        run_initial_cleaning()