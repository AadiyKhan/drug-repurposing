import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Load dataset
df = pd.read_csv("drugbank_small_molecules_features.csv")

# CHANGE this if your column name is different
SMILES_COLUMN = "SMILES_String"

# Morgan fingerprint parameters
RADIUS = 2          # ECFP4
N_BITS = 2048       # fingerprint length

def smiles_to_morgan(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=RADIUS, nBits=N_BITS
    )
    return np.array(fp)

# Generate fingerprints
fingerprints = []
valid_rows = []

for idx, smiles in enumerate(df[SMILES_COLUMN]):
    fp = smiles_to_morgan(smiles)
    if fp is not None:
        fingerprints.append(fp)
        valid_rows.append(idx)

# Keep only valid molecules
df = df.iloc[valid_rows].reset_index(drop=True)

# Convert fingerprints to DataFrame
fp_df = pd.DataFrame(
    fingerprints,
    columns=[f"FP_{i}" for i in range(N_BITS)]
)

# Combine original data + fingerprints
final_df = pd.concat([df, fp_df], axis=1)

# Save output
final_df.to_csv("drugbank_morgan_fingerprints.csv", index=False)

print("✅ Morgan fingerprints generated!")
print("📁 Output saved as: drugbank_morgan_fingerprints.csv")
print("🧬 Fingerprint shape:", fp_df.shape)
