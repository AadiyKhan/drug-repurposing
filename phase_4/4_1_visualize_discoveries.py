import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FILE = 'repurposing_candidates_final.csv'
OUTPUT_IMAGE = 'repurposing_discovery_map.png'

def create_discovery_heatmap():
    print("Generating Discovery Map...")
    df = pd.read_csv(INPUT_FILE)

    if df.empty:
        print("No candidates to visualize.")
        return

    # Pivot the data: Rows = Drugs, Columns = Diseases, Values = Scores
    # We'll take the top 20 drugs overall to keep the map readable
    top_drugs = df.groupby('Candidate_Name')['Repurposing_Score'].max().sort_values(ascending=False).head(20).index
    df_filtered = df[df['Candidate_Name'].isin(top_drugs)]

    heatmap_data = df_filtered.pivot_table(
        index='Candidate_Name', 
        columns='Disease', 
        values='Repurposing_Score'
    ).fillna(0)

    # Plotting
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_data, annot=True, cmap="YlGnBu", fmt=".1f", cbar_kws={'label': 'Repurposing Score'})
    
    plt.title('Top 20 Drug Repurposing Candidates across Diseases', fontsize=16, pad=20)
    plt.xlabel('Disease Category', fontsize=12)
    plt.ylabel('Candidate Drug Name', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    plt.savefig(OUTPUT_IMAGE, dpi=300)
    print(f"✅ Discovery Map saved to: {OUTPUT_IMAGE}")

if __name__ == "__main__":
    create_discovery_heatmap()