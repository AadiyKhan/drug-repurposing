#pip install gensim nltk
import pandas as pd
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from nltk.tokenize import word_tokenize
from pathlib import Path
import re
import numpy as np

# Note: You may need to run 'pip install gensim nltk' and run the following line once:
# import nltk; nltk.download('punkt')

# --- Configuration ---
# Use the file that contains the Description_Clean column and the new features
INPUT_FILE = 'drugbank_morgan_fingerprints.csv' 
OUTPUT_FILE = 'drugbank_final_features_with_nlp.csv'

# Doc2Vec Model Parameters
VECTOR_SIZE = 100 # Size of the final description vector (embedding)
WINDOW_SIZE = 5   
MIN_COUNT = 2     
EPOCHS = 20       

# --- Preprocessing and Tokenization ---

def preprocess_text(text):
    """Clean and tokenize the text for Doc2Vec."""
    if pd.isna(text) or text == 'no description':
        # Handle missing or empty descriptions gracefully
        return ["no_description_placeholder"]
    
    # Remove special characters, multiple spaces, and convert to lowercase
    text = re.sub(r'[^\w\s]', '', text.lower())
    # Tokenize the text (split into words)
    return word_tokenize(text)

# --- Main Processing Script ---

def run_nlp_feature_generation():
    
    if not Path(INPUT_FILE).exists():
        print(f"Error: Required input file '{INPUT_FILE}' not found.")
        print("Please ensure 'drugbank_morgan_fingerprints.csv' is in the directory.")
        return

    print(f"Loading data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE)
    
    print("Preprocessing and tokenizing drug descriptions...")
    
    # Apply preprocessing to the description column
    tokenized_data = df['Description_Clean'].apply(preprocess_text)
    
    # Tagging: Preparing data for Gensim (Each document needs a tag/label)
    # The tag is simply the row index 'i'
    tagged_data = [TaggedDocument(doc, [i]) for i, doc in enumerate(tokenized_data)]

    # 1. Initialize and Train the Doc2Vec Model
    print(f"Training Doc2Vec model (Vector Size: {VECTOR_SIZE}, Epochs: {EPOCHS})...")
    
    model = Doc2Vec(vector_size=VECTOR_SIZE, 
                    window=WINDOW_SIZE, 
                    min_count=MIN_COUNT, 
                    workers=4, 
                    epochs=EPOCHS)
    
    model.build_vocab(tagged_data)
    model.train(tagged_data, total_examples=model.corpus_count, epochs=model.epochs)
    
    # 2. Infer Vectors (Generate the embedding for each document)
    print("Inferring document vectors (NLP embeddings)...")
    
    # Create a numerical array to store the vectors
    nlp_vectors = np.zeros((len(df), VECTOR_SIZE))
    
    for i, tokens in enumerate(tokenized_data):
        # The infer_vector method generates the final embedding for the tokens
        vector = model.infer_vector(tokens, epochs=model.epochs)
        nlp_vectors[i] = vector

    # 3. Add the NLP features to the DataFrame
    # Convert the array into a DataFrame
    df_nlp_features = pd.DataFrame(
        nlp_vectors, 
        columns=[f'NLP_Vec_{i}' for i in range(VECTOR_SIZE)]
    )
    
    # Concatenate the new NLP features with the existing DataFrame
    # Since we used df.index for the tokens, we can concatenate directly
    df_final = pd.concat([df, df_nlp_features], axis=1)

    # 4. Save the new, feature-rich DataFrame
    df_final.to_csv(OUTPUT_FILE, index=False)
    
    print(f"\n✅ Final feature dataset with {VECTOR_SIZE} NLP embeddings saved to {OUTPUT_FILE}")
    
    # Display summary of new data structure
    print("\nSample of the Final Feature-Rich Dataset (showing new NLP columns):")
    
    # Display the new columns to confirm they are numerical
    display_cols = ['DrugBank ID', 'Name', 'MW', 'AlogP'] + [f'NLP_Vec_{i}' for i in range(5)]
    print(df_final[display_cols].head().to_markdown(index=False, numalign="left", stralign="left"))

    print("\nData Structure of the Final Dataset:")
    print(df_final.info(verbose=False, memory_usage="deep"))

if __name__ == "__main__":
    run_nlp_feature_generation()