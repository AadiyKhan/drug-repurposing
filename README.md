![Language](https://img.shields.io/badge/Language-Python-blue)
![License](https://img.shields.io/badge/License-MIT-green)

# Drug Repurposing Project

This project leverages computational methods to identify potential new uses for existing drugs, accelerating the discovery of treatments for various diseases. By integrating chemical properties, biological targets, and disease associations, it builds a comprehensive knowledge graph and employs machine learning to predict novel drug-disease indications.

## 📝 Table of Contents

- [✨ Features](#-features)
- [🛠️ Tech Stack](#️-tech-stack)
- [🚀 Installation](#-installation)
- [📊 Usage](#-usage)
- [📂 Project Structure](#-project-structure)
- [⚙️ Configuration](#️-configuration)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [👤 Author](#-author)

## ✨ Features

*   **Automated Data Ingestion & Cleaning**: Scripts to fetch, parse, and clean drug and target data from public databases.
*   **Comprehensive Feature Engineering**: Calculates molecular properties, generates chemical fingerprints (e.g., Morgan), and extracts insights from drug descriptions using Natural Language Processing (NLP).
*   **Knowledge Graph Construction**: Assembles a unified knowledge graph connecting drugs, targets, and diseases, enabling sophisticated relationship analysis.
*   **Predictive Candidate Identification**: Utilizes graph-based methods and computational models to predict novel drug-disease repurposing candidates.
*   **Interactive Visualization Dashboard**: Provides a Streamlit-powered web application for exploring discovered repurposing candidates and their associated data.

## 🛠️ Tech Stack

*   **Python**: The core language for all data processing, analysis, and model building.
    *   **NumPy**: Essential for high-performance numerical operations, especially with large datasets.
    *   **Pandas**: Used extensively for data manipulation and analysis of structured data.
    *   **Streamlit**: Powers the interactive web dashboard for visualizing results.
    *   **PubChemPy**: Used for programmatic access and retrieval of chemical information from PubChem.
    *   **Altair / PyDeck**: Libraries likely used for generating interactive and visually rich plots and maps within the dashboard.
*   **Express.js / JavaScript**: While not central to the core logic, their presence (e.g., in `pydeck` and `altair` entry points) suggests they are used for interactive web components or underlying visualization frameworks.

## 🚀 Installation

Follow these steps to set up the project locally.

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/your-username/drug-repurposing.git
    cd drug-repurposing

    ```

2.  **Set up a Python virtual environment:**

    It is highly recommended to use a virtual environment to manage project dependencies.

    ```bash
    python -m venv .venv

    ```

3.  **Activate the virtual environment:**

    *   **On Windows:**
        ```bash
        .venv\Scripts\activate

        ```
    *   **On macOS/Linux:**
        ```bash
        source .venv/bin/activate

        ```

4.  **Install dependencies:**

    This project uses several Python libraries. Since no `requirements.txt` was detected, you might need to create one or install key packages manually.

    To generate a `requirements.txt` from your active environment (if you have installed packages globally):
    ```bash
    pip freeze > requirements.txt

    ```

    Alternatively, manually install common packages inferred from the project:

    ```bash
    pip install pandas numpy scikit-learn matplotlib seaborn streamlit rdkit networkx pubchempy altair pydeck

    ```

    If `requirements.txt` exists, install dependencies as follows:

    ```bash
    pip install -r requirements.txt

    ```

## 📊 Usage

The project is structured into several phases, designed to be run sequentially, culminating in a dashboard for interactive exploration.

### Running the Data Processing Pipeline

Navigate through the `phase_X` directories and execute the Python scripts in order.

1.  **Phase 1: Initial Data Preparation**
    ```bash
    python phase_1/1_1_initial_cleaning.py
    python phase_1/1_2_target_normalization.py
    python phase_1/1_3_fetch_smiles.py

    ```

2.  **Phase 2: Feature Engineering**
    ```bash
    python phase_2/2_1_calculate_properties.py
    python phase_2/2_2_generate_fingerprints.py
    python phase_2/2_3_nlp_feature_engineering.py

    ```

3.  **Phase 3: Knowledge Graph & Candidate Prediction**
    ```bash
    python phase_3/3_1_target_enrichment_analysis.py
    python phase_3/3_2_build_kg_nodes.py
    python phase_3/3_3_build_kg_edges.py
    python phase_3/3_4_predict_candidates.py

    ```

4.  **Phase 4: Visualization & Recommendations**
    ```bash
    python phase_4/4_1_visualize_discoveries.py
    python phase_4/4_2_final_recommendation.py

    ```

### Running the Discovery Dashboard

After processing all phases, you can launch the interactive Streamlit dashboard:

```bash
streamlit run phase_5/discovery_dashboard.py

```

This command will open the dashboard in your default web browser, typically at `http://localhost:8501`.

## 📂 Project Structure

```
├── README.md
├── datset/
│   ├── drugbank_cleaned_step1.csv
│   ├── drugbank_data.csv
│   ├── drugbank_morgan_fingerprints.csv
│   ├── drugbank_normalized_targets.csv
│   ├── drugbank_small_molecules_features.csv
│   ├── drugbank_vocabulary.csv
│   └── drugbank_with_smiles.csv
├── disease_target_enrichment.csv
├── disgenet_gda_subset.csv
├── drug_repurposing_graph.pkl
├── drug_repurposing_graph_final.pkl
├── drugbank_final_features_with_nlp.csv
├── final_recommendation_shortlist.csv
├── phase_1/
│   ├── 1_1_initial_cleaning.py
│   ├── 1_2_target_normalization.py
│   └── 1_3_fetch_smiles.py
├── phase_2/
│   ├── 2_1_calculate_properties.py
│   ├── 2_2_generate_fingerprints.py
│   └── 2_3_nlp_feature_engineering.py
├── phase_3/
│   ├── 3_1_target_enrichment_analysis.py
│   ├── 3_2_build_kg_nodes.py
│   ├── 3_3_build_kg_edges.py
│   └── 3_4_predict_candidates.py
├── phase_4/
│   ├── 4_1_visualize_discoveries.py
│   └── 4_2_final_recommendation.py
├── phase_5/
│   └── discovery_dashboard.py
├── repurposing_candidates_final.csv
├── repurposing_discovery_map.png
└── uniprot_target_mapping.pkl

```

*   **`dataset/`**: Contains raw and intermediate datasets used throughout the project.
*   **`phase_1/`**: Scripts for initial data cleaning, normalization, and fetching essential drug information.
*   **`phase_2/`**: Scripts for advanced feature engineering, including molecular property calculation, fingerprint generation, and NLP features.
*   **`phase_3/`**: Focuses on building the knowledge graph, performing target enrichment analysis, and predicting drug candidates.
*   **`phase_4/`**: Scripts for visualizing the discoveries and generating final recommendation lists.
*   **`phase_5/`**: Contains the Streamlit dashboard application for interactive data exploration.
*   **`.csv`, `.pkl`, `.png` files**: Various output files from different processing stages, including the final recommendations and visualizations.

## ⚙️ Configuration

No specific configuration files (e.g., `.env`, `config.ini`) were detected. All configurations, such as API keys or file paths, are currently embedded within the Python scripts. For customization, you will need to modify the relevant script files directly.

## 🤝 Contributing

Contributions are welcome! If you have suggestions for improving the project, please open an issue or submit a pull request.

## 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## 👤 Author

[Aadiy Khan / student @VIT Bhopal University]
[Ghaziah Shoeb / student @VIT Bhopal University]
