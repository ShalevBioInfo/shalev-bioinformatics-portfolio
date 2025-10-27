# main_analysis.py
# Demo script for Cilia Clusters project (synthetic data example)
# Author: Shalev Yaacov
# Date: 2025-10-27

import numpy as np
import pandas as pd

# --- Step 1: Create synthetic demo matrix ---
genes = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
species = ["Sp1", "Sp2", "Sp3", "Sp4", "Sp5"]
np.random.seed(42)
demo_matrix = pd.DataFrame(np.random.randn(5, 5), index=genes, columns=species)

# --- Step 2: Save demo matrix to CSV ---
demo_matrix.to_csv("results/demo_matrix.csv")

# --- Step 3: Print summary ---
print("Demo matrix created and saved to results/demo_matrix.csv")
