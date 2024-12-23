# Correlations-with-physico-chemical-properties-of-a-curated-database
An RDKit project to apply Python in data analysis, correlations and predictions for cheminformatics 
README: Molecular Analysis and Correlation Study
This project performs molecular analysis on a dataset of compounds using RDKit, focusing on calculating molecular descriptors, such as the Synthetic Accessibility Score (SAScore) and logP. The analysis also includes calculating the Tanimoto similarity between molecules, exploring correlations between key molecular descriptors, and visualizing these relationships using heatmaps and scatter plots.

Before starting the analysis, you need to install the necessary libraries. You can do this by running the following command:

pip install rdkit pandas matplotlib seaborn
RDKit: Used to handle and process molecular data.
pandas: For data manipulation and analysis.
matplotlib & seaborn: For plotting and visualization.
Loading Data
You can load your molecular data from a CSV file containing SMILES representations of compounds:

python
import pandas as pd

# Load the CSV file containing SMILES
df = pd.read_csv('Syntons_5567.csv')

# Ensure the SMILES column exists
if 'SMILES' not in df.columns:
    raise ValueError("CSV file must have a column named 'SMILES'.")
The data is loaded into a pandas DataFrame, and the SMILES strings are converted into RDKit molecule objects for further analysis.

Molecular Descriptors
In this project, two important molecular descriptors are calculated for each compound:

logP (Partition Coefficient): A measure of the compound’s lipophilicity.
Synthetic Accessibility Score (SAScore): A rough estimate of how synthetically accessible a compound is. This is calculated using a custom formula based on molecular weight, number of rotatable bonds, and number of rings:
python
Copy code
from rdkit.Chem import Descriptors

# Function to calculate the SAScore
def calculate_sascore(mol):
    if mol is None:
        return None
    mol_weight = Descriptors.MolWt(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    num_rings = Descriptors.RingCount(mol)
    sascore = (mol_weight + num_rotatable_bonds + num_rings) / 3.0
    return sascore
Calculating Tanimoto Similarity
The Tanimoto similarity is used to quantify how similar two molecular structures are based on their fingerprints. In this case, it is calculated using RDKit's GetMorganFingerprint method, and the similarity is calculated as follows:

from rdkit.Chem import DataStructs
from rdkit.Chem import RDKFingerprint

# Function to calculate Tanimoto similarity between two molecules
def calculate_tanimoto_similarity(mol1, mol2):
    fp1 = RDKFingerprint(mol1)
    fp2 = RDKFingerprint(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)
The resulting Tanimoto similarity score ranges from 0 to 1, where 1 indicates identical structures.

Correlation Analysis
To explore the relationships between molecular descriptors, the correlation coefficient between logP and SAScore is calculated:

# Calculate the correlation coefficient
correlation = df['logP'].corr(df['Synthetic Accessibility Score'])
print(f"Correlation between logP and Synthetic Accessibility Score: {correlation:.2f}")
Visualizations
Several visualizations are created to help analyze the data:

Scatter plot: A scatter plot of logP vs. Synthetic Accessibility Score helps visualize the relationship between the two descriptors.
python
Copy code
import matplotlib.pyplot as plt
import seaborn as sns

# Create a scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x='logP', y='Synthetic Accessibility Score', data=df)

# Add a regression line (optional)
sns.regplot(x='logP', y='Synthetic Accessibility Score', data=df, scatter=False, color='r', line_kws={'linewidth': 2})

# Set labels and title
plt.title("Correlation between logP and Synthetic Accessibility Score", fontsize=14)
plt.xlabel("logP", fontsize=12)
plt.ylabel("Synthetic Accessibility Score", fontsize=12)

# Display the plot
plt.show()
Heatmap: A heatmap of correlations between all numerical descriptors in the dataset is created to observe broader trends in the data.

# Create a heatmap of correlations
correlation_matrix = df.corr()
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", cbar=True)
plt.title("Correlation Heatmap of Molecular Descriptors", fontsize=14)
plt.show()
Conclusion
In this project, we:

Loaded a dataset of compounds and converted SMILES strings to RDKit molecule objects.
Calculated molecular descriptors such as logP and Synthetic Accessibility Score.
Analyzed the correlation between logP and SAScore.
Calculated Tanimoto similarity between pairs of molecules.
Created visualizations to explore these relationships using scatter plots and heatmaps.
This analysis can be extended further by including more molecular descriptors and performing additional statistical analyses.

Notes
Ensure that the CSV file you're working with includes the necessary SMILES column.
The custom SAScore formula is a rough approximation. More complex models for synthetic accessibility can be used if needed.
