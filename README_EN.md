# Antibody Sequence Extraction and Clustering

## Description
This script is designed to extract amino acid sequences of antibodies from PDB files using the SAbDab database. 

## Installation and Environment Setup
1. Install Anaconda (if not already installed) from the official website.
2. Clone the repository using the command `git clone https://github.com/SergeiNikolenko/ml.git`.
3. Create a new environment in Anaconda using the command `conda env create -f environment.yml`.
4. Activate the environment using the command `conda activate ml4`.

## Usage
1. Download the SAbDab database from the official website (https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/archive/all/) and extract it.
2. Place all PDB files in the `all_structures/chothia` folder.
3. Launch `ml.ipynb`, which contains all the scripts for processing and instructions.

## Functionality
1. Extraction of amino acid sequences for heavy and light chains of antibodies from PDB files. The sequences are saved in separate FASTA files in the `heavy_chains` and `light_chains` folders.
2. Removal of empty files. The script checks the FASTA files in the `heavy_chains` and `light_chains` folders for empty files and deletes them.
3. Clustering of amino acid sequences using DBSCAN, K-mean, and hierarchical clustering algorithms.
4. Analysis of clustering.

## Interpretation of Results
The clustering results can be used to analyze the structure and properties of antibodies in different samples. You can also use the code from the script for your own research.

## Usage Instructions
To use the script, follow these steps:

1. Download the SAbDab database from https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/archive/all/ and extract it. The script uses the `all_structures/chothia` folder, so make sure it is present after extraction.
2. Launch the `ml.ipynb` script in Jupyter Notebook or Jupyter Lab.
3. Run the cell with the desired function.
4. Analyze the results and create visualizations using the generated data.

## Additional Notes
This script can process a large number of PDB files with some delay. If you encounter problems with processing large files, it is recommended to split the files into smaller parts and run the script on each part separately.

It is also worth noting that this script was written as part of an educational project and can be improved and expanded. If you have any suggestions or feedback, please contact the author.
