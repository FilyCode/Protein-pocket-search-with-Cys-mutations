# Protein Pocket Search with Cysteine Mutations

## Overview
This repository contains Python scripts designed for *in silico* protein engineering and property prediction, specifically focusing on identifying potential binding pockets through systematic Cysteine mutations. The pipeline utilizes the Schrödinger software suite via its API to perform structural modifications and predict various physicochemical characteristics of mutated proteins. This methodology supports the preliminary stages of drug discovery by providing data for subsequent ligand docking and molecular dynamics (MD) simulations.

## Features
*   **Systematic Cysteine Mutagenesis:** Iterates through each amino acid residue of a target protein, mutating it to Cysteine.
*   **In Silico Property Prediction:** For each mutant, the pipeline predicts key properties using Schrödinger's capabilities:
    *   **Solvent Accessible Surface Area (SASA):** Both absolute and relative SASA, including calculations around specific residues and within an 8Å radius.
    *   **Energetic Properties:** Potential energy, Prime energy, and interaction energies.
    *   **Protein Stability:** Changes in stability (e.g., d Stability).
    *   **pKa Values:** Prediction of pKa for mutated residues.
*   **Structural Minimization:** Automatically performs energy minimization on each mutated structure to ensure realistic conformations.
*   **Data Export:** Generates various data formats including `.mae` (Schrödinger structure files), `.txt`, `.h5` (HDF5), and `.csv` files for quantitative results.
*   **Automated Visualization:** Produces a range of plots to visualize predicted properties, including:
    *   Heatmaps of SASA and docking scores.
    *   Bar plots for individual residue properties.
    *   Combined plots integrating multiple property analyses.
*   **Binding Pocket Identification:** Incorporates a residue selection mechanism to filter potential binding sites based on predefined thresholds of SASA, energy, and stability changes.
*   **Docking Score Correction:** Includes functionality to correct ligand docking scores based on the heavy atom count of the ligands, improving comparability.

## Requirements
*   **Schrödinger Software Suite:** This pipeline heavily relies on the Schrödinger software, specifically its Python API. A licensed installation of Schrödinger is required.

## Usage
The scripts within this repository are designed to be executed in a Python environment where the Schrödinger API is accessible.

1.  **Configuration:** Update the file paths and other parameters within the scripts to match your specific project directory, input protein structure (`.mae` file), and desired output locations.
2.  **Input Structure:** Provide a prepared and minimized `.mae` file of your target protein (e.g., a KRAS-G12D structure as shown in the example code).
3.  **Execution:** Run the Python script(s) from your command line or an IDE. The scripts will:
    *   Load the reference protein structure.
    *   Calculate reference properties (e.g., SASA) for the wild-type/starting structure.
    *   Iterate through specified residues, mutate them to Cysteine, minimize the new structure, and calculate properties.
    *   Export all calculated data to various file formats.
    *   Generate and save plots visualizing the results.
    *   Perform residue selection based on defined criteria.
    *   Process and visualize docking scores (if relevant data is provided), including corrections.

**Note on Schrödinger API:** Ensure your environment is correctly configured to import `schrodinger` modules (e.g., by sourcing the Schrödinger `activate` script).

## Scope
These scripts focus solely on the *in silico* mutation process, prediction of protein properties, and initial identification of potential binding pockets. This repository **does not cover** the actual ligand docking simulations or subsequent molecular dynamics simulations, though it generates crucial input data and analyses for those downstream steps which were done in subsequent simulations.

## Authorship
This pipeline and its associated scripts were solely developed by Philipp Trollmann during his Bachelor's Thesis in Chemistry. The work was conducted at Boehringer Ingelheim within the group of Dr. Dirk Kessler and was academically supervised by Christian Schröder at the University of Vienna.
