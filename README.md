# Protein-Analysis-tool-Using-Biopython-package-and-Streamlit
Protein Structure & Sequence Analysis Pipeline

This project is a full end-to-end pipeline for analyzing protein structures from PDB files. It extracts sequences, aligns them, detects motifs, builds phylogenetic trees, and visualizes everything — including 3D protein structures — inside a notebook environment.

What This Pipeline Can Do

This tool allows you to:

Upload a PDB protein structure file

Extract amino-acid sequences from chains

Align the sequences using reference-based multiple sequence alignment

Identify conserved residues and generate motif logos

Build UPGMA and Neighbor-Joining phylogenetic trees

Perform bootstrap tree analysis

Highlight motifs directly on the 3D structure

Map hydrophobic residues in 3D

Optionally fetch additional records from GenBank and SwissProt

View a full visual summary of the analysis

Requirements

You’ll need Python and the following libraries installed:

biopython
nglview
matplotlib
ipywidgets
requests


To install dependencies:

pip install biopython nglview matplotlib ipywidgets requests

How To Run

Open the notebook or Python environment where you will run the script

Run the script

When prompted, upload a .pdb file

The pipeline will automatically:

extract sequences

align them

detect motifs

construct phylogenetic trees

display structural visualizations

Output Files

The pipeline will create a results/ folder containing:

motif_logo.png

upgma_tree.newick

nj_tree.newick

upgma_tree.png

nj_tree.png

You’ll also see all visualizations inside the notebook.

Code Structure (Simple Explanation)

The script consists of functions for:

Uploading the PDB file

Extracting peptide sequences

Multiple sequence alignment

Motif analysis

Phylogenetic tree construction

Bootstrap analysis

3D visualization

Hydrophobicity mapping

GenBank/SwissProt fetching

Summary display

At the end, the main_pipeline() function runs the entire workflow step by step.

Notes

This pipeline is intended for exploratory protein analysis and teaching purposes

Works best inside Jupyter Notebook / Google Colab due to interactive widgets

You can easily modify it to add more analyses

Example Use Case

Upload a PDB protein structure such as 1A3N.pdb, and the script will:

Extract its amino acid chains

Identify conserved motifs

Highlight them onto the actual 3D structure

Show evolutionary relationships among chains

Provide visual and data-driven results
