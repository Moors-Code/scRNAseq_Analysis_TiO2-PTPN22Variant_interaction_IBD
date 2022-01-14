# TiO2 nanoparticles abrogate the protective effect ot the Crohn´s disease-associated variation within the PTPN22 gene locus

This repository contains code from single-cell RNAseq data analysis. 

# Citation 

The code in this repository pertains to publication: 

TiO2 nanoparticles abrogate the protective effect ot the Crohn´s disease-associated variation within the PTPN22 gene locus

Marlene Schwarzfischer, Anna Niechcial, Kristina Handler, Yasser Morsy, Marcin Wawrzyniak, Andrea Laimbacher, Kirstin Atrott, Roberto Manzini, Katharina Baebler, Larissa Hering, Egle Katkeviciutė, Janine Häfliger, Silvia Lang, Maja Keller, Jérome Woodtli, Lisa Eisenbeiss, Thomas Krämer, Elisabeth M Schraner, Mahesa Wiesendanger, Sebastian Zeissig, Gerhard Rogler, Andreas Moor, Michael Scharl and Marianne R Spalinger


# Data
Data for the above publication can be found at SRA. Access number: 


# Code description: 

R code: 

1_Packages_and_functions.R: This code lists all the packages used for R analysis and functions used in other R scripts. 

2_Preprocessing_and_clustering.R: Code for pre-processing and clustering of 12 samples together. Takes output gene count matrices from BD Seven Bridges platform from 12 samples as input and converts them to a merged Seurat object with proper quality control cut-offs. 

3_Annotation_and_clustering_comparison.R: Code for annotation of cell clusters and comparison of cell cluster abundance between different conditions. 

4_Subset_CD8_T_cells.R : Code for closer analysis of CD8 T cell cluster. Cells from this cluster are re-clustered and subtypes are annotated. Then differences in cell cluster abundance, gene expression and gene set enrichments between different conditions are analyzed.  

5_Subset_CD4_T_cells.R : Code for closer analysis of CD4 T cell cluster. Cells from this cluster are re-clustered and subtypes are annotated. Then differences in cell cluster abundance and gene expression between different conditions are analyzed.  

6_Subset_Macrophages.R : Code for closer analysis of macrophages cell cluster. Cells from this cluster are re-clustered and subtypes are annotated. Then differences in cell cluster abundance and gene expression between different conditions are analyzed.  

7_CellPHoneDB_input_files_generation.R: This code produces input files for the shell script 8_CellPhoneDB.sh; a gene count matrix from single cells and a meta data file containing the annotation of single cells are produced from each of the 4 different conditions. It also converts mouse gene ensembl symbols to human gene ensembl symbols because CellPhoneDB is based on human ligand-receptor databases. 

Shell script: 

8_CellPHoneDB.sh: This code runs the python package CellPhoneDB within a conda environment for all 4 different conditions separately. 


