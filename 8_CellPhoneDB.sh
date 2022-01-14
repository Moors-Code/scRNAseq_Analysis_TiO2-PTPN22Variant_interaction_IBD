#!/bin/bash
#Activate cellphoneDB conda environment 
source /mnt/khandler/cell_phone_db/bin/activate

#define output directory (WD) and directory where the input files for CellPhoneDB are stored (CELLFPHONE_DIR)
WD="/mnt/khandler/RStudio_Data/IBD_zurich_marlene/code_for_github/OUTPUTS"
CELLPHONE_DIR="/mnt/khandler/RStudio_Data/IBD_zurich_marlene/code_for_github"

#Create directories for outputs 
mkdir ${WD}/ouputs_wtDSS
mkdir ${WD}/ouputs_wtTiO2
mkdir ${WD}/ouputs_mtDSS
mkdir ${WD}/ouputs_mtTiO2

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_wtDSS.txt ${CELLPHONE_DIR}/cellphonedb_count_wtDSS.txt --output-path=${WD}/ouputs_wtDSS --project-name="IBD_wtDSS" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_wtTiO2.txt ${CELLPHONE_DIR}/cellphonedb_count_wtTiO2.txt --output-path=${WD}/ouputs_wtTiO2 --project-name="IBD_wtTiO2" --counts-data=ensembl 
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_mtDSS.txt ${CELLPHONE_DIR}/cellphonedb_count_mtDSS.txt --output-path=${WD}/ouputs_mtDSS --project-name="IBD_mtDSS" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_mtTiO2.txt ${CELLPHONE_DIR}/cellphonedb_count_mtTiO2.txt --output-path=${WD}/ouputs_mtTiO2 --project-name="IBD_mtTiO2" --counts-data=ensembl 
