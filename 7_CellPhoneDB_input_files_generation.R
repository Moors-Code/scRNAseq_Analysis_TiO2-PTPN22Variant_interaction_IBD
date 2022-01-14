##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

##### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### load R object after preprocessing, lustering and annotation ##### 
merged_ann <- readRDS(file = "merged_ann.Rds")

##### Include annotation of CD8 T cell subclusters to analyse L-R interactions of cytotoxic T cells with other cell types #####
#subset all cell types and then merge while using the re-annotated CD8 T cell object (cd8 )
Idents(merged_ann) <- "annotation"
Macs <- subset(merged_ann, idents = "Macrophages")
Endo <- subset(merged_ann, idents = "Endothelial")
DC1 <- subset(merged_ann, idents = "CD103+ DC")
DC2 <- subset(merged_ann, idents = "CD11b+ DC")
ILC2 <- subset(merged_ann, idents = "ILC2")
Fib <- subset(merged_ann, idents = "Fibroblasts")
Pl <- subset(merged_ann, idents = "Plasma cells")
CD4 <- subset(merged_ann, idents = "CD4+ T")
Eos <- subset(merged_ann, idents = "Eosinophils")
B <- subset(merged_ann, idents = "B")
ILC3 <- subset(merged_ann, idents = "ILC3")
Neut <- subset(merged_ann, idents = "Neutrophils")

#load cd8 seurat object 
CD8 <- readRDS(file = "CD8_seurat_object_reclustered.Rds")

#merge all objects 
merged_ann2 <- merge(CD8, c(Macs,Endo,DC1,DC2,ILC2,Fib,Pl,CD4,Eos,B, ILC3,Neut))
#control presence of all cell types 
table(merged_ann2$annotation, merged_ann2$condition)

##### Generation of input files for each condition ##### 
#load human and mouse ensemble symbols 
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

Input_files_CellPhoneDB_generation(merged_ann2, "wtDSS", "cellphonedb_count_wtDSS.txt","cellphonedb_meta_wtDSS.txt")
Input_files_CellPhoneDB_generation(merged_ann2, "mtDSS", "cellphonedb_count_mtDSS.txt","cellphonedb_meta_mtDSS.txt")
Input_files_CellPhoneDB_generation(merged_ann2, "mtDSS-TiO2", "cellphonedb_count_mtTiO2.txt","cellphonedb_meta_mtTiO2.txt")
Input_files_CellPhoneDB_generation(merged_ann2, "wtDSS-TiO2", "cellphonedb_count_wtTiO2.txt","cellphonedb_meta_wtTiO2.txt")


