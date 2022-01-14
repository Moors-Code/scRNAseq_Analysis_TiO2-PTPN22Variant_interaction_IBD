##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

##### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### Create Seurat objects from each condition and sample #####
#wtDSS samples 
wtDSS_1 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data", "IBD_2_SampleTag06_mm_wtDSS_1", "IBD_2_SampleTag06_mm_wtDSS_1_Expression_Data.st"), 
  project = "wtDSS_1", condition = "wtDSS")

wtDSS_2 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag02_mm_wtDSS_2", "IBD_2_SampleTag02_mm_wtDSS_2_Expression_Data.st"), 
  project = "wtDSS_2", condition = "wtDSS")

wtDSS_3 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag11_mm_wtDSS_3", "IBD_2_SampleTag11_mm_wtDSS_3_Expression_Data.st"), 
  project = "wtDSS_3", condition = "wtDSS")

#wtDSS-TiO2 samples 
wtTiO2_1 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data", "IBD_2_SampleTag10_mm_wtTiO2_1", "IBD_2_SampleTag10_mm_wtTiO2_1_Expression_Data.st"), 
  project = "wtTiO2_1", condition = "wtDSS-TiO2")

wtTiO2_2 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag07_mm_wtTiO2_2", "IBD_2_SampleTag07_mm_wtTiO2_2_Expression_Data.st"), 
  project = "wtTiO2_2", condition = "wtDSS-TiO2")

wtTiO2_3 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag12_mm_wtTiO2_3", "IBD_2_SampleTag12_mm_wtTiO2_3_Expression_Data.st"), 
  project = "wtTiO2_3", condition = "wtDSS-TiO2")

#mtDSS samples 
mtDSS_1 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag03_mm_mtDSS_1", "IBD_2_SampleTag03_mm_mtDSS_1_Expression_Data.st"), 
  project = "mtDSS_1", condition = "mtDSS")

mtDSS_2 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag09_mm_mtDSS_2", "IBD_2_SampleTag09_mm_mtDSS_2_Expression_Data.st"), 
  project = "mtDSS_2", condition = "mtDSS")

mtDSS_3 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag05_mm_mtDSS_3", "IBD_2_SampleTag05_mm_mtDSS_3_Expression_Data.st"), 
  project = "mtDSS_3", condition = "mtDSS")

#mtDSS-TiO2 samples 
mtTiO2_1 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag04_mm_mtTiO2_1", "IBD_2_SampleTag04_mm_mtTiO2_1_Expression_Data.st"), 
  project = "mtTiO2_1", condition = "mtDSS-TiO2")

mtTiO2_2 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag01_mm_mtTiO2_2", "IBD_2_SampleTag01_mm_mtTiO2_2_Expression_Data.st"), 
  project = "mtTiO2_2", condition = "mtDSS-TiO2")

mtTiO2_3 <- create_seurat_from_condition(
  path_to_st_file = file.path("/mnt","khandler","RStudio_Data", "IBD_zurich_marlene", "Data","IBD_2_SampleTag08_mm_mtTiO2_3", "IBD_2_SampleTag08_mm_mtTiO2_3_Expression_Data.st"), 
  project = "mtTiO2_3", condition = "mtDSS-TiO2")

##### Merge Seurat objects #####
merged <- merge(wtDSS_1, c(wtDSS_2, wtDSS_3, wtTiO2_1, wtTiO2_2, wtTiO2_3, mtDSS_1, mtDSS_2, mtDSS_3, mtTiO2_1, mtTiO2_2, mtTiO2_3))
#control merged object if it contains all samples 
table(merged$orig.ident)

##### Quality checks and removal of low quality cells and doublets #####
#add meta data column with percentage of genes mapped to mitochondrial genes
mito.features <- grep(pattern = "^mt-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")

#QC measurements 
#set Ident to condition, this way differences in quality of cells between different conditions can be analysed
Idents(merged) <- "condition"
p1 <- VlnPlot(merged, feature = "percent.mito")
p2 <- VlnPlot(merged, feature = "nCount_RNA")
p3 <- VlnPlot(merged, feature = "nFeature_RNA")
p4 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mito")
p5 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
median_percent.mito <- median(merged@meta.data$percent.mito)
median_nCount_RNA <- median(merged@meta.data$nFeature_RNA)
median_nFeature_RNA <- median(merged@meta.data$nFeature_RNA)
plot_grid(p1, p2,p3,p4,p5, labels = c(paste("median=",median_percent.mito, sep = ""), paste("median=",median_nCount_RNA, sep = ""),
                                paste("median=",median_nFeature_RNA, sep = ""))) + ggtitle("Quality measurements") +
  ggsave("Quality_measures.pdf", width = 15, height = 10)

#subset data to remove cells with a higher percentage of mitochondrial genes (> 0.1) --> removes low quality cells 
#also remove cells with > 8000 Feature_RNA counts --> removes doublets
merged_sub <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 0.1 )

##### Clustering of merged Seurat object ##### 
merged_sub <- SCTransform(merged_sub, vars.to.regress = c("percent.mito", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
merged_sub <- RunPCA(object = merged_sub, features = VariableFeatures(object =merged_sub), npcs = 20, verbose = FALSE)
print(ElbowPlot(merged_sub)) 
merged_sub <- FindNeighbors(object = merged_sub, dims = 1:10)
#Resolution set to 0.8 which gives a good resolution of clustering 
#algorithm =2 --> Louvain clustering with multilevel refinement
merged_sub <- FindClusters(merged_sub, resolution = 0.8, random.seed = 5, algorithm = 2)
merged_sub <- RunUMAP(merged_sub, dims = 1:12, seed.use = 5)

#plot feature umap space
p <- DimPlot(merged_sub, reduction = "umap", label = TRUE)
p + theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Feature space") + ggsave("Feature_umap_space_merged_sample.pdf", width = 6, height = 5)

##### save clustered Seurat object to go back later #####
saveRDS(merged_sub, file = "merged_sub.Rds")

#clustered Seurat object is then used for annotation in R code 3_Annotation_and_cluster_comparison.R 





