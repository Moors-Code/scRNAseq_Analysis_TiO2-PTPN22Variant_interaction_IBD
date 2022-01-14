##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

###### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### load R object after preprocessing, clustering and annotation ##### 
merged_ann <- readRDS(file = "merged_ann.Rds")

##### Re-clustering of CD4+ T cells #####
Idents(merged_ann) <- "annotation"
cd4 <- subset(merged_ann, idents = "CD4+ T")
cd4 <- SCTransform(cd4, vars.to.regress = c("percent.mito", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
cd4 <- RunPCA(object = cd4, features = VariableFeatures(object =cd4), npcs = 20, verbose = FALSE)
print(ElbowPlot(cd4))
cd4 <- FindNeighbors(object = cd4, dims = 1:15)
cd4 <- FindClusters(cd4, resolution = 0.1, random.seed = 5, algorithm = 2)
cd4 <- RunUMAP(cd4, dims = 1:15, seed.use = 5)

#plot umap space per condition to see differences in cluster abundance
p <- DimPlot(cd4, reduction = "umap", label = TRUE, split.by = "condition", pt.size = 3)
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("CD4+ T cells across conditions") + ggsave("CD4_T_cell_subsets_across_conditions.pdf", width = 10, height = 5)

#DGE analysis of differentially expressed genes per cluster
FindMarker_plus_heatmap_plotting(cd4,"seurat_clusters","Top 10 marker genes (DGE) CD4 T cells subclusters", 
                                 "Cluster_marker_CD4_T_cell_subclusters.pdf")

##### Annotation of CD4+ T cell subclusters #####
#T regulatory cells (Treg): Tnfrsf18, Foxp3, Cd4, Icos,Ccr5 -->  Cluster  0,3
#T helper cells: Bcl6, Il6ra, Cd3d,Stat6,Ahr -->  Cluster 1,2
p <- DotPlot(cd4, features = c("Cd4", "Tnfrsf18", "Foxp3","Ccr5","Icos",
                               "Bcl6","Il6ra","Cd3d", "Ahr","Stat6"))
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("CD4 + T cell subsets marker") + theme(axis.text.x = element_text(angle = 90)) + 
  ggsave("DotPlot_MarkerGenes_CD4_T_cell_subsets.pdf", width = 6, height = 5)

#add annotation of cells to Seurat object into a new meta.data column  
current.cluster.ids <- c(0,1,2,3)
new.cluster.ids <- c( "Treg",  "Treg","T helper cells", "T helper cells")
cd4$annotation <- plyr::mapvalues(x = cd4$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(cd4, group.by = "annotation")

#Dotplot of marker used for annotation in annotated Seurat object
Idents(cd4) <- "annotation"
p <- DotPlot(cd4, features = c("Cd4", "Tnfrsf18", "Foxp3","Ccr5","Icos",
                               "Bcl6","Il6ra","Cd3d", "Ahr","Stat6"))
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("CD4 + T cell subsets marker") + theme(axis.text.x = element_text(angle = 90)) + 
  ggsave("DotPlot_MarkerGenes_CD4_T_cell_subtypes.pdf", width = 6, height = 5)

#DGE analysis of top differentially expressed genes between annotated CD4+ T cell subsets 
FindMarker_plus_heatmap_plotting_and_csv(cd4,"annotation","Top 10 marker genes (DGE) CD4 T cells subtypes", 
                                         "Cluster_marker_CD4_T_cell_subtypes.csv","Cluster_marker_CD4_T_cell_subtypes.pdf")

##### Plot annotated umap plots for each condition separately #####
create_annotated_umap_plot_per_condition_noL_noCol(cd4, "wtDSS", "WT DSS","CD4_WT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd4, "wtDSS-TiO2", "WT DSS TiO2","CD4_WT_DSS_TiO2.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd4, "mtDSS", "R619W DSS","CD4_MT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd4, "mtDSS-TiO2", "R619W DSS TiO2","CD4_MT_DSS_TiO2.pdf")

##### Proportions of cell types per condition/sample #####
#these values were used for statistical analysis in Graphpad prism 
create_table_cell_type_prop(cd4, "CD4T_Proportions_cell_types_per_condition_sample.csv")


