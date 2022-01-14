##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

##### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### load R object after preprocessing, clustering and annotation #####
merged_ann <- readRDS(file = "merged_ann.Rds")

##### Re-clustering of CD8+ T cells #####
Idents(merged_ann) <- "annotation"
cd8 <- subset(merged_ann, idents = "CD8+ T")
cd8 <- SCTransform(cd8, vars.to.regress = c("percent.mito", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
cd8 <- RunPCA(object = cd8, features = VariableFeatures(object =cd8), npcs = 20, verbose = FALSE)
print(ElbowPlot(cd8))
cd8 <- FindNeighbors(object = cd8, dims = 1:10)
cd8 <- FindClusters(cd8, resolution = 0.1, random.seed = 5, algorithm = 2)
cd8 <- RunUMAP(cd8, dims = 1:10, seed.use = 5)

#plot umap space per condition to see differences in cluster abundance
p <- DimPlot(cd8, reduction = "umap", label = TRUE, split.by = "condition", pt.size = 3)
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("CD8+ T cells across conditions") + ggsave("CD8_T_cell_subsets_across_conditions.pdf", width = 10, height = 5)

#DGE analysis per cluster
FindMarker_plus_heatmap_plotting_and_csv(cd8,"seurat_clusters","Top 10 marker genes (DGE) CD8 T cells subclusters", 
                                         "Cluster_marker_CD8_T_cell_subclusters.csv","Cluster_marker_CD8_T_cell_subclusters.pdf")

##### Annotation of CD8+ T cell subclusters #####
#cytotoxic T cells: GzmB, Ifng, Tbet = Tbx21, Klrg1, Tim3 = Havcr2, Blimp1 = Prdm1, Zeb2, Irf4, Batf --> Cluster 0,2 (Cluster 2 is more in G2M/S --> see analysis underneath)
#memory T cells: Ccr7, Cd62l = Sell, Tcf1 = Tcf7, Nsg2, Xcl1, Bcl6, ID2, Nr4a1 --> Cluster 1
DotPlot(cd8, features = c("Gzmb", "Ifng", "Tbx21", "Klrg1", "Havcr2", "Prdm1", "Zeb2", "Irf4", "Batf", 
                          "Ccr7", "Sell", "Tcf7","Nsg2", "Xcl1", "Bcl6","Id2","Nr4a1","Il10")) 

#Visualize cell cycle of cells in umap space depending on the expression of know cell cycle marker genes
#load cell cycle gene list 
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
cyclingCD8T <-CellCycleScoring(object = cd8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
p <- DimPlot(cyclingCD8T, pt.size =3)
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 15))+ theme(axis.text = element_text(size = 8)) +
  ggtitle("Cycling score in cells of CD8+ T cells") + ggsave("Cell_cycle_umap.pdf", width = 6, height = 5)

#add annotation of cells to Seurat object into a new meta.data column  
Idents(cd8) <- "seurat_clusters"
current.cluster.ids <- c(0,1,2)
new.cluster.ids <- c("CD8+ Cytotoxic T-cells", "CD8+ Memory T-cells","CD8+ Cytotoxic T-cells" )
cd8$annotation <- plyr::mapvalues(x = cd8$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#DGE analysis of top differentially expressed genes between annotated CD8+ T cell subsets 
FindMarker_plus_heatmap_plotting_and_csv(cd8,"annotation","Top 10 marker genes (DGE) CD8 T cells subctypes", 
                                         "Cluster_marker_CD8_T_cell_subtypes.csv","Cluster_marker_CD8_T_cell_subtypes.pdf")

#DotPlot of marker genes used for annotation 
Idents(cd8) <- "annotation"
p <- DotPlot(cd8, features = c("Gzmb", "Ifng", "Tbx21", "Klrg1", "Havcr2", "Prdm1", "Zeb2", "Irf4", "Batf", 
                               "Ccr7", "Sell", "Tcf7","Nsg2", "Xcl1", "Bcl6","Id2","Nr4a1","Il10")) 
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("CD8 + T cell subsets marker") + theme(axis.text.x = element_text(angle = 90)) + 
  ggsave("DotPlot_MarkerGenes_CD8_T_cell_subtypes.pdf", width = 6, height = 5)

##### Plot annotated umap plots for each condition separately #####
create_annotated_umap_plot_per_condition_noL_noCol(cd8, "wtDSS", "WT DSS","CD8_T_WT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd8, "wtDSS-TiO2", "WT DSS TiO2","CD8_WT_DSS_TiO2.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd8, "mtDSS", "R619W DSS","CD8_MT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(cd8, "mtDSS-TiO2", "R619W DSS TiO2","CD8_MT_DSS_TiO2.pdf")

##### Proportions of cell types per condition/sample #####
#these values were used for statistical analysis in Graphpad prism 
create_table_cell_type_prop(cd8, "CD8T_Proportions_cell_types_per_condition_sample.csv")

##### DGE analysis between conditions in cytotoxic CD8+ T cells (most difference in cell abundance found) #####
#Subset cytotoxic CD8+ T cells
Idents(cd8) <- "annotation"
cd8_sub <- subset(cd8, idents = "CD8+ Cytotoxic T-cells")

#DGE MT DSS vs. MT DSS TiO2
DGE_between_cond_csv(cd8_sub, "mtDSS","mtDSS-TiO2",
                     "cd8_clusterCytotoxixT_DGE_mtDSS_vs_mtTiO2.csv","cd8_clusterCytotoxicT_DGE_mtTiO2_vs_mtDSS.csv")

#DGE WT DSS vs. MT DSS
DGE_between_cond_csv(cd8_sub, "wtDSS","mtDSS",
                     "cd8_clusterCytotoxicT_DGE_wtDSS_vs_mtDSS.csv","cd8_clusterCytotoxicT_DGE_mtDSS_vs_wtDSS.csv")

#plot average expression of significant genes that are present in all 4 comparisons 
#(overlapping genes from Venn diagram analysis of above generated DEG lists done by Yasser Morsy)
sig_genes <- c("Mmd","Fabp5","Snx9","Ddit4","Ets2","Cpd","Tnfrsf4","Tmem64","Foxp3",
               "Smco4","Dusp4","Il10","Il2ra","Tnfrsf9",
               "Odc1","Capg","Havcr2","Vim","Socs2","Zfp36l1","Maf","Gem","Nfil3","Samsn1","S100a4","Ass1",
               "H2-K1","Nkg7","Pde3b","Cd8a","Plac8","Ccl5","B2m","Gzmk","Gramd3","Ptprc","Gimap6","Sh2d2a","H2-Q7","Ms4a4b","Trbc2")

cd8_sub$condition_reorder <- factor(x = cd8_sub$condition, levels = c("mtDSS","wtDSS", "mtDSS-TiO2", "wtDSS-TiO2"))
Idents(cd8_sub) <- "condition_reorder"
av_cd8_cytoT <- AverageExpression(cd8_sub, return.seurat = TRUE)
a <- DoHeatmap(av_cd8_cytoT, features = sig_genes)  
a + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 10))+ 
  theme(axis.text = element_text(size = 10)) + ggtitle("40 sig differentially expressed genes CD8+ cytotoxic T cells")  +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggsave("DGE_40sig_cytoT_averageExpression_Uniform_colours.pdf", width = 6, height = 5)

##### Analysis of Infg expression in CD8+ cytotoxic T cells across different conditions ##### 
#slot counts is the raw data 
#slot data is the log normalized data if LogNormalization was done beforehand
#SCT is scaled and normalized 
#log normalized values are the best option for proper comparison 
#normalize data before comparing the average expression, because slot data is not normalized 
#these values were used for statistical analysis in Graphpad prism
cd8_sub_norm <- NormalizeData(cd8_sub, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
average_exp_infg_cd8cyto <- AverageExpression(cd8_sub_norm, features = "Ifng", assays = "RNA", slot = "data", group.by = "orig.ident") #log normalized 
write.csv(average_exp_infg_cd8cyto, "average_exp_infg_cd8cyto.csv")
#AverageExpression(cd8_sub_norm, features = "Ifng", assays = "RNA", slot = "counts", group.by = "orig.ident") #raw data
#AverageExpression(cd8_sub_norm, features = "Ifng", assays = "SCT", slot = "data", group.by = "orig.ident") #SCT normalized

##### GSEA #####
#takes lists of differentially expressed genes from DGE analysis between conditions in cytotoxic CD8+ T cells as input 
#the Hallmarks gene set was used for analysis 
#Enrichment scores were considered significant with p-value adjusted ≤0.5
m_df      <- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#mtDSS vs. mtTiO2
#no significant Hallmarks
gsea_DGE_between_conditions("cd8_clusterCytotoxixT_DGE_mtDSS_vs_mtTiO2.csv",
                            Hallmarks,"Hallmarks_cd8_cytoT_DGE_mtDSS_vs_mtTiO2.csv")

#mtTiO2 vs. mtDSS
#no significant Hallmarks
gsea_DGE_between_conditions("cd8_clusterCytotoxicT_DGE_mtTiO2_vs_mtDSS.csv",
                            Hallmarks,"Hallmarks_cd8_cytoT_DGE_mtTiO2_vs_mtDSS.csv")

#mtDSS vs. wtDSS
gsea_DGE_between_conditions("cd8_clusterCytotoxicT_DGE_mtDSS_vs_wtDSS.csv",
                            Hallmarks,"Hallmarks_cd8_cytoT_DGE_mtDSS_vs_wtDSS.csv")

#wtDSS vs. mtDSS
gsea_DGE_between_conditions("cd8_clusterCytotoxicT_DGE_wtDSS_vs_mtDSS.csv",
                            Hallmarks,"Hallmarks_cd8_cytoT_DGE_wtDSS_vs_mtDSS.csv")

##### save re-clusterd CD8 T cell object #####
saveRDS(cd8, file = "CD8_seurat_object_reclustered.Rds")

