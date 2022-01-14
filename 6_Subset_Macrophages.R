##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

##### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### load R object after preprocessing, clustering and annotation #####
merged_ann <- readRDS(file = "merged_ann.Rds")

##### Re-clustering of Macrophages #####
Idents(merged_ann) <- "annotation"
macs <- subset(merged_ann, idents = "Macrophages")
macs <- SCTransform(macs, vars.to.regress = c("percent.mito", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
macs <- RunPCA(object = macs, features = VariableFeatures(object =macs), npcs = 20, verbose = FALSE)
print(ElbowPlot(macs))
macs <- FindNeighbors(object = macs, dims = 1:10)
macs <- FindClusters(macs, resolution = 0.1, random.seed = 2, algorithm = 2)
macs <- RunUMAP(macs, dims = 1:10, seed.use = 5)

#plot umap space per condition to see differences in cluster abundance
p <- DimPlot(macs, reduction = "umap", label = TRUE, split.by = "condition", pt.size = 3)
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Macrophages across conditions") + ggsave("Macrophages_subsets_across_conditions.pdf", width = 10, height = 5)

#DGE analysis of differentially expressed genes per cluster
FindMarker_plus_heatmap_plotting_and_csv(macs,"seurat_clusters","Top 10 marker genes (DGE) Macrophages", 
                                         "Cluster_marker_Macrophages_subsets.csv","Cluster_marker_Macrophages_subsets.pdf")

##### Annotation of Macrophages subclusters #####
#pro-inflammatory macrophages: Ly6c = Ly6c1, Il1b, Il6, Tnf, Il12a, Il23=Il12b --> cluster 0
#Anti-inflammatory macrophages: MhcII = H2-Ab1, Il10, Septin1, Lsp1, Cd206 = Mrc1 --> cluster 1,2
Idents(macs) <- "seurat_clusters"
p <- DotPlot(macs, features = c("Ly6c1", "Il1b", "Il6", "Tnf", "Il12a", "Il12b",
                                "Septin1", "Lsp1", "Il10", "Mrc1", "H2-Ab1"))
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 20))+ 
  theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Macrophages for annotation of subclusters") + ggsave("DotPlot_Macs_subset_marker.pdf", width = 6, height = 5)

#add annotation of cells to Seurat object into a new meta.data column  
current.cluster.ids <- c(0,1,2)
new.cluster.ids <- c("pro-infl.", "anti-infl.", "anti-infl.")
macs$annotation <- plyr::mapvalues(x = macs$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#plot average expression per annotated cluster in a Dotplot 
Idents(macs) <- "annotation"
p <- DotPlot(macs, features = c("Ly6c1", "Il1b", "Il6", "Il12a", "Il12b", "Lsp1",
                                "Septin1", "Il10", "Mrc1", "H2-Ab1", "Tnf"))
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) +
  ggtitle("Macrophages subclusters marker genes") + ggsave("Dotplot_marker_Macs_anti_pro.pdf", width = 6, height = 5)

#DEG analysis of differentially expressed genes of annotated subclusters and plotting in heatmap 
FindMarker_plus_heatmap_plotting_and_csv(macs,"annotation","Top 10 marker genes (DGE) subclusters Macrophages", 
                                         "Macs_annotated_clusters_marker.csv","Macs_annotated_clusters_marker.pdf")

##### Plot annotated umap plots for each condition separately #####
create_annotated_umap_plot_per_condition_noL_noCol(macs, "wtDSS", "WT DSS","Macs_WT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(macs, "wtDSS-TiO2", "WT DSS TiO2","Macs_WT_DSS_TiO2.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(macs, "mtDSS", "R619W DSS","Macs_MT_DSS.pdf")
create_annotated_umap_plot_per_condition_noL_noCol(macs, "mtDSS-TiO2", "R619W DSS TiO2","Macs_MT_DSS_TiO2.pdf")

##### Proportions of cell types per condition/sample #####
#these values were used for statistical analysis in Graphpad prism 
create_table_cell_type_prop(macs, "Macs_Proportions_cell_types_per_condition_sample.csv")

##### DGE analysis between conditions in Macrophage subtypes #####
#Subset into two different objects 
Idents(macs) <- "annotation"
macs_pro <- subset(macs, idents = "pro-infl.")
Idents(macs) <- "annotation"
macs_anti <- subset(macs, idents = "anti-infl.")

###Pro inflammatory macrophages 
#DGE MT DSS vs. MT DSS TiO2
DGE_between_cond_csv(macs_pro, "mtDSS","mtDSS-TiO2",
                     "Macs_pro_DGE_mtDSS_vs_mtTiO2.csv","Macs_pro_DGE_mtTiO2_vs_mtDSS.csv")

#DGE WT DSS vs. MT DSS
DGE_between_cond_csv(macs_pro, "wtDSS","mtDSS",
                     "Macs_pro_DGE_wtDSS_vs_mtDSS.csv","Macs_pro_DGE_mtDSS_vs_wtDSS.csv")

#plot average expression of significant genes (p_val_adj ≤0.05) in pro-inflammatory macrophages across different conditions
sig_genes_wtDSS_mtDSS <- c("Gbp2","Aif1","H2-Aa","H2-Eb1","H2-Ab1", "Cd74","Stat1","Rplp2","Ifitm6")
sig_genes_mtDSS_mtTiO2 <- c("CT010467.1","Gda","F5","Gbp2","Aif1","Hsp90b1","Stat1")
macs_pro$condition_reorder <- factor(x = macs_pro$condition, 
                                     levels = c("mtDSS","wtDSS", "mtDSS-TiO2", "wtDSS-TiO2"))
Idents(macs_pro) <- "condition_reorder"

#DEG genes wtDSS vs. mtDSS
av_macs_pro <- AverageExpression(macs_pro, return.seurat = TRUE)
p <- DoHeatmap(av_macs_pro, features = sig_genes_wtDSS_mtDSS)  
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 10))+ 
  theme(axis.text = element_text(size = 5)) + ggtitle("Sig DE genes pro-inflammatory macrophage wtDSS vs. mtDSS") +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggsave("DGE_significant_macs_pro_averageExpression_wtDSS_mtDSS.pdf", width = 6, height = 5)

#DEG genes mtDSS vs. mtTiO2
p <- DoHeatmap(av_macs_pro, features = sig_genes_mtDSS_mtTiO2)  
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 10))+ 
  theme(axis.text = element_text(size = 5)) + ggtitle("Sig DE genes pro-inflammatory macrophage mtDSS vs. mtTiO2") +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggsave("DGE_significant_macs_pro_averageExpression_mtDSS_mtTiO2.pdf", width = 6, height = 5)

###Anti inflammatory macrophages
#DGE MT DSS vs. MT DSS TiO2
DGE_between_cond_csv(macs_anti, "mtDSS","mtDSS-TiO2",
                     "Macs_anti_DGE_mtDSS_vs_mtTiO2.csv","Macs_anti_DGE_mtTiO2_vs_mtDSS.csv")

#DGE WT DSS vs. MT DSS
DGE_between_cond_csv(macs_anti, "wtDSS","mtDSS",
                     "Macs_anti_DGE_wtDSS_vs_mtDSS.csv","Macs_anti_DGE_mtDSS_vs_wtDSS.csv")

#plot average expression of significant genes (p_val_adj ≤0.05) in anti-inflammatory macrophages across different conditions
sig_genes2 <- c("B2m")
macs_anti$condition_reorder <- factor(x = macs_anti$condition, 
                                      levels = c("mtDSS","wtDSS", "mtDSS-TiO2", "wtDSS-TiO2"))
Idents(macs_anti) <- "condition_reorder"
av_macs_anti <- AverageExpression(macs_anti, return.seurat = TRUE)
p <- DoHeatmap(av_macs_anti, features = sig_genes2)  
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + theme(title = element_text(size = 10))+ 
  theme(axis.text = element_text(size = 5)) + ggtitle("Sig. differentially expressed genes anti-inflammatory macrophages") +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggsave("DGE_significant_macs_anti_averageExpression.pdf", width = 6, height = 5)


##### Analysis of Infgr1/Infgr2 expression in pro- and anti-inflammatory macrophages across different conditions #####  
#log normalized values are the best option for proper comparison 
#normalize counts first 
#these values were used for statistical analysis in Graphpad prism

#Pro-inflammatory macrophages
mac_pro_norm <- NormalizeData(macs_pro, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
Pro_av_exp_infgr1 <- as.data.frame(AverageExpression(mac_pro_norm, features = "Ifngr1", assays = "RNA", slot = "data", group.by = "orig.ident"))
rownames(Pro_av_exp_infgr1) <- "Infgr1"
Pro_av_exp_infgr2 <- as.data.frame(AverageExpression(mac_pro_norm, features = "Ifngr2", assays = "RNA", slot = "data", group.by = "orig.ident"))
rownames(Pro_av_exp_infgr2) <- "Infgr2"
Pro_av_exp_infgr <- rbind(Pro_av_exp_infgr1,Pro_av_exp_infgr2)
write.csv(Pro_av_exp_infgr, "Pro_infl_macs_infgr_expression.csv")

#anti-inflammatory macrophages 
macs_anti_norm <- NormalizeData(macs_anti, normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                margin = 1, assay = "RNA")
Anti_av_exp_infgr1 <- as.data.frame(AverageExpression(macs_anti_norm, features = "Ifngr1", assays = "RNA", slot = "data", group.by = "orig.ident"))
rownames(Anti_av_exp_infgr1) <- "Infgr1"
Anti_av_exp_infgr2 <- as.data.frame(AverageExpression(macs_anti_norm, features = "Ifngr2", assays = "RNA", slot = "data", group.by = "orig.ident"))
rownames(Anti_av_exp_infgr2) <- "Infgr2"
Anti_av_exp_infgr <- rbind(Pro_av_exp_infgr1,Pro_av_exp_infgr2)
write.csv(Anti_av_exp_infgr, "Anti_infl_macs_infgr_expression.csv")




