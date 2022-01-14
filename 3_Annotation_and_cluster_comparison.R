##### set working directory relatively to where the scripts are #####
getCurrentFileLocation()
setwd(getCurrentFileLocation())

##### connect to the R code 1_Packages_and_functions.R #####
source('1_Packages_and_functions.R', echo=TRUE)

##### load R object after preprocessing and clustering #####
merged_sub <- readRDS(file = "merged_sub.Rds")


##### automatic immune cell annotation with SingleR #####
#ImmGenData from the Immunological Genome Project (Heng et al. 2008) used as a reference 
#This analysis is used for a broad annotation into main immune cell types

#load ImmGen reference data 
mouse.se <- celldex::ImmGenData()

#SingleR needs conversion of the Seurat object to a SingleCellExperiment
#before doing that run DietSeurat to only keep the umap graph in addition to the counts and data slots
merged_sub_diet <- DietSeurat(merged_sub)
results <- SingleR(test = as.SingleCellExperiment(merged_sub_diet), ref = mouse.se, labels = mouse.se$label.main)

#this plots all the cells with a cell type specific score in a heatmap 
plotScoreHeatmap(results)

#look at the heatmap and write a list with all the potential cell types (only take these, where heatmap gives high expression signal (yellow))
cell.types <- c("ILC", "NK cells", "NKT", "T cells", "Tgd", "Fibroblasts", 
                "Stromal cells", "Eosinophils", "Neutrophils", "Endothelial cells", "Mast cells", "B cells", 
                "DC", "Monocytes", "Macrophages")

#loop through the list of cell types defined before 
#and print umap space with cells expressing specific cell markers highlighted in red = broad annotation
lapply(cell.types, function(x) project_annotation_to_umap(x, results, merged_sub))
#Cluster 12 = Macrophages/Monocytes, Cluster 10,13 = DCs, Cluster 0,6,4,16,11 = B cells, Cluster 17 = Endothelial cells 
#Cluster 14 = Neutrophils, Cluster 1 = Eosinophils, Cluster 7,15 = Fibroblasts, Cluster 9,3,2 = T cells,  
#Cluster 5,8 = ILCs

##### Fine annotation of immune cell subtypes #####

###DGE analysis of genes differentially expressed in each cluster
markers <- FindAllMarkers(object = merged_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#save top 10 DEGs each cluster in an excel sheet 
write.csv((markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)), "Marker_genes_clustering_not_annotated.csv")

###known marker gene analysis to define cell subtypes 
#gene expression of known marker genes used from https://www.rndsystems.com/resources/intestinal-lamina-propria-cell-markers#intraepithelial

#DCs (cluster 10,13)
#CD103+DC: H2-Ab1, H2-Eb1, Itgae+, Itgam- --> cluster13
#CD11b+DC: Itgam+, Cxcr1+ = cluster10
DotPlot(merged_sub, features = c("Cx3cr1", "Itgae", "Itgam", "Itgax", "H2-Ab1", "H2-Eb1"))

#Monocytes/Macrophages (cluster 12)
#Cd11b = Itgam,Adgre1 = F4/80 
DotPlot(merged_sub, features = c("Lyz2", "Itgam", "Cx3cr1", "Adgre1", "Il10"))

#Fibroblasts (Cluster 7,15)
DotPlot(merged_sub, features = c("Col3a1", "Dcn", "Igfbp4", "Ereg", "Adamdec1"))

#Endothelial cells (cluster 17)
DotPlot(merged_sub, features = c("Apold1", "Fabp4", "Timp3"))

#B cells (Cluster 0,6,4,11,16)
#Immature B cells: Cd19 high, Ms4a1 --> cluster 0,6
#B plasma cells: Cd38 high, Cd138 = Syndecan1+, Cd19 low, Igkc --> cluster 4,11,16
DotPlot(merged_sub, features = c("Cd19", "Ms4a1", "Igha", "Sdc1", "Cd38", "Cd5", "Cr2", "Il10", "Igkc"))

#T cells (Cluster 3,9,2)
#CD8 T cells: Cd8a  --> cluster 3,9
#CD4 T cells: --> cluster 2
DotPlot(merged_sub, features = c("Cd3d","Cd4", "Bcl6", "Cd8a", "Icos", "Pdcd1", "Il2", "Cxcr5", "Foxp3", "Tnfrsf18"))

#ILCs (Cluster 5,12,14)
#ILC3: Rorc, Il23r --> cluster 5
#ILC2: Gata3, Il2r --> cluster 8
DotPlot(merged_sub, features = c("Il7r", "Il2ra", "Gata3", "Klrg1", "Il23r", "Rorc"))

#Eosinophils (cluster 1)
DotPlot(merged_sub, features = c("Ccr3", "Cd69", "Itga4", "Itgam", "Retnlg"))

#Neutrophils (cluster 14)
DotPlot(merged_sub, features = c("Sell", "S100a9", "S100a8", "Cxcl2"))

#add annotation of cells to Seurat object into a new meta.data column  
Idents(merged_sub) <- "seurat_clusters"
current.cluster.ids <- c(0:17)
new.cluster.ids <- c("B", "Eosinophils", "CD4+ T", "CD8+ T", "Plasma cells", "ILC3", "B","Fibroblasts", "ILC2", "CD8+ T", 
                     "CD11b+ DC", "Plasma cells", "Macrophages","CD103+ DC", "Neutrophils", "Fibroblasts", "Plasma cells", "Endothelial")
merged_sub$annotation <- plyr::mapvalues(x = merged_sub$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

merged_ann <- merged_sub

#colors used for cell clusters: 
#B cells: Plasma cells: #507E9B, B: #5BCAD8             
#T cells: CD4: #EF6572, CD8: #EF65B1  
#ILC2: #8E5E78, ILC3: #EA8695  
#Eosinophils: #878481
#Stromal cells: Fibroblasts: #AA62C4, Endothelial: #C782F4
#Macrophages: #E0872F
#Neutrophils: #E2E20C
#DC: CD11b: #71BA56, CD103: #81DD17

p <- DimPlot(merged_ann, label = TRUE, label.size = 3, pt.size = 0.2, group.by = "annotation",
             cols = c("#5BCAD8", "#878481", "#EF6572", "#EF65B1","#507E9B", "#EA8695", "#AA62C4", 
                      "#8E5E78","#71BA56","#E0872F", "#81DD17", "#E2E20C", "#C782F4"))
p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 30))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Cell type annotation") +
  ggsave("Annotated_umap.pdf", width = 6, height = 5)

##### Plot annotated umap plots for each condition separately #####
create_annotated_umap_plot_per_condition(merged_ann, "wtDSS", "WT DSS","WT_DSS.pdf",
                                         c("#5BCAD8", "#878481", "#EF6572", "#EF65B1","#507E9B", "#EA8695", "#AA62C4", 
                                           "#8E5E78","#71BA56","#E0872F", "#81DD17", "#E2E20C", "#C782F4"))
create_annotated_umap_plot_per_condition(merged_ann, "wtDSS-TiO2", "WT DSS TiO2","WT_DSS_TiO2.pdf",
                                         c("#5BCAD8", "#878481", "#EF6572", "#EF65B1","#507E9B", "#EA8695", "#AA62C4", 
                                           "#8E5E78","#71BA56","#E0872F", "#81DD17", "#E2E20C", "#C782F4"))
create_annotated_umap_plot_per_condition(merged_ann, "mtDSS", "R619W DSS","MT_DSS.pdf",
                                         c("#5BCAD8", "#878481", "#EF6572", "#EF65B1","#507E9B", "#EA8695", "#AA62C4", 
                                           "#8E5E78","#71BA56","#E0872F", "#81DD17", "#E2E20C", "#C782F4"))
create_annotated_umap_plot_per_condition(merged_ann, "mtDSS-TiO2", "R619W DSS TiO2","MT_DSS_TiO2.pdf",
                                         c("#5BCAD8", "#878481", "#EF6572", "#EF65B1","#507E9B", "#EA8695", "#AA62C4", 
                                           "#8E5E78","#71BA56","#E0872F", "#81DD17", "#E2E20C", "#C782F4"))

##### Plot top 10 marker genes per annotated cluster in heatmap ##### 
#DGE analysis of annotated cluster
FindMarker_plus_heatmap_plotting_and_csv(merged_ann,"annotation","Top 10 marker genes (DGE)", 
                                 "Annotated_cluster_top10_marker_genes.csv","Annotated_cluster_top10_marker_genes.pdf")

##### Proportions of cell types per condition/sample #####
#these values were used for statistical analysis in Graphpad prism 
create_table_cell_type_prop(merged_ann, "Proportions_cell_types_per_condition_sample.csv")

##### save annotated Seurat object to go back later #####
saveRDS(merged_ann, file = "merged_ann.Rds")

#clustered and annotated object is then used for further analysis

