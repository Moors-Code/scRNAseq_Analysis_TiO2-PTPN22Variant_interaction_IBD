##### PACKAGES #####
library(dplyr)
library(tidyverse)
library(devtools)
library(ggplot2)
library(stringdist)
library(BiocManager)
library(ShortRead)
library(Seurat)
library(biomaRt)
library(ggraph)
library(data.table)   
library(Matrix)
library(tidyr)
library(SingleR)
library(celldex)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(cowplot)


#### FUNCTIONS ####

###Function to set working directory relatively to where the scripts are
#from https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

###Function to read in data generated from BD Seven Bridges platform (by Simona Baghai)
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Gene, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Gene")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}

###Function to generate Seurat objects for each sample from the dgMatrix downloaded form Seven Bridges BD Platform (by Erkin Acar)
#takes path from dgMatrix file, project name and condition 
#gives seurat object with first quality cut-off, min.cells =3, min.features = 200
create_seurat_from_condition <- function(
  path_to_st_file,
  project,
  condition,
  min.cells = 3,
  min.features = 200
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                               project = project,
                               min.cells = min.cells,
                               min.features = min.features)
  condition_sample$condition <- condition
  
  return(condition_sample)
}


###Function to project the cell type from SingleR result to the umap space to identify which cluster represents which cell type 
project_annotation_to_umap <- function(cell.type, singleResult, seurat_object) {
  temp <- as.data.frame(singleResult[5])
  # singleResult[5] is extracting pruned.labels from data frame
  temp$cell <- rownames(temp)
  temp <- temp%>%filter(temp$pruned.labels %in% cell.type)
  temp <- temp$cell
  print(DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 10, pt.size = 2, cells.highlight = temp, sizes.highlight = 2) + NoLegend() + ggtitle(cell.type))
}


###Function to create annotated umap for each condition 
create_annotated_umap_plot_per_condition <- function(
  seurat_object,
  ident_condition,
  title,
  output_file_name,
  colours
) {
  Idents(seurat_object) <- "condition"
  condition_subset <- subset(seurat_object, idents = ident_condition)
  p <- DimPlot(condition_subset, label = TRUE, label.size = 3, pt.size = 0.2, group.by = "annotation",
               cols = colours)
  p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 30))+ theme(axis.text = element_text(size = 10)) +
    ggtitle(title) +
    ggsave(output_file_name, width = 6, height = 5)
}

###Function to create annotated umap for each condition with no legend and without specified colors  
create_annotated_umap_plot_per_condition_noL_noCol <- function(
  seurat_object,
  ident_condition,
  title,
  output_file_name
) {
  Idents(seurat_object) <- "condition"
  condition_subset <- subset(seurat_object, idents = ident_condition)
  p <- DimPlot(condition_subset, label = TRUE, label.size = 5, pt.size = 3, group.by = "annotation") + NoLegend()
  p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 20)) + theme(title = element_text(size = 30))+ theme(axis.text = element_text(size = 10)) +
    ggtitle(title) +
    ggsave(output_file_name, width = 6, height = 5)
}

###Function to create tables of cell type proportions per condition/samples
create_table_cell_type_prop <- function(
  seurat_object, 
  output_file_name_csv
) {
  ann_tab<- table(seurat_object$orig.ident, seurat_object$annotation)
  ann_tab <- cbind(ann_tab, Total = rowSums(ann_tab))
  ann_tab <- as.data.frame(ann_tab)
  ann_tab_pct = lapply(ann_tab[,], function(x) {
    x/ann_tab$Total})
  ann_tab_pct <- as.data.frame(ann_tab_pct)
  rownames(ann_tab_pct) <- rownames(ann_tab)
  write.csv(ann_tab_pct, output_file_name_csv)
}


###Function to find differentially expressed genes between conditions 
##saves output as csv files
DGE_between_cond_csv <- function(
  seurat_object, 
  ident1,
  ident2,
  output_file_name_csv1,
  output_file_name_csv2) {
  Idents(seurat_object) <- "condition"
  marker1 <- FindMarkers(object = seurat_object, ident.1 = ident1, ident.2 = ident2, 
                         min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(marker1, file = output_file_name_csv1)
  marker2 <- FindMarkers(object = seurat_object, ident.1 = ident2, ident.2 = ident1, 
                         min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(marker2, file = output_file_name_csv2)
}


###Function for GSEA from differentially expressed genes between different conditions 
#takes list of genes generated by function DGE_between_cond_csv as input
gsea_DGE_between_conditions <- function(
  path_to_csv_file,
  gsea_pathways,
  output_file_name
){
  dge <- read.csv(path_to_csv_file)
  #genes are ranked taking both p_val_adj and avg_log2FC into account 
  ranks <- dge  %>% na.omit()%>% mutate(ranking=-log10(p_val_adj)/sign(avg_log2FC))
  ranks <- ranks$ranking
  names(ranks) <- dge$X
  Hallmarks_dge <- fgsea(pathways = Hallmarks, 
                                               stats = ranks,
                                               minSize=10,
                                               maxSize=500,
                                               nperm=1000000)
  Hallmarks_dge %>% filter(abs(NES)>1.5 & pval<0.5)  #only when this applies it is significant, 
  #if no significant hallmarks csv file is empty 
  Hallmarks_dge<- apply(Hallmarks_dge,2,as.character)
  write.csv(Hallmarks_dge, file = output_file_name)
}


###Function for DGE analysis of Seurat clusters or annotated clusters (FindMarker function) 
#saves output in a csv file and plots top 10 genes per cluster in a heatmap 
FindMarker_plus_heatmap_plotting_and_csv <- function(
  seurat_object,
  ident_seurat_object,
  title,
  output_file_name1,
  output_file_name2
){
  Idents(seurat_object) <- ident_seurat_object
  #DGE analysis of differentially expressed genes per cluster
  markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv((markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)), output_file_name1)
  #extract top 10 genes per annotated cluster and plot in heatmap 
  top10_markers <- markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)
  #plot top 10 markers in heatmap 
  p <- DoHeatmap(seurat_object, features = top10_markers$gene, size = 2) +theme(text = element_text(size=5)) 
  p + theme(legend.title = element_text(size = 8), legend.text = element_text(size = 5)) + theme(title = element_text(size = 10))+ 
    ggtitle(title) + scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
    ggsave(output_file_name2, width = 6, height = 5)
}


###Function for generation of input files for CellPhoneDB for separate conditions 
#first gene symbols have to be converted from mouse to human because CellPhoneDB database only contains human L-R interactions   
#Conversion of mouse to human genes partially from: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md
#and: https://www.cellphonedb.org/faq-and-troubleshooting
#output are two text tiles: gene counts and cell annotations (meta)
Input_files_CellPhoneDB_generation <- function(
  seurat_object,
  ident_condition,
  output_file_name_counts,
  output_file_name_meta
){
  Idents(seurat_object) <- "condition"
  #subset seurat object based on condition
  subset_cond <- subset(seurat_object, idents = ident_condition)
  #generating counts file 
  # take raw data and normalize it
  count_raw_meta <- GetAssayData(object = subset_cond, slot = "counts")[,colnames(x = subset_cond)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(count_norm_meta) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
  print(head(genesV2))
  matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=F),]
  matrixB <- matrixA
  matrixB$gene <- genesV2$Gene.stable.ID
  rownames(matrixA) <- matrixB$gene
  #save count matrix as text file 
  write.table(matrixA, output_file_name_counts, sep='\t', quote=F, row.names = T)
  # generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(subset_cond@meta.data), subset_cond@meta.data[,'annotation', drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, output_file_name_meta, sep='\t', quote=F, row.names=F)
}




