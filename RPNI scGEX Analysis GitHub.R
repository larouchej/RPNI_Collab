# RPNI Analysis.R
# Jacqueline Larouche
# Winter 2023

# Workflow: 
# 1. Initialize variables
# 2. Load Data and Create Seurat Object
# 3. Integrate using Seurat
# 4. Celltype Annotation
# 5. Differential Celltype Abundance
# 6. Differential Gene Expression

#### 1. Initialize Environment & Variables ####
library(Seurat) #v3.2.2
library(dplyr) #v1.0.2
library(ggplot2) # v3.2.1
library(tidyverse) # v1.3.0
library(tibble) # v3.0.1
library(Matrix)
library(data.table)
library(RColorBrewer)
library(tibble)
library(cowplot)
library(reshape)
library(dittoSeq)
library(scales) # for hue_pal
library(NCmisc) #v1.1.5 for p.toZ
library(patchwork) # plot_layout
library(ggpubr)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/RPNI_Collab")

cells.combined <- readRDS(file = "ProcessedData/RPNI_merged_integrated_02062023.RDS")

#### 2. Load Data and Create Seurat Object ####
tibial.data <- Read10X_h5(paste0(getwd(), "/RawData/7438-JM-1/filtered_feature_bc_matrix.h5"))
femoral.data <- Read10X_h5(paste0(getwd(), "/RawData/7438-JM-2/filtered_feature_bc_matrix.h5"))

tib <- CreateSeuratObject(counts = tibial.data, project = "Tibial", min.cells = 3, min.features = 200)
fem <- CreateSeuratObject(counts = femoral.data, project = "Femoral", min.cells = 3, min.features = 200)

cells.combined <- merge(tib, fem)
table(cells.combined$orig.ident)

cells.combined <- PercentageFeatureSet(cells.combined, pattern = "^Mt", col.name = "percent.mt")
VlnPlot(cells.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cells.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cells.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(file="Plots/QC_violin.pdf", height = 4, width = 8)
plot1 + plot2
dev.off()


cells.combined <- subset(cells.combined, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)

cells.combined <- cells.combined %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
cells.combined <- FindClusters(cells.combined, resolution = 0.2)

DimPlot(cells.combined, group.by = c('orig.ident', 'seurat_clusters')) #Needs integration

#### 3. Integrate using Seurat ####
cells.list <- SplitObject(cells.combined, split.by = "orig.ident")
cells.list <- lapply(X = cells.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = cells.list)
anchors <- FindIntegrationAnchors(object.list = cells.list, anchor.features = features)
cells.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(cells.combined) <- "integrated"
cells.combined <- cells.combined %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

cells.combined <- FindClusters(cells.combined, resolution = 0.1)
DimPlot(cells.combined, group.by = c('orig.ident', 'seurat_clusters')) 
table(cells.combined$orig.ident, cells.combined$seurat_clusters)

pdf(file="Plots/UMAP_sample_cluster.pdf", height = 4, width = 10)
DimPlot(cells.combined, group.by = c('orig.ident', 'seurat_clusters')) 
dev.off()

saveRDS(cells.combined, "ProcessedData/RPNI_merged_integrated_02062023.RDS")

#### 4. Celltype Annotation ####
Idents(cells.combined) <- 'seurat_clusters'
DefaultAssay(cells.combined) <- 'RNA'
cluster.degs <- FindAllMarkers(cells.combined, logfc.threshold = 0.5, only.pos = TRUE)
write.table(cluster.degs, file = 'seurat_cluster_degs.csv', sep = ',')

top10 <- cluster.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file="Plots/Heatmap_cluster_degs.pdf", height = 12, width = 8.5)
DoHeatmap(cells.combined, features = top10$gene) + NoLegend()
dev.off()

Idents(cells.combined) <- 'seurat_clusters'
current.cluster.ids <- 0:13
celltype.cluster.ids <- c('FAP-A', 'Macrophage', 'TCell', 'FAP-B', 'Endothelial-A', 'FAP-C', 'Fibroblast', 'NMJ', 'FAP-D', 'MuSC', 'Neutrophil', 'Endothelial-B', 'BCell', 'MastCell')
cells.combined$celltype <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = celltype.cluster.ids)

DefaultAssay(cells.combined) <- "RNA"

pdf(file="Plots/Dotplot_celltype_markers.pdf", height = 5, width = 8)
DotPlot(cells.combined, features = c('Pdgfra', 'Col3a1', 'Cd68', 'Adgre1', 'Cd3e', 'Nkg7', 'Fabp4', 'Aspn', 'Acta2', 'Cadm2', 'Mpz', 'Pax7', 'Myod1', 'Clec4d', 'Csf3r', 'Mmrn1', 'Cd79a', 'Ighm', 'Ms4a2', 'Mcpt9'), group.by = 'celltype') + RotatedAxis() + coord_flip()
dev.off()

saveRDS(cells.combined, file = "ProcessedData/maresin_merged_sct_09282022.RDS")

pdf(file="Plots/UMAP_celltype.pdf", height = 5, width = 8)
DimPlot(cells.combined, group.by = 'celltype', reduction = "umap", label = TRUE, pt.size = 2, split.by = 'orig.ident') + NoLegend()
dev.off()

#### 5. Differential Celltype Abundance ####
# Fold change differences in abundance
t1 <- table(cells.combined$orig.ident, cells.combined$celltype)
t2 <- table(cells.combined$orig.ident)
t3 <- sweep(t1, 1, t2, FUN = '/')
fc <- (t3[1,]-t3[2,])/t3[2,]
fc
df <- data.frame(fc)
df$celltype <- rownames(df)

p1 <- ggplot(df, aes(x = celltype, y = fc, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Fold Change") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-2,6)) + 
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
p1

pdf("Plots/Barplot_abundance_fold_change.pdf", width = 4, height = 4)
p1
dev.off()
# Statistical Test shows significance between samples for all celltypes
for (i in 1:ncol(t1)){
  print(colnames(t1)[i])
  res <- prop.test(x = t1[,i], n = t2, correct = FALSE)
  print(res)
}

#### 6. Differential Gene Expression & GOTerm ####
DefaultAssay(cells.combined) <- "integrated"
Idents(cells.combined) <- 'orig.ident'
sample.degs <- FindAllMarkers(cells.combined, only.pos = TRUE)

top10 <- sample.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(file="Plots/Heatmap_sample_degs.pdf", height = 12, width = 8.5)
DoHeatmap(cells.combined, features = top10$gene) + NoLegend()
dev.off()

library(clusterProfiler)
library(org.Rn.eg.db)
sample.degs$entrez <- mapIds(org.Rn.eg.db, rownames(sample.degs),'ENTREZID', 'SYMBOL')
write.table(sample.degs, file = 'RPNI_sample_degs.csv', sep = ',')

## GOTerm analysis
go_tibial <- enrichGO(sample.degs$entrez[sample.degs$cluster == 'Tibial'],
                      OrgDb = "org.Rn.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)
go_femoral <- enrichGO(sample.degs$entrez[sample.degs$cluster == 'Femoral'],
                      OrgDb = "org.Rn.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)

pdf(file="Plots/Dotplot_enrichGO.pdf", width = 12, height = 5)
dotplot(go_tibial, showCategory = 20) + dotplot(go_femoral, showCategory = 20)
dev.off()
