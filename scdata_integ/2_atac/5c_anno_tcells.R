#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R
  

library(Seurat)
library(Signac)
library(preprocessCore)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(viridis)


atac_obj <- readRDS('seurat_objects/merged_atac.tcells.init2.rds')

##========================================================
## assign celltype
# DN(Q) - 13, 20, 27
# abT(entry) - 9, 46, 77
# CD8aa - 10, 68, 74
# CD8+T - 4, 26, 30, 33, 50, 59, 80, 83
# CD8+Tmem - 1, 6, 18, 22, 36, 43, 49, 76, 81, 86, 95, 96, 97
# CD4+Tmem - 21, 51, 71, 72, 78, 85
# Treg - 7, 14, 48, 55, 84, 89
# T(agonist) - 42
# rdT - 16, 34, 90
# NK - 3, 32, 40, 41, 52, 53, 54, 92, 93
# ILC3 - 11, 23, 35

celltype_vec <- as.character(atac_obj$celltype)
celltype_vec[atac_obj$kmeans_cluster %in% c(13, 20, 27)] <- 'DN(Q)'
celltype_vec[atac_obj$kmeans_cluster %in% c(9, 46, 77)] <- 'abT(entry)'
celltype_vec[atac_obj$kmeans_cluster %in% c(10, 68, 74)] <- 'CD8aa'
celltype_vec[atac_obj$kmeans_cluster %in% c(4, 26, 30, 33, 50, 59, 80, 83)] <- 'CD8+T'
celltype_vec[atac_obj$kmeans_cluster %in% c(1, 6, 18, 22, 36, 43, 49, 76, 81, 86, 95, 96, 97)] <- 'CD8+Tmem'
celltype_vec[atac_obj$kmeans_cluster %in% c(21, 51, 71, 72, 78, 85)] <- 'CD4+Tmem'
celltype_vec[atac_obj$kmeans_cluster %in% c(7, 14, 48, 55, 84, 89)] <- 'Treg'
celltype_vec[atac_obj$kmeans_cluster %in% c(42)] <- 'T(agonist)'
celltype_vec[atac_obj$kmeans_cluster %in% c(16, 34, 90)] <- 'rdT'
celltype_vec[atac_obj$kmeans_cluster %in% c(3, 32, 40, 41, 52, 53, 54, 92, 93)] <- 'NK'
celltype_vec[atac_obj$kmeans_cluster %in% c(11, 23, 35)] <- 'ILC3'

celltype_vec <- factor(celltype_vec, levels=levels(atac_obj$celltype))
atac_obj$celltype <- celltype_vec
print(table(atac_obj$celltype))

##
my_theme <- theme(
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
)

celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(atac_obj$celltype)))

##
png('plots/tcell_atac.celltype.clean.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, label=T, raster=F) + my_theme + ggtitle('')
dev.off()

pdf('plots/tcell_atac.celltype.clean.wnn_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, label=T, raster=F) + my_theme + ggtitle('')
dev.off()

##========================================================
##

atac_obj[[c('peak_weight','gene_weight','kmeans_cluster')]] <- NULL
saveRDS(atac_obj, 'seurat_objects/merged_atac.tcells.anno.rds')

##========================================================


