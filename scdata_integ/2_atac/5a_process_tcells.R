#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(rtracklayer)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(viridis)


set.seed(123)


##=======================================
##
diff_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/data/marker_genes.tcell_subsets.txt', header=T)
marker_genes <- unique(diff_data$gene)
print(length(marker_genes))

##=======================================
##
atac_obj <- readRDS('seurat_objects/merged_atac.tcells.rds')

features <- intersect(marker_genes, rownames(atac_obj@assays$activity@counts))
print(length(features))

## peak by css
DefaultAssay(atac_obj) <- 'peaks'
max_dim <- 30

atac_obj <- RunTFIDF(atac_obj, assay='peaks')
atac_obj <- FindTopFeatures(atac_obj, min.cutoff='q90')
print(length(atac_obj@assays$peaks@var.features))
print(summary(atac_obj@assays$peaks@meta.features[atac_obj@assays$peaks@var.features,'count']))
atac_obj <- RunSVD(atac_obj)

cor_vec1 <- data.frame(t(cor(atac_obj$atac_fragments, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec2 <- data.frame(t(cor(atac_obj$TSS.enrichment, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec3 <- data.frame(t(cor(atac_obj$pct_reads_in_peaks, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))

max_cor <- 0.6
exc_dims <- c(rownames(cor_vec1)[abs(cor_vec1) > max_cor], rownames(cor_vec2)[abs(cor_vec2) > max_cor], rownames(cor_vec3)[abs(cor_vec3) > max_cor])
print(exc_dims)

atac_obj <- cluster_sim_spectrum(atac_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='spearman', use_dr='lsi', dims_use=setdiff(1:ncol(Embeddings(atac_obj, 'lsi')), c(1)), reduction.key='peakcss_', reduction.name='peak_css')
atac_obj <- RunUMAP(atac_obj, reduction='peak_css', dims=1:ncol(Embeddings(atac_obj, 'peak_css')), reduction.name='peak_css_umap')

## gene by css
features <- intersect(marker_genes, rownames(atac_obj@assays$activity@counts))
print(length(features))

DefaultAssay(atac_obj) <- 'activity'
max_dim <- 30

atac_obj <- NormalizeData(atac_obj, normalization.method='LogNormalize', scale.factor=10000)
VariableFeatures(atac_obj) <- features
atac_obj <- ScaleData(atac_obj, features=VariableFeatures(atac_obj))
atac_obj <- RunPCA(atac_obj, features=VariableFeatures(atac_obj), npcs=max_dim, verbose=F)

cor_vec1 <- data.frame(t(cor(atac_obj$atac_fragments, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))
cor_vec2 <- data.frame(t(cor(atac_obj$TSS.enrichment, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))
cor_vec3 <- data.frame(t(cor(atac_obj$pct_reads_in_peaks, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))

max_cor <- 0.6
exc_dims <- c(rownames(cor_vec1)[abs(cor_vec1) > max_cor], rownames(cor_vec2)[abs(cor_vec2) > max_cor], rownames(cor_vec3)[abs(cor_vec3) > max_cor])
print(exc_dims)

atac_obj <- cluster_sim_spectrum(atac_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='spearman', use_dr='pca', dims_use=setdiff(1:ncol(Embeddings(atac_obj, 'pca')), c(2)), reduction.key='genecss_', reduction.name='gene_css')
atac_obj <- RunUMAP(atac_obj, reduction='gene_css', dims=1:ncol(Embeddings(atac_obj, 'gene_css')), reduction.name='gene_css_umap')

##
atac_obj <- FindMultiModalNeighbors(
        atac_obj,
        reduction.list = list('peak_css','gene_css'),
        dims.list = list(1:ncol(Embeddings(atac_obj, 'peak_css')), 1:ncol(Embeddings(atac_obj, 'gene_css'))),
        modality.weight.name = c('peak_weight', 'gene_weight')
)

atac_obj <- RunUMAP(atac_obj, nn.name='weighted.nn', reduction.name='wnn_umap', reduction.key='wnnUMAP_')

##
nclusters <- 100
atac_obj$kmeans_cluster <- kmeans(atac_obj@reductions$wnn_umap@cell.embeddings, centers=nclusters)$cluster
print(table(atac_obj$kmeans_cluster))
print(summary(as.vector(table(atac_obj$kmeans_cluster))))

saveRDS(atac_obj, 'seurat_objects/merged_atac.tcells.init.rds')

##=======================================
## theme:
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

source_colors <- brewer.pal(12, 'Paired')[c(2,4)]

##
#
png('plots/tcell_atac.kmeans_cluster.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='kmeans_cluster', label=T, raster=F) + my_theme + ggtitle('')
dev.off()

pdf('plots/tcell_atac.kmeans_cluster.wnn_umap.pdf')
DimPlot(atac_obj, reduction='wnn_umap', group.by='kmeans_cluster', label=T, raster=F) + my_theme + ggtitle('')
dev.off()

## source
#
pdf('plots/tcell_atac.source.peak_css_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(atac_obj, reduction='peak_css_umap', group.by='source', cols=source_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/tcell_atac.source.peak_css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='peak_css_umap', group.by='source', cols=source_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

#
pdf('plots/tcell_atac.source.gene_css_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(atac_obj, reduction='gene_css_umap', group.by='source', cols=source_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/tcell_atac.source.gene_css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='gene_css_umap', group.by='source', cols=source_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

##=======================================



