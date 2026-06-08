#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(ggplot2)
library(RColorBrewer)
library(viridis)


dir.create('plots')

set.seed(123)

##============================================
##

rna_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/seurat_objects/thymus_rna.anno.rds')
marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/data/cluster_markers.txt', header=T)
marker_genes <- unique(marker_data$gene)

##
atac_obj <- readRDS('seurat_objects/thymus_atac.rds')
atac_obj@meta.data <- atac_obj@meta.data[, c('passed_filters','peak_region_fragments','pct_reads_in_peaks','TSS.enrichment','sample')]

## label transfer
features <- intersect(marker_genes, rownames(atac_obj@assays$activity@counts))
print(length(features))

transfer.anchors <- FindTransferAnchors(reference=rna_obj, query=atac_obj, reference.assay='RNA', query.assay='activity', features=features, reduction='cca', dims=1:30)
pred_res <- TransferData(anchorset=transfer.anchors, refdata=rna_obj$br_cluster, weight.reduction='cca', dims=1:30)

atac_obj <- AddMetaData(atac_obj, metadata=pred_res[,c('predicted.id','prediction.score.Tcell','prediction.score.max')])
atac_obj$predicted.id <- factor(atac_obj$predicted.id, levels=levels(rna_obj$br_cluster))
print(table(atac_obj$predicted.id))

##============================================
## process

## peak by css
DefaultAssay(atac_obj) <- 'peaks'
max_dim <- 30

atac_obj <- RunTFIDF(atac_obj, assay='peaks')
atac_obj <- FindTopFeatures(atac_obj, min.cutoff='q50')
atac_obj <- RunSVD(atac_obj)

cor_vec1 <- data.frame(t(cor(atac_obj$passed_filters, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec2 <- data.frame(t(cor(atac_obj$TSS.enrichment, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec3 <- data.frame(t(cor(atac_obj$pct_reads_in_peaks, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))

max_cor <- 0.6
exc_dims <- c(rownames(cor_vec1)[abs(cor_vec1) > max_cor], rownames(cor_vec2)[abs(cor_vec2) > max_cor], rownames(cor_vec3)[abs(cor_vec3) > max_cor])
print(exc_dims)

atac_obj <- cluster_sim_spectrum(atac_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='spearman', use_dr='lsi', dims_use=setdiff(1:ncol(Embeddings(atac_obj, 'lsi')), c(1,6,7)), reduction.key='peakcss_', reduction.name='peak_css')
atac_obj <- RunUMAP(atac_obj, reduction='peak_css', dims=1:ncol(Embeddings(atac_obj, 'peak_css')), reduction.name='peak_css_umap')

## gene by css
DefaultAssay(atac_obj) <- 'activity'
max_dim <- 30

atac_obj <- NormalizeData(atac_obj, normalization.method='LogNormalize', scale.factor=10000)
atac_obj <- FindVariableFeatures(atac_obj, selection='vst', nfeatures=4000)
atac_obj <- ScaleData(atac_obj, features=VariableFeatures(atac_obj))
atac_obj <- RunPCA(atac_obj, features=VariableFeatures(atac_obj), npcs=max_dim, verbose=F)

cor_vec1 <- data.frame(t(cor(atac_obj$passed_filters, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))
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

##============================================
##
my_theme <- theme(
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
)

## celltype
celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(levels(atac_obj$predicted.id)))

#
pdf('plots/thymus_atac.celltype.peak_css_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(atac_obj, group.by='predicted.id', reduction='peak_css_umap', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.celltype.peak_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='predicted.id', reduction='peak_css_umap', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()


#
pdf('plots/thymus_atac.celltype.gene_css_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(atac_obj, group.by='predicted.id', reduction='gene_css_umap', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.celltype.gene_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='predicted.id', reduction='gene_css_umap', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()


#
pdf('plots/thymus_atac.celltype.wnn_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(atac_obj, group.by='predicted.id', reduction='wnn_umap', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.celltype.wnn_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='predicted.id', reduction='wnn_umap', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

#
png('plots/thymus_atac.tcell_score.wnn_umap.png', width=7, height=7, res=1000, units='in')
FeaturePlot(atac_obj, features='prediction.score.Tcell', reduction='wnn_umap', min.cutoff='q0', max.cutoff='q100', order=T, raster=F)+ scale_color_viridis(option='rocket', direction=-1) + my_theme + ggtitle('')
dev.off()

png('plots/thymus_atac.pred_score.wnn_umap.png', width=7, height=7, res=1000, units='in')
FeaturePlot(atac_obj, features='prediction.score.max', reduction='wnn_umap', min.cutoff='q0', max.cutoff='q100', order=T, raster=F)+ scale_color_viridis(option='rocket', direction=-1) + my_theme + ggtitle('')
dev.off()


## tss enrich
atac_obj$tss_label <- ifelse(atac_obj$TSS.enrichment > 6, 'highQ', 'lowQ')
print(table(atac_obj$tss_label))

png('plots/thymus_atac.tss_enrich.peak_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='tss_label', reduction='peak_css_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.tss_enrich.gene_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='tss_label', reduction='gene_css_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.tss_enrich.wnn_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='tss_label', reduction='wnn_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

## pct peak
atac_obj$pctpeak_label <- ifelse(atac_obj$pct_reads_in_peaks > 40, 'highQ', 'lowQ')
print(table(atac_obj$pctpeak_label))

png('plots/thymus_atac.pctpeak_qc.peak_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='pctpeak_label', reduction='peak_css_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.pctpeak_qc.gene_css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='pctpeak_label', reduction='gene_css_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/thymus_atac.pctpeak_qc.wnn_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='pctpeak_label', reduction='wnn_umap', label=F, order=T, raster=F) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

##============================================
## select t cells

atac_obj <- FindClusters(atac_obj, graph.name='wsnn', resolution=2)

##
png('plots/thymus_atac.seurat_cluster.wnn_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='seurat_clusters', reduction='wnn_umap', order=T, raster=F, label=T) + my_theme + ggtitle('')
dev.off()

tcell_clusters <- c(0, 1, 2, 3, 5, 7, 8, 9, 10, 12, 13, 14, 15, 20, 22, 25, 28, 31, 34, 35, 36, 37, 38, 41)
atac_obj <- subset(atac_obj, subset=seurat_clusters %in% tcell_clusters)
atac_obj <- subset(atac_obj, subset=predicted.id=='Tcell')

#
png('plots/thymus_atac.tcell_cluster.wnn_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(atac_obj, group.by='seurat_clusters', reduction='wnn_umap', order=T, raster=F, label=T) + my_theme + ggtitle('')
dev.off()

##
atac_obj[[c('tss_label','pctpeak_label','peak_weight','gene_weight','wsnn_res.2','seurat_clusters')]] <- NULL
saveRDS(atac_obj, 'seurat_objects/thymus_atac.tcell.rds')

##============================================



