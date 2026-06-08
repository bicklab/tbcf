#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(rtracklayer)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(viridis)


set.seed(1234)

dir.create('seurat_objects')
dir.create('plots')

##
marrow_ref_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/seurat_objects/merged_integ.anno.subset.rds')
DefaultAssay(marrow_ref_obj) <- 'peaks'
marrow_ref_obj[['RNA']] <- NULL

celltype_list <- c('HSC','MPP_MyLy','MPP_MkEry','LMPP','MLP','CLP','PreProB','Cycling_ProB','VDJ_ProB','Large_PreB','Small_PreB','ImmatureB', 'MatureB')
marrow_ref_obj <- subset(marrow_ref_obj, subset=celltype %in% celltype_list)
marrow_ref_obj$celltype <- factor(marrow_ref_obj$celltype, levels=celltype_list)
print(table(marrow_ref_obj$celltype))

#
marrow_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/2_merge/seurat_objects/marrow_atac.rds')
DefaultAssay(marrow_obj) <- 'peaks'

valid_cells <- colnames(marrow_ref_obj)[colnames(marrow_ref_obj) %in% colnames(marrow_obj)]
marrow_obj <- subset(marrow_obj, cells=valid_cells)
marrow_obj$celltype <- marrow_ref_obj$celltype[colnames(marrow_obj)]
marrow_obj$source <- 'marrow'

##
tcell_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/2_merge/seurat_objects/merged_atac.tcells.anno.rds')
DefaultAssay(tcell_obj) <- 'peaks'

tcell_vec <- as.character(tcell_obj$celltype)
tcell_vec[tcell_vec %in% c('DN(Q)','DN(P)')] <- 'DN'
tcell_vec[tcell_vec %in% c('DP(Q)','DP(P)')] <- 'DP'

tcell_vec <- factor(tcell_vec, levels=c('DN','DP','abT(entry)','CD8aa','CD8+T','CD8+Tmem','CD4+T','CD4+Tmem','Treg','T(agonist)','rdT','NK','ILC3'))
tcell_obj$celltype <- tcell_vec

print(table(tcell_obj$celltype))
print(table(tcell_obj$celltype, tcell_obj$source))

## merge
print(identical(rownames(marrow_obj), rownames(tcell_obj)))
peaks_list <- intersect(rownames(marrow_obj), rownames(tcell_obj))
print(length(peaks_list))

DefaultAssay(marrow_obj) <- 'activity'
DefaultAssay(tcell_obj) <- 'activity'
print(identical(rownames(marrow_obj), rownames(tcell_obj)))
genes_list <- intersect(rownames(marrow_obj), rownames(tcell_obj))
print(length(genes_list))

atac_obj <- merge(x=marrow_obj, y=tcell_obj, merge.data=F)
atac_obj@meta.data <- atac_obj@meta.data[, c('atac_fragments','peak_region_fragments','pct_reads_in_peaks','TSS.enrichment','sample','celltype','source')]

## remove peaks with low read counts
DefaultAssay(atac_obj) <- 'peaks'
print(dim(atac_obj))

counts <- atac_obj@assays$peaks$counts
count_vec <- rowSums(counts)
active_peaks <- names(count_vec[count_vec >= 5])
print(length(active_peaks))

atac_obj[['peaks']] <- subset(atac_obj[['peaks']], features=active_peaks)

##
tmp_obj <- subset(tcell_obj, subset=source=='thymus')
atac_obj$sample <- factor(atac_obj$sample, levels=c(levels(marrow_obj$sample), names(table(tmp_obj$sample))[table(tmp_obj$sample) > 0]))
atac_obj$celltype <- factor(atac_obj$celltype, levels=c(levels(marrow_obj$celltype), levels(tcell_obj$celltype)))
atac_obj$source <- factor(atac_obj$source, levels=c('marrow','thymus'))

print(sum(table(atac_obj$celltype)))
print(table(atac_obj$celltype))
print(table(atac_obj$sample))
print(table(atac_obj$source))
print(table(atac_obj$celltype, atac_obj$source))

##=====================================
## get markers

tcell_marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/data/marker_genes.tcell_subsets.txt', header=T)
tcell_markers <- unique(tcell_marker_data$gene)
print(length(tcell_markers))

bcell_marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/data/marker_genes.bcell_subsets.txt', header=T)
bcell_markers <- unique(bcell_marker_data$gene)
print(length(bcell_markers))

marker_genes <- unique(c(tcell_markers, bcell_markers))
print(length(marker_genes))

features <- intersect(marker_genes, rownames(atac_obj@assays$activity@counts))
print(length(features))

##=====================================
## process

## peak by css
DefaultAssay(atac_obj) <- 'peaks'
max_dim <- 30

atac_obj <- RunTFIDF(atac_obj, assay='peaks')
atac_obj <- FindTopFeatures(atac_obj, min.cutoff='q80')
atac_obj <- RunSVD(atac_obj)

cor_vec1 <- data.frame(t(cor(atac_obj$atac_fragments, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec2 <- data.frame(t(cor(atac_obj$TSS.enrichment, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))
cor_vec3 <- data.frame(t(cor(atac_obj$pct_reads_in_peaks, Embeddings(atac_obj, 'lsi')[,1:ncol(Embeddings(atac_obj, 'lsi'))])))

max_cor <- 0.6
exc_dims <- c(rownames(cor_vec1)[abs(cor_vec1) > max_cor], rownames(cor_vec2)[abs(cor_vec2) > max_cor], rownames(cor_vec3)[abs(cor_vec3) > max_cor])
print(exc_dims)

atac_obj <- cluster_sim_spectrum(atac_obj, label_tag='sample', spectrum_type='corr_ztransform', corr_method='spearman', use_dr='lsi', dims_use=setdiff(1:ncol(Embeddings(atac_obj, 'lsi')), c(1)), reduction.key='peakcss_', reduction.name='peak_css')
atac_obj <- RunUMAP(atac_obj, reduction='peak_css', dims=1:ncol(Embeddings(atac_obj, 'peak_css')), reduction.name='peak_css_umap')

## gene by css
DefaultAssay(atac_obj) <- 'activity'
max_dim <- 30

atac_obj <- NormalizeData(atac_obj, normalization.method='LogNormalize', scale.factor=10000)
#atac_obj <- FindVariableFeatures(atac_obj, selection='vst', nfeatures=3000)
VariableFeatures(atac_obj) <- features
atac_obj <- ScaleData(atac_obj, features=VariableFeatures(atac_obj))
atac_obj <- RunPCA(atac_obj, features=VariableFeatures(atac_obj), npcs=max_dim, verbose=F)

cor_vec1 <- data.frame(t(cor(atac_obj$atac_fragments, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))
cor_vec2 <- data.frame(t(cor(atac_obj$TSS.enrichment, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))
cor_vec3 <- data.frame(t(cor(atac_obj$pct_reads_in_peaks, Embeddings(atac_obj, 'pca')[,1:ncol(Embeddings(atac_obj, 'pca'))])))

max_cor <- 0.6
exc_dims <- c(rownames(cor_vec1)[abs(cor_vec1) > max_cor], rownames(cor_vec2)[abs(cor_vec2) > max_cor], rownames(cor_vec3)[abs(cor_vec3) > max_cor])
print(exc_dims)

atac_obj <- cluster_sim_spectrum(atac_obj, label_tag='sample', spectrum_type='corr_ztransform', corr_method='spearman', use_dr='pca', dims_use=setdiff(1:ncol(Embeddings(atac_obj, 'pca')), c()), reduction.key='genecss_', reduction.name='gene_css')
atac_obj <- RunUMAP(atac_obj, reduction='gene_css', dims=1:ncol(Embeddings(atac_obj, 'gene_css')), reduction.name='gene_css_umap')

##
atac_obj <- FindMultiModalNeighbors(
        atac_obj,
        reduction.list = list('peak_css','gene_css'),
        dims.list = list(1:ncol(Embeddings(atac_obj, 'peak_css')), 1:ncol(Embeddings(atac_obj, 'gene_css'))),
        modality.weight.name = c('peak_weight', 'gene_weight')
)

atac_obj <- RunUMAP(atac_obj, nn.name='weighted.nn', reduction.name='wnn_umap', reduction.key='wnnUMAP_')

##=====================================
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

source_colors <- brewer.pal(12, 'Paired')[c(2,4)]
celltype_colors <- colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(atac_obj$celltype)))

## source
#
pdf('plots/lympho_atac.source.wnn_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(atac_obj, reduction='wnn_umap', group.by='source', cols=source_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_atac.source.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='source', cols=source_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

## celltype
#
pdf('plots/lympho_atac.celltype.peak_css_umap.pdf', width=6, height=8, pointsize=3)
DimPlot(atac_obj, reduction='peak_css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_atac.celltype.peak_css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='peak_css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

#
pdf('plots/lympho_atac.celltype.gene_css_umap.pdf', width=6, height=8, pointsize=3)
DimPlot(atac_obj, reduction='gene_css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_atac.celltype.gene_css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='gene_css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

#
pdf('plots/lympho_atac.celltype.wnn_umap.pdf', width=6, height=8, pointsize=3)
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_atac.celltype.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

png('plots/lympho_atac.celltype.label.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T) + my_theme + ggtitle('')
dev.off()

##=====================================
## save 

saveRDS(atac_obj, 'seurat_objects/integ_atac.lympho.rds')

##=====================================
