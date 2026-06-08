#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(ggplot2)
library(RColorBrewer)
library(viridis)

set.seed(1234)

dir.create('plots')
dir.create('seurat_objects')

##
marrow_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/seurat_objects/merged_integ.anno.subset.rds')

celltype_list <- c('HSC','MPP_MyLy','MPP_MkEry','LMPP','MLP','CLP','PreProB','Cycling_ProB','VDJ_ProB','Large_PreB','Small_PreB','ImmatureB','MatureB')
marrow_obj <- subset(marrow_obj, subset=celltype %in% celltype_list)
marrow_obj$celltype <- factor(marrow_obj$celltype, levels=celltype_list)
print(table(marrow_obj$celltype))

marrow_obj[['activity']] <- NULL

#
tcell_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/seurat_objects/merged_rna.tcells.rds')

tcell_vec <- as.character(tcell_obj$celltype)
tcell_vec[tcell_vec %in% c('DN(Q)','DN(P)')] <- 'DN'
tcell_vec[tcell_vec %in% c('DP(Q)','DP(P)')] <- 'DP'

tcell_vec <- factor(tcell_vec, levels=c('DN','DP','abT(entry)','CD8aa','CD8+T','CD8+Tmem','CD4+T','CD4+Tmem','Treg','T(agonist)','rdT','NK','ILC3'))
tcell_obj$celltype <- tcell_vec

print(table(tcell_obj$celltype))
print(table(tcell_obj$celltype, tcell_obj$source))

## merge
print(identical(rownames(marrow_obj), rownames(tcell_obj)))

genes_list <- intersect(rownames(marrow_obj), rownames(tcell_obj))
print(length(genes_list))

marrow_obj <- subset(marrow_obj, features=genes_list)
tcell_obj <- subset(tcell_obj, features=genes_list)

marrow_obj$source <- 'marrow'

rna_obj <- merge(x=marrow_obj, y=tcell_obj, merge.data=F)
rna_obj[['RNA']] <- JoinLayers(rna_obj[['RNA']])

rna_obj@meta.data <- rna_obj@meta.data[,c('nCount_RNA','nFeature_RNA','sample','source','celltype')]
tmp_obj <- subset(tcell_obj, subset=source=='thymus')

rna_obj$sample <- factor(rna_obj$sample, levels=c(levels(marrow_obj$sample), names(table(tmp_obj$sample))[table(tmp_obj$sample) > 0]))
rna_obj$celltype <- factor(rna_obj$celltype, levels=c(levels(marrow_obj$celltype), levels(tcell_obj$celltype)))
rna_obj$source <- factor(rna_obj$source, levels=c('marrow','thymus'))

print(sum(table(rna_obj$celltype)))
print(table(rna_obj$celltype))
print(table(rna_obj$sample))
print(table(rna_obj$source))
print(table(rna_obj$celltype, rna_obj$source))

##=============================
## markers

tcell_marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/data/marker_genes.tcell_subsets.txt', header=T)
tcell_markers <- unique(tcell_marker_data$gene)
print(length(tcell_markers))

bcell_marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/data/marker_genes.bcell_subsets.txt', header=T)
bcell_markers <- unique(bcell_marker_data$gene)
print(length(bcell_markers))

features <- unique(c(tcell_markers, bcell_markers))
print(length(features))

## process
max_dim <- 40

rna_obj <- NormalizeData(rna_obj, normalization.method='LogNormalize', scale.factor=10000)
#rna_obj <- FindVariableFeatures(rna_obj, selection='vst', nfeatures=4000)
VariableFeatures(rna_obj) <- features
rna_obj <- ScaleData(rna_obj, features=VariableFeatures(rna_obj))
rna_obj <- RunPCA(rna_obj, features=VariableFeatures(rna_obj), npcs=max_dim, verbose=F)
rna_obj <- cluster_sim_spectrum(rna_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='spearman')
rna_obj <- RunUMAP(rna_obj, reduction='css', dims=1:ncol(Embeddings(rna_obj, 'css')), reduction.name='css_umap')

##=============================
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
celltype_colors <- colorRampPalette(brewer.pal(12, 'Paired'))(length(levels(rna_obj$celltype)))

## source
pdf('plots/lympho_rna.source.css_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(rna_obj, reduction='css_umap', group.by='source', cols=source_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_rna.source.css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='source', cols=source_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

## celltype
pdf('plots/lympho_rna.celltype.css_umap.pdf', width=6, height=8, pointsize=3)
DimPlot(rna_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/lympho_rna.celltype.css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

png('plots/lympho_rna.celltype.label.css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T) + my_theme + ggtitle('')
dev.off()

##=============================
## save

saveRDS(rna_obj, 'seurat_objects/integ_rna.lympho.rds')

##=============================



