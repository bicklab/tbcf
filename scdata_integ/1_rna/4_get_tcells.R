#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(ggplot2)
library(RColorBrewer)
library(viridis)


set.seed(123)

dir.create('data')

##
rna_obj <- readRDS('seurat_objects/marrow_rna.rds')

print(sum(table(rna_obj$celltypist_label)))
print(table(rna_obj$celltypist_label))

## curation
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='CD4+CTL')] <- 'CD4+T'
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='Th17')] <- 'CD4+T'
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='CD8αα(I)')] <- 'CD8aa(I)'
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='CD8αα(II)')] <- 'CD8aa(II)'
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='αβT(entry)')] <- 'abT(entry)'
rna_obj$celltypist_label[which(rna_obj$celltypist_label=='γδT')] <- 'rdT'

#
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('ETP','DN(early)'))] <- 'DN(Q)'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('CD8aa(I)','CD8aa(II)'))] <- 'CD8aa'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('Treg','Treg(diff)'))] <- 'Treg'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('B_pro/pre','B_naive','B_memory','B_plasma'))] <- 'Bcell'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('DC1','DC2'))] <- 'DC'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('aDC1','aDC2','aDC3'))] <- 'aDC'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('Fb_1','Fb_2'))] <- 'Fb'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('Endo','Lymph'))] <- 'Endo'
rna_obj$celltypist_label[which(rna_obj$celltypist_label %in% c('cTEC','mcTEC','mTEC(I)','mTEC(II)','mTEC(III)','TEC(neuro)','Epi_GCM2'))] <- 'TEC'

rna_obj <- subset(rna_obj, subset= celltypist_label %in% c('NMP', 'Tfh'), invert=T)

print(sum(table(rna_obj$celltypist_label)))
print(table(rna_obj$celltypist_label))

#
celltype_list <- c(
'DN(Q)','DN(P)','DP(Q)','DP(P)','abT(entry)',
'CD8aa','CD8+T','CD8+Tmem',
'CD4+T','CD4+Tmem',
'Treg','T(agonist)','rdT','NK','ILC3',
'Bcell','aDC','pDC','DC',
'Mono','Mac','Mast',
'Ery','Fb','Endo','VSMC','TEC'
)

rna_obj$celltype <- factor(rna_obj$celltypist_label, levels=celltype_list)
rna_obj[['celltypist_label']] <- NULL

#
print(sum(table(rna_obj$celltype)))
print(table(rna_obj$celltype))

tcell_types <- c('DN(Q)','DN(P)','DP(Q)','DP(P)','abT(entry)','CD8aa','CD8+T','CD8+Tmem','CD4+T','CD4+Tmem','Treg','T(agonist)','rdT','NK','ILC3')
rna_obj <- subset(rna_obj, subset=celltype %in% tcell_types)

print(table(rna_obj$celltype))

##
ref_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/seurat_objects/merged_integ.anno.clean.rds')
valid_cells <- colnames(ref_obj)[which(ref_obj$celltype %in% c('CD4T','CD8T','NK'))]
valid_cells <- colnames(rna_obj)[colnames(rna_obj) %in% valid_cells]

rna_obj <- subset(rna_obj, cells=valid_cells)
rna_obj$celltype <- factor(rna_obj$celltype, levels=tcell_types)

print(sum(table(rna_obj$celltype)))
print(table(rna_obj$celltype))

##
saveRDS(rna_obj, 'seurat_objects/marrow_rna.tcells.rds')

##======================================
## merge

thymus_obj <- readRDS('seurat_objects/thymus_rna.anno.rds')

print(table(thymus_obj$celltype))
thymus_obj[[c('RNA_snn_res.5','seurat_clusters','br_cluster')]] <- NULL

thymus_obj <- subset(thymus_obj, subset=celltype %in% tcell_types)
thymus_obj$celltype <- factor(thymus_obj$celltype, levels=tcell_types)
thymus_obj$source <- 'thymus'

#
rna_obj$source <- 'marrow'

#
merged_obj <- merge(x=rna_obj, y=thymus_obj, merge.data=F)
merged_obj[['RNA']] <- JoinLayers(merged_obj[['RNA']])

merged_obj$source <- factor(merged_obj$source, levels=c('marrow','thymus'))
merged_obj$celltype <- factor(merged_obj$celltype, levels=tcell_types)
merged_obj$sample <- factor(merged_obj$sample, levels=c(levels(rna_obj$sample), levels(thymus_obj$sample)))

print(table(merged_obj$source))
print(table(merged_obj$celltype, merged_obj$source))

##
saveRDS(merged_obj, 'seurat_objects/merged_rna.tcells.rds')

##======================================
## marker genes for tcell subsets

print(table(merged_obj$celltype))
merged_obj <- NormalizeData(merged_obj, normalization.method='LogNormalize', scale.factor=10000)

#
max_cell <- 1500
cells_list <- c()
for (celltype in levels(merged_obj$celltype)){
        print(celltype)

        tmp_cells <- colnames(merged_obj)[which(merged_obj$celltype==celltype)]
        if (length(tmp_cells) > max_cell){
                tmp_cells <- sample(tmp_cells, size=max_cell, replace=F)
        }
        cells_list <- c(cells_list, tmp_cells)
}

tmp_obj <- subset(merged_obj, cells=cells_list)
print(table(tmp_obj$celltype))

#
min_pct <- 0.1
min_q <- 0.001
min_log2fc <- 1

Idents(tmp_obj) <- 'celltype'
DefaultAssay(tmp_obj) <- 'RNA'

diff_data <- FindAllMarkers(object=tmp_obj, test.use='wilcox', only.pos=T, log2fc.threshold=min_log2fc, min.pct=min_pct)
diff_data <- diff_data[which(diff_data$p_val_adj < min_q & diff_data$avg_log2FC > min_log2fc), ]
print(table(diff_data$cluster))

marker_genes <- unique(diff_data$gene)
print(length(marker_genes))

write.table(diff_data, 'data/marker_genes.tcell_subsets.txt', col.names=T, row.names=F, quote=F, sep='\t')

##======================================
##
max_dim <- 40

merged_obj <- NormalizeData(merged_obj, normalization.method='LogNormalize', scale.factor=10000)
#merged_obj <- FindVariableFeatures(merged_obj, selection='vst', nfeatures=4000)
VariableFeatures(merged_obj) <- marker_genes
merged_obj <- ScaleData(merged_obj, features=VariableFeatures(merged_obj))
merged_obj <- RunPCA(merged_obj, features=VariableFeatures(merged_obj), npcs=max_dim, verbose=F)
merged_obj <- cluster_sim_spectrum(merged_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='pearson')
merged_obj <- RunUMAP(merged_obj, reduction='css', dims=1:ncol(Embeddings(merged_obj, 'css')), reduction.name='css_umap')

#
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

celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(merged_obj$celltype)))

pdf('plots/tcell_rna.celltype.css_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(merged_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/tcell_rna.celltype.css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(merged_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

png('plots/tcell_rna.celltype.label.css_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(merged_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T) + my_theme + ggtitle('')
dev.off()

##======================================



