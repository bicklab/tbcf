#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(simspec)
library(ggplot2)
library(RColorBrewer)
library(viridis)

set.seed(123)

dir.create('plots')
dir.create('data')

rna_obj <- readRDS('seurat_objects/thymus_rna.rds')

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

##=============================
## removing cell cycle, ribosome, and tcr genes (NOT USED)

#data1 <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/cellcycle_genes.txt')
#cellcycle_genes <- data1$V1

#
#mito_genes <- rownames(rna_obj)[grep('^MT-', (rownames(rna_obj)))]

#a <- rownames(rna_obj)[grep('^RPS', (rownames(rna_obj)))]
#b <- rownames(rna_obj)[grep('^RPL', (rownames(rna_obj)))]
#ribo_genes <- unique(c(a, b))

#tcr_genes <- rownames(rna_obj)[grep('^TR[AB][VDJ]|^IG[HKL][VDJC]', (rownames(rna_obj)))]

#
#excluded_genes <- unique(c(cellcycle_genes, mito_genes, ribo_genes, tcr_genes))
#print(length(excluded_genes))

##=============================
##
max_dim <- 40

rna_obj <- NormalizeData(rna_obj, normalization.method='LogNormalize', scale.factor=10000)
rna_obj <- FindVariableFeatures(rna_obj, selection='vst', nfeatures=4000)
rna_obj <- ScaleData(rna_obj, features=VariableFeatures(rna_obj))
rna_obj <- RunPCA(rna_obj, features=VariableFeatures(rna_obj), npcs=max_dim, verbose=F)
rna_obj <- cluster_sim_spectrum(rna_obj, label_tag='sample', spectrum_type='corr_kernel', corr_method='pearson')
rna_obj <- RunUMAP(rna_obj, reduction='css', dims=1:ncol(Embeddings(rna_obj, 'css')), reduction.name='css_umap')

options(future.globals.maxSize = 1000 * 1024^2 * 4) # maximum memory of 4GB 
rna_obj <- FindNeighbors(rna_obj, reduction='css', features=VariableFeatures(rna_obj), dims=1:ncol(Embeddings(rna_obj, 'css')))
rna_obj <- FindClusters(rna_obj, resolution=5)

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

## seurat cluster
png('plots/seurat_cluster.css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='seurat_clusters', label=T, raster=F) + my_theme + ggtitle('')
dev.off()

## celltype
celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(rna_obj$celltype)))

pdf('plots/celltye.css_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(rna_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/celltype.css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='celltype', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

##=============================
##
markers_list <- list(
	tcell=c('CD3D','CD3E','CD3G','CD8A','CD8B','CD4','TRAC','RORC'),
	bcell=c('CD19','VPREB1','MS4A1','IGHD','IGHA1','FCER2'), # pre/pro, naive, memory, plasma b cells

	dc=c('CLEC9A','CLEC10A','LAMP3','LILRA4'), 
	mono=c('S100A9'), 
	mac=c('C1QA','CD68','CD14'), 
	mast=c('TPSB2','TPSAB1'),
	mgk=c('ITGA2B','ITGB3','GP1BA'), 
	ery=c('GYPA','GATA1','HBM','ALAS2'),
	fibro = c('PDGFRA','COLEC11','DCN','COL1A1'),
	vsmc=c('RGS5','ACTA2','TAGLN'),
	endo=c('CDH5','CLDN5'),
	tec=c('EPCAM','KRT1','KRT8','KRT14')
)

## expression
#for (i in names(markers_list)){
#        for (gene in markers_list[[i]]){
#                print(paste0(i, '; ', gene))
#
#                p1 <- FeaturePlot(rna_obj, features=gene, reduction='css_umap', min.cutoff='q0', max.cutoff='q100', order=T, raster=F) + scale_color_viridis_c() + my_theme + ggtitle('')
#
#                png(paste0('plots/expression.', i, '.', gene, '.png'), width=6, height=6, res=1000, units='in')
#                plot(p1)
#                dev.off()
#        }
#}

## signature
nbin <- 100
nctrl <- 200
nsize <- 1000

module_obj <- AddModuleScore(rna_obj, features=markers_list, bin=bin, ctrl=nctrl, name=names(markers_list), slot='data')
idx <- 1
for (i in names(markers_list)){
        print(i)

        p1 <- FeaturePlot(module_obj, features=paste0(i, idx), reduction='css_umap', min.cutoff='q0', max.cutoff='q100', order=T, raster=F) + scale_color_viridis(option='rocket', direction=-1) + my_theme + ggtitle('')

        png(paste0('plots/signatures.', i, '.png'), width=7, height=7, res=1000, units='in')
        plot(p1)
        dev.off()

        idx <- idx+1
}

##=============================
## broad cluster label
label_vec <- as.character(rna_obj$celltype)
label_vec[label_vec %in% c('DN(Q)','DN(P)','DP(Q)','DP(P)','abT(entry)','CD8aa','CD8+T','CD8+Tmem','CD4+T','CD4+Tmem','Treg','T(agonist)','rdT','NK','ILC3')] <- 'Tcell'
label_vec[label_vec %in% c('aDC','pDC','DC')] <- 'DC'
label_vec[label_vec %in% c('Mono','Mac')] <- 'Mono_Mac'

rna_obj$br_cluster <- factor(label_vec, levels=c('Tcell','Bcell','DC','Mono_Mac','Mast','Ery','Fb','Endo','VSMC','TEC'))

##
celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(rna_obj$br_cluster)))

pdf('plots/br_cluster.css_umap.pdf', width=5, height=6, pointsize=3)
DimPlot(rna_obj, reduction='css_umap', group.by='br_cluster', cols=celltype_colors, raster=F, label=T, label.size=2) + my_theme + ggtitle('') + theme(legend.position='bottom')
dev.off()

png('plots/br_cluster.css_umap.png', width=7, height=7, res=1000, units='in')
DimPlot(rna_obj, reduction='css_umap', group.by='br_cluster', cols=celltype_colors, raster=F, label=F) + my_theme + ggtitle('')
dev.off()

## marker
markers_res <- FindAllMarkers(rna_obj, group.by='br_cluster', test.use='wilcox', min.pct = 0.1, logfc.threshold = 0.3, only.pos=T)
markers_res <- markers_res[which(markers_res$p_val_adj < 0.001),]
print(table(markers_res$cluster))
print(length(unique(markers_res$gene)))

write.table(markers_res, 'data/cluster_markers.txt', col.names=T, row.names=F, quote=F, sep='\t')

##=============================
## save
saveRDS(rna_obj, 'seurat_objects/thymus_rna.anno.rds')


