#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)
library(viridis)


set.seed(123)


##
marrow_obj <- readRDS('seurat_objects/marrow_atac.rds')
marrow_obj[[colnames(marrow_obj@meta.data)[!colnames(marrow_obj@meta.data) %in% c('atac_fragments','peak_region_fragments','pct_reads_in_peaks','TSS.enrichment','sample')]]] <- NULL

##
ref_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/2_integ_css/seurat_objects/merged_integ.anno.clean.rds')
valid_cells <- colnames(ref_obj)[which(ref_obj$celltype %in% c('CD4T','CD8T','NK'))]
valid_cells <- colnames(marrow_obj)[colnames(marrow_obj) %in% valid_cells]

marrow_obj <- subset(marrow_obj, cells=valid_cells)


##
thymus_obj <- readRDS('seurat_objects/thumus_atac.tcell.rds')
thymus_obj$atac_fragments <- thymus_obj$passed_filters
thymus_obj[[c('passed_filters','predicted.id','prediction.score.Tcell','prediction.score.max')]] <- NULL

##
marrow_obj$source <- 'marrow'
thymus_obj$source <- 'thymus'
merged_obj <- merge(x=marrow_obj, y=thymus_obj, merge.data=F)

#
merged_obj$source <- factor(merged_obj$source, levels=c('marrow','thymus'))
merged_obj$sample <- factor(merged_obj$sample, levels=c(levels(marrow_obj$sample), levels(thymus_obj$sample)))

print(table(merged_obj$source))
print(table(merged_obj$sample, merged_obj$source))

##
saveRDS(merged_obj, 'seurat_objects/merged_atac.tcells.rds')



