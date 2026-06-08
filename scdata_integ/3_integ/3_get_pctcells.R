#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)

dir.create('data')


rna_obj <- readRDS('seurat_objects/integ_rna.lympho.rds')
atac_obj <- readRDS('seurat_objects/integ_atac.lympho.rds')

##=============================
## pct cell rna

nbreaks <- 5
get_cellcount <- function(x){sum(x > 0)}

pctcell_mat <- NULL
for (celltype in levels(rna_obj$celltype)){
        print(celltype)

        tmp_cells <- colnames(rna_obj)[which(rna_obj$celltype==celltype)]
        subset_obj <- subset(rna_obj, cells=tmp_cells)
        tmp_breaks <- cut(1:length(tmp_cells), breaks=nbreaks)
	idx <- 1
        tmp_mat <- NULL
        for (i in levels(tmp_breaks)){
		print(idx)

                tmp_obj <- subset(subset_obj, cells=tmp_cells[which(tmp_breaks==i)])
                vec <- apply(tmp_obj@assays$RNA$counts, 1, get_cellcount)
                tmp_mat <- cbind(tmp_mat, vec)
		idx <- idx+1
        }

        #
        tmp_mat <- as.matrix(tmp_mat)
        cellcount_vec <- rowSums(tmp_mat)
        pctcell_vec <- cellcount_vec / table(rna_obj$celltype)[celltype]

        pctcell_mat <- cbind(pctcell_mat, pctcell_vec)
}

colnames(pctcell_mat) <- levels(rna_obj$celltype)
out_df <- data.frame(gene=rownames(pctcell_mat), pctcell_mat)
write.table(out_df, 'data/pctcells.rna.txt', row.names=F, col.names=T, quote=F, sep='\t')

##=============================
## pct cell atac

nbreaks <- 10
get_cellcount <- function(x){sum(x > 0)}

pctcell_mat <- NULL
for (celltype in levels(atac_obj$celltype)){
        print(celltype)

        tmp_cells <- colnames(atac_obj)[which(atac_obj$celltype==celltype)]
        subset_obj <- subset(atac_obj, cells=tmp_cells)
        tmp_breaks <- cut(1:length(tmp_cells), breaks=nbreaks)
	idx <- 1
        tmp_mat <- NULL
        for (i in levels(tmp_breaks)){
		print(idx)

                tmp_obj <- subset(subset_obj, cells=tmp_cells[which(tmp_breaks==i)])
                vec <- apply(tmp_obj@assays$peaks@counts, 1, get_cellcount)
                tmp_mat <- cbind(tmp_mat, vec)
		idx <- idx+1
        }

        #
        tmp_mat <- as.matrix(tmp_mat)
        cellcount_vec <- rowSums(tmp_mat)
        pctcell_vec <- cellcount_vec / table(atac_obj$celltype)[celltype]

        pctcell_mat <- cbind(pctcell_mat, pctcell_vec)
}

colnames(pctcell_mat) <- levels(atac_obj$celltype)
out_df <- data.frame(peak=rownames(pctcell_mat), pctcell_mat)
write.table(out_df, 'data/pctcells.peak.txt', row.names=F, col.names=T, quote=F, sep='\t')

##=============================


