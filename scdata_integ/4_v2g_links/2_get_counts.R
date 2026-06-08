#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(ggplot2)
library(preprocessCore) #BiocManager::install("preprocessCore", configure.args="--disable-threading")


rna_obj <- readRDS('seurat_objects/integ_rna.lympho.subclusters.rds')
atac_obj <- readRDS('seurat_objects/integ_atac.lympho.subclusters.rds')

link_data <- read.table('data/cluster_links.txt', header=T)

## active peaks
minpct <- 0.02
pctcell_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_lympho/1_integ/data/pctcells.peak.txt', header=T, row.names=1)
pctcell_vec <- apply(pctcell_data, 1, max)

active_peaks <- names(pctcell_vec)[which(pctcell_vec > minpct)]
print(length(active_peaks))

## expressed genes
minpct <- 0.05
pctcell_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_lympho/1_integ/data/pctcells.rna.txt', header=T, row.names=1)
pctcell_vec <- apply(pctcell_data, 1, max)

expressed_genes <- names(pctcell_vec)[which(pctcell_vec > minpct)]
print(length(expressed_genes))

##==============================================
## rna
# counts
rna_mat <- NULL
for (cluster in link_data$rna_cluster){
        print(cluster)

	tmp_cells <- colnames(rna_obj)[which(rna_obj$subcluster==cluster)]
        tmp_obj <- subset(rna_obj, cells=tmp_cells)
        sumvec <- rowSums(tmp_obj@assays$RNA$counts)

        rna_mat <- cbind(rna_mat, sumvec)
}

colnames(rna_mat) <- link_data$rna_cluster
rna_mat <- rna_mat[expressed_genes,]
print(dim(rna_mat))

# qq
colsum_vec <- colSums(rna_mat)
print(colsum_vec)

max_idx <- which.max(colsum_vec)
target_vec <- rna_mat[, max_idx]

rna_qq <- normalize.quantiles.use.target(x=as.matrix(rna_mat), target=target_vec, copy=F)
rownames(rna_qq) <- rownames(rna_mat); colnames(rna_qq) <- colnames(rna_mat)

# log2 norm
log2norm <- function(x){log2(x+1)}
log2mat <- t(apply(rna_qq, 1, log2norm))

out_df <- data.frame(gene=rownames(log2mat), log2mat)
write.table(out_df, gzfile('data/rna.log2norm.txt.gz'), row.names=F, col.names=T, quote=F, sep='\t')


## atac
# counts
atac_mat <- NULL
for (cluster in link_data$atac_cluster){
        print(cluster)

	tmp_cells <- colnames(atac_obj)[which(atac_obj$subcluster==cluster)]
        tmp_obj <- subset(atac_obj, cells=tmp_cells)
        val_mat <- tmp_obj@assays$peaks@counts
        sumvec <- rowSums(val_mat)

        atac_mat <- cbind(atac_mat, sumvec)
}

colnames(atac_mat) <- link_data$atac_cluster
atac_mat <- atac_mat[active_peaks,]
print(dim(atac_mat))

# qq
colsum_vec <- colSums(atac_mat)
print(colsum_vec)

max_idx <- which.max(colsum_vec)
target_vec <- atac_mat[, max_idx]

atac_qq <- normalize.quantiles.use.target(x=as.matrix(atac_mat), target=target_vec, copy=F)
rownames(atac_qq) <- rownames(atac_mat); colnames(atac_qq) <- colnames(atac_mat)

# log2 norm
log2norm <- function(x){log2(x+1)}
log2mat <- t(apply(atac_qq, 1, log2norm))

out_df <- data.frame(peak=rownames(log2mat), log2mat)
write.table(out_df, gzfile('data/atac.log2norm.txt.gz'), row.names=F, col.names=T, quote=F, sep='\t')

##==============================================



