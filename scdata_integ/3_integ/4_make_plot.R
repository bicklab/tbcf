#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(preprocessCore)
library(stringr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)


set.seed(123)

rna_obj <- readRDS('seurat_objects/integ_rna.lympho.rds')
atac_obj <- readRDS('seurat_objects/integ_atac.lympho.rds')

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

##================================
## get counts

## atac
atac_mat <- NULL
for (celltype in levels(atac_obj$celltype)){
        print(celltype)

        tmp_cells <- colnames(atac_obj)[which(atac_obj$celltype==celltype)]
        tmp_obj <- subset(atac_obj, cells= tmp_cells)
        sumvec <- rowSums(tmp_obj@assays$activity@counts)

        atac_mat <- cbind(atac_mat, sumvec)
}

colnames(atac_mat) <- paste0('atac_', levels(atac_obj$celltype))
atac_mat <- atac_mat[features, ]
print(dim(atac_mat))

# qq
colsum_vec <- colSums(atac_mat)
print(colsum_vec)

max_idx <- which.max(colsum_vec)
target_vec <- atac_mat[, max_idx]
atac_qqmat <- normalize.quantiles.use.target(x=as.matrix(atac_mat), target=target_vec, copy=F)
rownames(atac_qqmat) <- rownames(atac_mat); colnames(atac_qqmat) <- colnames(atac_mat)

log2norm <- function(x){log2(x+1)}
log2atac <- t(apply(atac_qqmat, 1, log2norm))

scaled_atac <- t(scale(t(log2atac)))

## ref counts
rna_mat <- NULL
for (celltype in levels(rna_obj$celltype)){
        print(celltype)

        tmp_cells <- colnames(rna_obj)[which(rna_obj$celltype==celltype)]
        tmp_obj <- subset(rna_obj, cells=tmp_cells)
        sumvec <- rowSums(tmp_obj@assays$RNA$counts)

        rna_mat <- cbind(rna_mat, sumvec)
}

colnames(rna_mat) <- paste0('rna_', levels(rna_obj$celltype))
rna_mat <- rna_mat[features, ]

# qq
colsum_vec <- colSums(rna_mat)
print(colsum_vec)

max_idx <- which.max(colsum_vec)
target_vec <- rna_mat[, max_idx]
rna_qqmat <- normalize.quantiles.use.target(x=as.matrix(rna_mat), target=target_vec, copy=F)
rownames(rna_qqmat) <- rownames(rna_mat); colnames(rna_qqmat) <- colnames(rna_mat)

log2rna <- t(apply(rna_qqmat, 1, log2norm))

scaled_rna <- t(scale(t(log2rna)))

##================================
## cor

rna_data <- scaled_rna
atac_data <- scaled_atac

cor_mat <- NULL
for (i in colnames(atac_data)){
        cor_vec <- c()
        for (j in colnames(rna_data)){
                cor_val <- cor(rna_data[,j], atac_data[,i], method='spearman')
                cor_vec <- c(cor_vec, cor_val)
        }
        names(cor_vec) <- colnames(rna_data)
        cor_mat <- rbind(cor_mat, cor_vec)
}
rownames(cor_mat) <- colnames(atac_data)

## heatmap
nlength <- 100
heatmap_colors <- colorRampPalette(c(rev(brewer.pal(9, 'Blues')), 'white', brewer.pal(9, 'Reds')))(nlength)

min_val <- -0.8
max_val <- 0.8
cor_mat[cor_mat < min_val] <- min_val
cor_mat[cor_mat > max_val] <- max_val

pdf('plots/rna_atac.cor.heatmap.pdf', height=6, width=6, pointsize=3)
pheatmap(
        cor_mat, color= heatmap_colors,
        cluster_rows=F, cluster_cols=F,
        show_colnames=T, show_rownames=T,
        breaks=seq(min_val, max_val, length.out=nlength),
        border_color = NA
)
dev.off()

##================================


