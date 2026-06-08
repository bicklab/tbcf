#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(preprocessCore)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(viridis)


rna_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/seurat_objects/merged_rna.tcells.rds')
atac_obj <- readRDS('seurat_objects/merged_atac.tcells.init.rds')

marker_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/2_merge/data/marker_genes.tcell_subsets.txt', header=T)
marker_genes <- unique(marker_data$gene)

#
features <- intersect(intersect(marker_genes, rownames(rna_obj)), rownames(atac_obj))
print(length(features))

##================================
## get counts

## atac
atac_mat <- NULL
for (cluster in 1:length(unique(atac_obj$kmeans_cluster))){
	print(cluster)

	tmp_obj <- subset(atac_obj, subset= kmeans_cluster==cluster)
	sumvec <- rowSums(tmp_obj@assays$activity@counts)

	atac_mat <- cbind(atac_mat, sumvec)
}

colnames(atac_mat) <- paste0('cluster_', 1:length(unique(atac_obj$kmeans_cluster)))
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

colnames(rna_mat) <- levels(rna_obj$celltype)
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

query_mat <- scaled_atac
ref_mat <- scaled_rna

#
df <- NULL
for (cluster in colnames(query_mat)){
        cor_vec <- c()
        for (celltype in colnames(ref_mat)){
                cor_val <- cor(ref_mat[,celltype], query_mat[,cluster], method='spearman')
                cor_vec <- c(cor_vec, cor_val)
        }
        names(cor_vec) <- colnames(ref_mat)
        df <- rbind(df, cor_vec)
}
rownames(df) <- colnames(query_mat)

idx_vec <- apply(df, 1, which.max)
celltype_vec <- colnames(ref_mat)[idx_vec]
names(celltype_vec) <- rownames(df)

#
atac_obj$celltype <- as.vector(celltype_vec[paste0('cluster_', atac_obj$kmeans_cluster)])
atac_obj$celltype <- factor(atac_obj$celltype, levels=levels(rna_obj$celltype))
print(table(atac_obj$celltype))

## save
saveRDS(atac_obj, 'seurat_objects/merged_atac.tcells.init2.rds')

##================================
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

celltype_colors <- colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(rna_obj$celltype)))

##
png('plots/tcell_atac.celltype.init.wnn_umap.png', width=6, height=6, res=1000, units='in')
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, label=T, raster=F) + my_theme + ggtitle('')
dev.off()

pdf('plots/tcell_atac.celltype.init.wnn_umap.pdf', width=6, height=7, pointsize=3)
DimPlot(atac_obj, reduction='wnn_umap', group.by='celltype', cols=celltype_colors, label=T, raster=F) + my_theme + ggtitle('')
dev.off()

## coef
for (celltype in levels(atac_obj$celltype)){
        print(celltype)

	cor_vec <- df[, celltype]
	atac_obj$scc <- as.vector(cor_vec[paste0('cluster_', atac_obj$kmeans_cluster)])

        p1 <- FeaturePlot(atac_obj, features='scc', reduction='wnn_umap', min.cutoff='q0', max.cutoff='q100', order=T, raster=F) + scale_color_viridis(option='rocket', direction=-1) + my_theme + ggtitle('')

        png(paste0('plots/tcell_atac.scc.', celltype, '.png'), width=6, height=6, res=1000, units='in')
        plot(p1)
        dev.off()
}

##================================



