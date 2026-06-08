#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(preprocessCore)
library(ggplot2)


dir.create('seurat_objects')
dir.create('data')


rna_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_lympho/1_integ/seurat_objects/integ_rna.lympho.rds')
atac_obj <- readRDS('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_lympho/1_integ/seurat_objects/integ_atac.lympho.rds')

print(table(rna_obj$celltype))
print(table(atac_obj$celltype))

##====================================
## kmeans clusters

#
k_list <- list(
'HSC'=5,
'MPP_MyLy'=5,
'MPP_MkEry'=5,
'LMPP'=5,
'MLP'=4,
'CLP'=2,
'PreProB'=3,
'Cycling_ProB'=3,
'VDJ_ProB'=5,
'Large_PreB'=4,
'Small_PreB'=3,
'ImmatureB'=5,
'MatureB'=5,
'DN'=5,
'DP'=5,
'abT(entry)'=5,
'CD8aa'=5,
'CD8+T'=5,
'CD8+Tmem'=4,
'CD4+T'=5,
'CD4+Tmem'=4,
'Treg'=5,
'T(agonist)'=3,
'rdT'=2,
'NK'=5,       
'ILC3'=5
)

##
min_cells <- 100
rna_clusters <- data.frame()
atac_clusters <- data.frame()
cluster_vec <- c()
for (celltype in names(k_list)){
	print(celltype)

	k <- k_list[[celltype]]

	#
	tmp_cells <- colnames(rna_obj)[which(rna_obj$celltype==celltype)]
	sub_obj <- subset(rna_obj, cells=tmp_cells)
	tmp1_data <- data.frame(kmc=kmeans(sub_obj@reductions$css_umap@cell.embeddings, centers=k)$cluster, celltype)
	tmp1_table <- table(tmp1_data$kmc)[order(table(tmp1_data$kmc), decreasing=T)]
	tmp1_valids <- names(tmp1_table)[tmp1_table > min_cells]

	#
	tmp_cells <- colnames(atac_obj)[which(atac_obj$celltype==celltype)]
	sub_obj <- subset(atac_obj, cells=tmp_cells)
	tmp2_data <- data.frame(kmc=kmeans(sub_obj@reductions$wnn_umap@cell.embeddings, centers=k)$cluster, celltype)
	tmp2_table <- table(tmp2_data$kmc)[order(table(tmp2_data$kmc), decreasing=T)]
	tmp2_valids <- names(tmp2_table)[tmp2_table > min_cells]

	#
	ncluster <- min(length(tmp1_valids), length(tmp2_valids)) 

	tmp1_data <- tmp1_data[as.character(tmp1_data$kmc) %in% tmp1_valids[1:ncluster],]
	tmp2_data <- tmp2_data[as.character(tmp2_data$kmc) %in% tmp2_valids[1:ncluster],]

	tmp1_vec <- as.character(1:ncluster); names(tmp1_vec) <- names(table(tmp1_data$kmc))[order(table(tmp1_data$kmc), decreasing=T)] 
	tmp2_vec <- as.character(1:ncluster); names(tmp2_vec) <- names(table(tmp2_data$kmc))[order(table(tmp2_data$kmc), decreasing=T)]

	tmp1_data$kmc_adj <- tmp1_vec[as.character(tmp1_data$kmc)]
	tmp2_data$kmc_adj <- tmp2_vec[as.character(tmp2_data$kmc)]

	#print(table(tmp1_data$kmc_adj))
	#print(table(tmp2_data$kmc_adj))

	#
	rna_clusters <- rbind(rna_clusters, tmp1_data)
	atac_clusters <- rbind(atac_clusters, tmp2_data)

	cluster_vec <- c(cluster_vec, paste0(celltype, '_', 1:ncluster))
}

rna_clusters$subcluster <- paste0(rna_clusters$celltype, '_', rna_clusters$kmc_adj)
atac_clusters$subcluster <- paste0(atac_clusters$celltype, '_', atac_clusters$kmc_adj)

rna_clusters$subcluster <- factor(rna_clusters$subcluster, levels=cluster_vec)
atac_clusters$subcluster <- factor(atac_clusters$subcluster, levels=cluster_vec)

print(length(cluster_vec))
print(table(unique(rna_clusters)$celltype))
print(table(rna_clusters$celltype))
print(table(rna_clusters$subcluster))
print(table(rna_clusters$subcluster))

#
rna_obj <- subset(rna_obj, cells=rownames(rna_clusters))
atac_obj <- subset(atac_obj, cells=rownames(atac_clusters))

rna_vec <- rna_clusters$subcluster; names(rna_vec) <- rownames(rna_clusters)
atac_vec <- atac_clusters$subcluster; names(atac_vec) <- rownames(atac_clusters)

rna_obj$subcluster <- as.character(rna_vec[colnames(rna_obj)])
atac_obj$subcluster <- as.character(atac_vec[colnames(atac_obj)])

rna_obj$subcluster <- factor(rna_obj$subcluster, levels=cluster_vec)
atac_obj$subcluster <- factor(atac_obj$subcluster, levels=cluster_vec)


##====================================
## get counts

## marker genes
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

## atac
atac_mat <- NULL
for (cluster in cluster_vec){
        #print(cluster)

	tmp_cells <- colnames(atac_obj)[which(atac_obj$subcluster==cluster)]
        tmp_obj <- subset(atac_obj, cells=tmp_cells)
        sumvec <- rowSums(tmp_obj@assays$activity@counts)

        atac_mat <- cbind(atac_mat, sumvec)
}

colnames(atac_mat) <- cluster_vec
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

## rna
rna_mat <- NULL
for (cluster in cluster_vec){
        #print(cluster)

        tmp_cells <- colnames(rna_obj)[which(rna_obj$subcluster==cluster)]
        tmp_obj <- subset(rna_obj, cells=tmp_cells)
        sumvec <- rowSums(tmp_obj@assays$RNA$counts)

        rna_mat <- cbind(rna_mat, sumvec)
}

colnames(rna_mat) <- cluster_vec
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

##====================================
## link by Hungarian algorithm - globally optimal 1:1 pairing 

library(clue)

link_data <- data.frame()
for (celltype in levels(rna_obj$celltype)){
	print(celltype)

	celltype_clusters <- unique(rna_obj$subcluster[which(rna_obj$celltype==celltype)])
	tmp_mat <- NULL
	for (rna_cluster in celltype_clusters){
		tmp_vec <- c()
		for (atac_cluster in celltype_clusters){

			cor_val <- cor(scaled_rna[,rna_cluster], scaled_atac[,atac_cluster], method='spearman')
			tmp_vec <- c(tmp_vec, cor_val)
		}
		tmp_mat <- rbind(tmp_mat, tmp_vec)
	}

	#
	colnames(tmp_mat) <- celltype_clusters
	rownames(tmp_mat) <- celltype_clusters

	min_val <- min(tmp_mat, na.rm = TRUE)
	tmp_mat <- tmp_mat - min_val

	#
	assignment <- solve_LSAP(tmp_mat, maximum=T)
	tmp_data  <- data.frame(rna_cluster = rownames(tmp_mat), atac_cluster = colnames(tmp_mat)[assignment], celltype)
	link_data <- rbind(link_data, tmp_data)
}

## save
write.table(link_data, 'data/cluster_links.txt', col.names=T, row.names=F, quote=F, sep='\t')

saveRDS(rna_obj, 'seurat_objects/integ_rna.lympho.subclusters.rds')
saveRDS(atac_obj, 'seurat_objects/integ_atac.lympho.subclusters.rds')

##====================================




