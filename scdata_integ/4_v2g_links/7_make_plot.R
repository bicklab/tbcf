#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(preprocessCore)
library(stringr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrastr)
library(scales)
library(pheatmap)
library(RColorBrewer)
library(viridis)


set.seed(123)

dir.create('plots')


rna_obj <- readRDS('seurat_objects/integ_rna.lympho.subclusters.rds')
atac_obj <- readRDS('seurat_objects/integ_atac.lympho.subclusters.rds')

# 
link_data <- read.table('data/cluster_links.txt', header=T)

link_data$subcluster <- link_data$rna_cluster
link_data$rna_cluster <- paste0('rna_', link_data$rna_cluster)
link_data$atac_cluster <- paste0('atac_', link_data$atac_cluster)

rna_cluster_vec <- link_data$celltype
names(rna_cluster_vec) <- link_data$rna_cluster
atac_cluster_vec <- link_data$celltype
names(atac_cluster_vec) <- link_data$atac_cluster

##=====================================
##
x_shift <- 1
rna_coords <- data.frame(
	x=rescale(rna_obj@reductions$css_umap@cell.embeddings[,2], to=c(0,1)), 
	y=rescale(rna_obj@reductions$css_umap@cell.embeddings[,1], to=c(0,1)), 
	subcluster=paste0('rna_', rna_obj$subcluster), celltype=rna_obj$celltype
)
atac_coords <- data.frame(
	x=rescale(atac_obj@reductions$wnn_umap@cell.embeddings[,2], to=c(0,1)) + x_shift, 
	y=rescale(atac_obj@reductions$wnn_umap@cell.embeddings[,1], to=c(0,1)), 
	subcluster=paste0('atac_', atac_obj$subcluster), celltype=atac_obj$celltype
)
nsample <- 20000
merged_coords <- rbind(rna_coords[sample(1:nrow(rna_coords), size=nsample, replace=F),], atac_coords[sample(1:nrow(atac_coords), size=nsample, replace=F),])
merged_coords$subcluster <- factor(merged_coords$subcluster, levels=c(link_data$rna_cluster, link_data$atac_cluster))

##
rna_means <- data.frame(
	x=tapply(rescale(rna_obj@reductions$css_umap@cell.embeddings[,2], to=c(0,1)), rna_obj$subcluster, mean),
	y=tapply(rescale(rna_obj@reductions$css_umap@cell.embeddings[,1], to=c(0,1)), rna_obj$subcluster, mean),
	subcluster=paste0('rna_', levels(rna_obj$subcluster))
)
atac_means <- data.frame(
	x=tapply(rescale(atac_obj@reductions$wnn_umap@cell.embeddings[,2], to=c(0,1)), atac_obj$subcluster, mean) + x_shift,
	y=tapply(rescale(atac_obj@reductions$wnn_umap@cell.embeddings[,1], to=c(0,1)), atac_obj$subcluster, mean),
	subcluster=paste0('atac_', levels(atac_obj$subcluster))
)
mean_coords <- rbind(rna_means, atac_means)
rownames(mean_coords) <- 1:nrow(mean_coords)
mean_coords$subcluster <- factor(mean_coords$subcluster, levels=c(link_data$rna_cluster, link_data$atac_cluster))

#
graph <- igraph::graph_from_data_frame(link_data)
graph <- tidygraph::as_tbl_graph(graph, directed=F)
graph <- tidygraph::activate(graph, nodes)
E(graph)$edge_type <- factor(E(graph)$subcluster, levels=link_data$subcluster) 

layout <- ggraph::create_layout(graph, mean_coords)

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

colors <- rep(hue_pal()(nrow(link_data)), 2)
p1 <- ggraph(layout) +
	geom_node_point(aes(fill=subcluster, color=subcluster), size=1.5, shape=21) +
	geom_edge_link(aes(color=edge_type), alpha=0.1, edge_width=0.4) +
	ggrastr::rasterise(geom_point(inherit.aes=FALSE, data=merged_coords, aes(x=x, y=y, color=subcluster), alpha=0.05, size=0.2), dpi=500) +

	scale_colour_manual(values=colors) +
	scale_fill_manual(values=colors) +
	scale_edge_colour_manual(values=colors, 'lightgrey') +
	my_theme + NoLegend()

pdf('plots/merged_umap.subclusters.nolegend.pdf', height=5, width=10, pointsize=3)
plot(p1)
dev.off()

##=====================================

