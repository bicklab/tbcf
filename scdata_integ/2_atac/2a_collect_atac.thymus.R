#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(Signac)
library(MASS)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(rtracklayer)
library(stringr)


dir.create('seurat_objects')

##
sample_list <- list(
	f8thy_tot='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/f8thy_tot/outs',
	p100thy='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p100thy/outs',
	p96thy1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p96thy1/outs',
	p96thy2='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p96thy2/outs',
	p97thy='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p97thy/outs',
	p98thy='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p98thy/outs',
	p99thy1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p99thy1/outs'
)

##=====================================
## load data

peak_list <- read.table('peaks/reproducible_peaks.merged.sorted.bed', col.names=c('chr','start','end'))
peak_list <- makeGRangesFromDataFrame(peak_list)

obj_list <- list()
for (sample in names(sample_list)){
        print(sample)

        input_path <- sample_list[[sample]]

        #
        peak_data <- Read10X_h5(paste0(input_path, '/filtered_peak_bc_matrix.h5'))
        cell_metadata <- read.csv(paste0(input_path, '/singlecell.csv'), header=T, row.names=1)
        cell_metadata <- cell_metadata[colnames(peak_data), ]

        #
        frag_file <- paste0(input_path, '/fragments.tsv.gz')
        frag_obj <- CreateFragmentObject(path=frag_file, cells=rownames(cell_metadata), validate.fragments=T, verbose=T)
        feat_mat <- FeatureMatrix(fragments=frag_obj, features=peak_list, cells=rownames(cell_metadata))

        chrom_assay <- CreateChromatinAssay(counts=feat_mat, fragments=frag_obj, sep=c(":", "-"))
        tmp_obj <- CreateSeuratObject(counts=chrom_assay, assay='peaks', meta.data=cell_metadata)
        tmp_obj$sample <- sample

	#
	amulet_path <- '/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/amulet_outs'
	amulet_data <- read.table(paste0(amulet_path, '/', sample, '/MultipletProbabilities.txt'), header=T)
	rownames(amulet_data) <- amulet_data$barcode
	amulet_data$amulet_label <- ifelse(amulet_data$q.value < 0.05, 'doublet', 'singlet')
	amulet_data <- amulet_data[, c('amulet_label','q.value')]; colnames(amulet_data)[2] <- 'amulet_q'
	print(table(amulet_data$amulet_label))

	tmp_obj <- AddMetaData(tmp_obj, metadata=amulet_data)

	#
        obj_list[[sample]] <- tmp_obj
}

atac_obj <- merge(x=obj_list[1][[1]], y=obj_list[2:length(names(obj_list))][1:length(names(obj_list))-1], add.cell.ids=names(obj_list))
atac_obj$sample <- factor(atac_obj$sample, levels=names(sample_list))

## gene annotation
annotations <- import('/lab-share/Hem-Sankaran-e2/Public/ajlee/genome/cellranger_arc/human_GRCh38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf')
genome(annotations) <- 'hg38'
seqlevelsStyle(annotations) <- 'UCSC'
annotations$gene_biotype <- annotations$gene_type
Annotation(atac_obj) <- annotations

#
atac_obj$peak_region_fragments <- colSums(atac_obj@assays$peaks@counts)
atac_obj$pct_reads_in_peaks <- atac_obj$peak_region_fragments / atac_obj$passed_filters * 100
atac_obj <- NucleosomeSignal(object=atac_obj)
atac_obj <- TSSEnrichment(object=atac_obj, fast=TRUE)

print(length(colnames(atac_obj)))
print(table(atac_obj$sample))

## filter
print(summary(atac_obj$peak_region_fragments))
print(table(atac_obj$peak_region_fragments > 1000))
print(table(atac_obj$peak_region_fragments < 100000))

print(summary(atac_obj$pct_reads_in_peaks))
print(table(atac_obj$pct_reads_in_peaks > 20))

print(summary(atac_obj$TSS.enrichment))
print(table(atac_obj$TSS.enrichment > 5))

print(table(atac_obj$sample, atac_obj$amulet_label))

atac_obj <- subset(x=atac_obj, subset = amulet_label == 'singlet' & peak_region_fragments > 1000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 20 & TSS.enrichment > 5 & nucleosome_signal < 10) ##

print(length(colnames(atac_obj)))
print(table(atac_obj$sample))

## gene activity 
gene.activities <- GeneActivity(
        atac_obj, assay='peaks', features=NULL,
        extend.upstream=2000, extend.downstream=0,
        biotypes='protein_coding', max.width=NULL, process_n=1000
)

atac_obj[['activity']] <- CreateAssayObject(counts = gene.activities)
atac_obj <- NormalizeData(object=atac_obj, assay='activity', normalization.method='LogNormalize', scale.factor=median(atac_obj$nCount_activity))

## save 
saveRDS(atac_obj, 'seurat_objects/thymus_atac.rds')



