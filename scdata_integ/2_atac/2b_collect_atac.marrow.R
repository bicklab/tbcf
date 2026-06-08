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
        # redeem
        young1t1bmmc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t1bmmc/outs',
	young1t1hspc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t1hspc/outs',
        young1t1hsc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t1hsc/outs',
        young1t2bmmc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t2bmmc/outs',
        young1t2hspc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t2hspc/outs',
        young1t2hsc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young1t2hsc/outs',
        young2bmmc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young2bmmc/outs',
        young2hspc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young2hspc/outs',
        young2hsc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/young2hsc/outs',
        aged1bmmc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/aged1bmmc/outs',
        aged1hspc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/aged1hspc/outs',
        aged2bmmc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/aged2bmmc/outs',
        aged2hspc='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_redeem/1_collect_data/cellranger/aged2hspc/outs',

       # cellarity
        Site1_Donor1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site1_Donor1/outs',
        Site1_Donor2='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site1_Donor2/outs',
        Site1_Donor3='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site1_Donor3/outs',
        Site2_Donor1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site2_Donor1/outs',
        Site2_Donor4='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site2_Donor4/outs',
        Site2_Donor5='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site2_Donor5/outs',
        Site3_Donor3='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site3_Donor3/outs',
        Site3_Donor6='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site3_Donor6/outs',
        Site3_Donor7='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site3_Donor7/outs',
        Site3_Donor10='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site3_Donor10/outs',
        Site4_Donor1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site4_Donor1/outs',
        Site4_Donor8='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site4_Donor8/outs',
        Site4_Donor9='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scmul_cellarity/1_collect_data/cellranger/Site4_Donor9/outs'
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
	peak_data <- Read10X_h5(paste0(input_path, '/filtered_feature_bc_matrix.h5'))[['Peaks']]
        cell_metadata <- read.csv(paste0(input_path, '/per_barcode_metrics.csv'), header=T, row.names=1)
        cell_metadata <- cell_metadata[colnames(peak_data), ]

        #
        frag_file <- paste0(input_path, '/atac_fragments.tsv.gz')
        frag_obj <- CreateFragmentObject(path=frag_file, cells=rownames(cell_metadata), validate.fragments=T, verbose=T)
        feat_mat <- FeatureMatrix(fragments=frag_obj, features=peak_list, cells=rownames(cell_metadata))

        chrom_assay <- CreateChromatinAssay(counts=feat_mat, fragments=frag_obj, sep=c(":", "-"))
        tmp_obj <- CreateSeuratObject(counts=chrom_assay, assay='peaks', meta.data=cell_metadata)
        tmp_obj$sample <- sample

	#
	amulet_path <- '/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1b_collect_data_marrow/amulet_outs'
	amulet_data <- read.table(paste0(amulet_path, '/', sample, '.txt'), header=T)
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
atac_obj$pct_reads_in_peaks <- atac_obj$peak_region_fragments / atac_obj$atac_fragments * 100
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
saveRDS(atac_obj, 'seurat_objects/marrow_atac.rds')



