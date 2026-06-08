#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(scrubletR)
library(reticulate)

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


##
obj_list <- list()
for (sample in names(sample_list)){
        print(sample)

        input_path <- sample_list[[sample]]

        #
        count_file <- Read10X_h5(paste0(input_path, '/filtered_feature_bc_matrix.h5'))[['Gene Expression']]
        tmp_obj <- CreateSeuratObject(counts=count_file, assay = 'RNA')
        tmp_obj$orig.ident <- NULL

	#
	prediction_data <- read.csv(paste0('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1b_collect_data_marrow/data/celltypist_labels.', sample, '.csv'), row.names=1)
	colnames(prediction_data) <- 'celltypist_label'
	tmp_obj <- AddMetaData(tmp_obj, metadata=prediction_data)

        # doublet
        doublet_data <- scrublet_R(seurat_obj=tmp_obj, python_home='/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/python')

        doublet_data$predicted_doublets[doublet_data$predicted_doublets=='TRUE'] <- 'doublet'
        doublet_data$predicted_doublets[doublet_data$predicted_doublets=='FALSE'] <- 'singlet'
        doublet_thresh <- min(doublet_data$doublet_scores[which(doublet_data$predicted_doublets=='doublet')])
        tmp_obj$doublet_score <- doublet_data$doublet_scores

        print(doublet_thresh)
        print(table(doublet_data$predicted_doublets))

        # pct mito/ribo genes
        pct_mito <- PercentageFeatureSet(tmp_obj, pattern = '^MT-')

        a <- PercentageFeatureSet(tmp_obj, pattern = '^RPS')
        b <- PercentageFeatureSet(tmp_obj, pattern = '^RPL')
        pct_ribo <- a+b

        tmp_obj$pct_mito <- pct_mito
        tmp_obj$pct_ribo <- pct_ribo

        #
        tmp_obj$sample <- sample

        #
        print(length(colnames(tmp_obj)))
        obj_list[[sample]] <- tmp_obj
}

rna_obj <- merge(x=obj_list[1][[1]], y=obj_list[2:length(names(obj_list))][1:length(names(obj_list))-1], add.cell.ids=names(obj_list), merge.data=F)
rna_obj[['RNA']] <- JoinLayers(rna_obj[['RNA']])
rna_obj$sample <- factor(rna_obj$sample, levels=names(sample_list))

print(length(colnames(rna_obj)))
print(table(rna_obj$sample))

## filter
rna_obj <- subset(rna_obj, subset= nFeature_RNA > 200 & nFeature_RNA < 5000 & pct_mito < 20 & doublet_score < 0.4)
print(length(colnames(rna_obj)))
print(table(rna_obj$sample))

##======================================
## save

dir.create('seurat_objects')

saveRDS(rna_obj, 'seurat_objects/marrow_rna.rds')

##======================================


