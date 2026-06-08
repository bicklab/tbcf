#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(Seurat)
library(scrubletR)
library(reticulate)

##
sample_list <- list(
f64thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/f64thy_tot1/outs',
f67thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/f67thy_tot1/outs',
t03thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/t03thy_tot1/outs',
t03thy_tot2='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/t03thy_tot2/outs',
t06thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/t06thy_tot1/outs',
t07thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/t07thy_tot1/outs',
c34thy_tot='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/c34thy_tot/outs',
c40thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/c40thy_tot1/outs',
c41thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/c41thy_tot1/outs',
a16thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/a16thy_tot1/outs',
a16thy_tot2='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/a16thy_tot2/outs',
a16thy_tot5='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/a16thy_tot5/outs',
a16thy_tot6='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/a16thy_tot6/outs',
a43thy_tot1='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/cellranger/a43thy_tot1/outs'
)


##
obj_list <- list()
for (sample in names(sample_list)){
        print(sample)

        input_path <- sample_list[[sample]]

        #
        count_file <- Read10X_h5(paste0(input_path, '/filtered_feature_bc_matrix.h5'))
        tmp_obj <- CreateSeuratObject(counts=count_file, assay = 'RNA')
        tmp_obj$orig.ident <- NULL

	#
	prediction_data <- read.csv(paste0('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scrna_thymus/1a_collect_data_thymus/data/celltypist_labels.', sample, '.csv'), row.names=1)
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

saveRDS(rna_obj, 'seurat_objects/thymus_rna.rds')

##======================================


