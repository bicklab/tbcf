#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


## active peaks
minpct <- 0.02
pctcell_data <- read.table('/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_lympho/1_integ/data/pctcells.peak.txt', header=T, row.names=1)
pctcell_vec <- apply(pctcell_data, 1, max)

a <- strsplit(names(pctcell_vec), '-')
b  <- matrix(unlist(a), ncol=3, byrow=TRUE)
c <- as.data.frame(b); colnames(c) <- c('chr','start','end')
pctcell_data <- data.frame(c, pctcell=pctcell_vec)
pctcell_data$peak <- paste0(pctcell_data$chr, ':', pctcell_data$start, '-', pctcell_data$end)
names(pctcell_vec) <- pctcell_data$peak

active_peaks <- pctcell_data$peak[which(pctcell_data$pctcell > minpct)]
print(length(active_peaks))

peaks_obj <- pctcell_data[which(pctcell_data$pctcell > minpct), c('chr', 'start', 'end')]
peaks_obj <- makeGRangesFromDataFrame(peaks_obj)

##
gtf <- import('/lab-share/Hem-Sankaran-e2/Public/ajlee/genome/cellranger_arc/human_GRCh38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf')
#gtf <- gtf[gtf$gene_type == 'protein_coding']

gene_map <- unique(data.frame(gene_id = gtf$gene_id, gene_name = gtf$gene_name))
gene_vec <- gene_map$gene_name
names(gene_vec) <- gene_map$gene_id

#
txdb <- makeTxDbFromGRanges(gtf)
peak_anno <- annotatePeak(
        peaks_obj,
        TxDb = txdb,
        tssRegion = c(-2500, 2500)
)

anno_df <- as.data.frame(peak_anno)
anno_df$gene_symbol <- gene_vec[anno_df$geneId]
anno_df$peak <- paste0(anno_df$seqnames, ':', anno_df$start, '-', anno_df$end)

anno_df <- anno_df[, c('peak','gene_symbol','geneId','transcriptId','geneStart','distanceToTSS','annotation')]
anno_df$annotation <- ifelse(grepl('Promoter',anno_df$annotation), 'prom', 'distal')

write.table(anno_df, 'data/peak_anno.txt', col.names=T, row.names=F, quote=F, sep='\t')



