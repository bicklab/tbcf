#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/R


win_size <- '250kb'
min_fdr <- 0.01


## prom pairs
data1 <- read.table('data/prom_pairs.cor.txt', header=T)
data1$fdr <- p.adjust(data1$pval, method='fdr')

data1$label <- ifelse(data1$fdr < min_fdr, 'sig', 'insig')
print(table(data1$label))

write.table(data1, 'data/prom_pairs.cor.sig.txt', col.names=T, row.names=F, quote=F, sep='\t')


## distal pairs
data1 <- read.table(paste0('data/distal_pairs.', win_size, '.cor.txt'), header=T)

data1$prom_fdr <- p.adjust(data1$prom_pval, method='fdr')
data1$gene_fdr <- p.adjust(data1$gene_pval, method='fdr')

prom_idx <- which(data1$prom_fdr < min_fdr)
gene_idx <- which(data1$gene_fdr < min_fdr)
sig_idx <- intersect(prom_idx, gene_idx)

data1$label <- 'insig'
data1$label[sig_idx] <- 'sig'
print(table(data1$label))

write.table(data1, paste0('data/distal_pairs.', win_size, '.cor.sig.txt'), col.names=T, row.names=F, quote=F, sep='\t')




