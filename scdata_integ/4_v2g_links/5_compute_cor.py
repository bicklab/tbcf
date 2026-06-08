#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/python


import os
import sys
import gzip
import numpy
from scipy.stats import pearsonr


DIR=os.getcwd()

rnacount_file=DIR+'/data/rna.log2norm.txt.gz'
peakcount_file=DIR+'/data/atac.log2norm.txt.gz'
win_size='250kb'


def load_counts():

    ## genes
    rna_dict={}
    input1=gzip.open(rnacount_file,'rt')
    for line in input1:
        each=line.strip().split('\t')
        if not each[0].count('gene') > 0:
            gene=each[0]

            vals=numpy.array(each[1:], dtype='float')
            rna_dict[gene]=vals
    input1.close()

    print('rna loaded')
    print('gene #: ', str(len(rna_dict)))

    ## peaks
    atac_dict={}
    input1=gzip.open(peakcount_file,'rt')
    for line in input1:
        each=line.strip().split('\t')
        if not each[0].count('peak') > 0:
            peak=each[0].replace(':', '-')

            vals=numpy.array(each[1:], dtype='float')
            atac_dict[peak]=vals
    input1.close()

    print('atac loaded')
    print('peak #: ', str(len(atac_dict)))

    return rna_dict, atac_dict


def compute_cor():
    rna_dict, atac_dict=load_counts()

    # prom pairs
    output1=open(DIR+'/data/prom_pairs.cor.txt','w')
    output1.write('prom_peak\tprom_gene\ttss\tdist\tpcc\tpval\n')

    input1=open(DIR+'/data/prom_pairs.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        [prom_peak, prom_gene, tss, dist]=line.strip().split('\t')

        rna_vals=rna_dict[prom_gene]
        prom_vals=atac_dict[prom_peak]

        pcc, pval = pearsonr(rna_vals, prom_vals, alternative='greater')

        new_line=[prom_peak, prom_gene, tss, dist, str(pcc), str(pval)]
        output1.write('\t'.join(new_line)+'\n')

    input1.close()
    output1.close()

    # distal pairs
    output1=open(DIR+'/data/distal_pairs.'+win_size+'.cor.txt','w')
    output1.write('prom_peak\tprom_gene\tdistal_peak\tdist\tprom_pcc\tprom_pval\tgene_pcc\tgene_pval\n')

    input1=open(DIR+'/data/distal_pairs.'+win_size+'.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        [prom_peak, prom_gene, distal_peak, dist]=line.strip().split('\t')

        rna_vals=rna_dict[prom_gene]
        prom_vals=atac_dict[prom_peak]
        enh_vals=atac_dict[distal_peak]

        prom_pcc, prom_pval=pearsonr(prom_vals, enh_vals, alternative='greater')
        gene_pcc, gene_pval=pearsonr(rna_vals, enh_vals, alternative='greater')

        new_line=[prom_peak, prom_gene, distal_peak, dist, str(prom_pcc), str(prom_pval), str(gene_pcc), str(gene_pval)]
        output1.write('\t'.join(new_line)+'\n')

    input1.close()
    output1.close()



compute_cor()



