#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/python


import os
import gzip
import random


DIR=os.getcwd()
rnacount_file=DIR+'/data/rna.log2norm.txt.gz'
peak_file=DIR+'/data/peak_anno.txt'

win_size=[250000,'250kb']


def load_genes():
    ## expressed genes
    gene_dict={}
    input1=gzip.open(rnacount_file,'rt')
    for line in input1:
        each=line.strip().split('\t')
        if not each[0].count('gene') > 0:
            gene=each[0]
            gene_dict[gene]=1
    input1.close()

    return gene_dict


def prom_links():
    gene_dict=load_genes()

    #
    output1=open(DIR+'/data/prom_pairs.txt','w')
    output1.write('prom_peak\tprom_gene\ttss\tdist\n')

    input1=open(peak_file,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peak=each[0].replace(':', '-')
        gene=each[1]
        tss=each[4]
        dist=str(abs(int(each[5])))
        peak_label=each[6]
        if gene in gene_dict:
            if peak_label=='prom':

                new_line=[peak, gene, tss, dist]
                output1.write('\t'.join(new_line)+'\n')

    input1.close()
    output1.close()


def distal_links():
    gene_dict=load_genes()

    #
    all_dict={}
    prom_dict={}
    input1=open(peak_file,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peak=each[0].replace(':', '-')
        chrom=peak.split('-')[0]
        gene=each[1]
        tss=each[4]
        peak_label=each[6]

        #
        if chrom in all_dict:
            tmp=all_dict[chrom]
            tmp.append([peak, gene])
            all_dict[chrom]=tmp
        else:
            all_dict[chrom]=[[peak, gene]]

        #
        if gene in gene_dict:
            if peak_label=='prom':
                prom_dict[peak]=[gene, int(tss)]

    input1.close()

    #
    output1=open(DIR+'/data/distal_pairs.'+win_size[1]+'.txt','w')
    output1.write('prom_peak\tprom_gene\tdistal_peak\tdist\n')
    for prom_peak in prom_dict:
        prom_gene=prom_dict[prom_peak][0]
        prom_tss=prom_dict[prom_peak][1]
        chrom=prom_peak.split('-')[0]
        for i in all_dict[chrom]:
            distal_peak=i[0]
            distal_gene=i[1]

            if distal_gene!=prom_gene:
                start=int(distal_peak.split('-')[1])
                end=int(distal_peak.split('-')[2])

                dist1= abs(prom_tss-start)
                dist2= abs(prom_tss-end)
                dist= min(dist1, dist2)

                if dist <= win_size[0]:
                    new_line=[prom_peak, prom_gene, distal_peak, str(dist)]
                    output1.write('\t'.join(new_line)+'\n')
    output1.close()


prom_links()
distal_links()



