#!/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/python


import os
import numpy


DIR=os.getcwd()
sample_dict={
    'f8thy_tot':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/f8thy_tot/outs',
    'p100thy':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p100thy/outs',
    'p96thy1':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p96thy1/outs',
    'p96thy2':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p96thy2/outs',
    'p97thy':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p97thy/outs',
    'p98thy':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p98thy/outs',
    'p99thy1':'/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/scatac_thymus/1a_collect_data_thymus/cellranger/p99thy1/outs'
}

os.system('mkdir '+DIR+'/peaks/')


def run_macs():
    macs_path='/lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env/bin/macs3'

    ncpu=1
    hr=120
    mem='128G'

    min_qval='0.05'
    for sample in sample_dict:
        print(sample)

        input_path=sample_dict[sample]

        output1=open(DIR+'/peakcall.'+sample+'.slurm','w')
        output1.write('#!/bin/bash\n')
        output1.write('\n')
        output1.write('#SBATCH --mem='+mem+'\n')
        output1.write('#SBATCH --cpus-per-task='+str(ncpu)+'\n')
        output1.write('#SBATCH --time='+str(hr)+':00:00\n')
        output1.write('#SBATCH --output=logs/peakcall.'+sample+'.out\n')
        output1.write('#SBATCH --error=logs/peakcall.'+sample+'.err\n')
        output1.write('\n')
        output1.write('source /lab-share/Hem-Sankaran-e2/Public/ajlee/anaconda3/etc/profile.d/conda.sh\n')
        output1.write('conda activate /lab-share/Hem-Sankaran-e2/Public/ajlee/conda_envs/sc_env\n')
        output1.write('\n')

        ##
        output1.write(macs_path+' callpeak -t '+input_path+'/possorted_bam.bam -f BAM -g hs --nomodel --extsize 200 --shift -100 -q '+min_qval+' --outdir '+DIR+'/peaks -n '+sample+'\n')

        output1.write('\n')
        output1.close()

        os.system('sbatch '+DIR+'/peakcall.'+sample+'.slurm')


def macs_peaks():

    min_qval=0.01
    for sample in sample_dict:
        input1=open(DIR+'/peaks/'+sample+'_peaks.narrowPeak','r')
        output1=open(DIR+'/peaks/'+sample+'.sig_peaks.bed','w')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            chrom=each[0]
            pt1=each[1]
            pt2=each[2]
            logq=each[8]

            if float(logq) > -numpy.log10(min_qval):
                new_line=[chrom, pt1, pt2, logq]
                output1.write('\t'.join(new_line)+'\n')
        input1.close()
        output1.close()

    ##
    run_command=['cat ']
    for sample in sample_dict:
        run_command.append(DIR+'/peaks/'+sample+'.sig_peaks.bed ')

    run_command.append('| cut -f1,2,3 | grep \'chr\' | bedtools sort -i stdin | bedtools merge -i stdin -d 200 | sort -V -k1,1 -k2,2n - > '+DIR+'/peaks/macs_peaks.bed')
    run_command=''.join(run_command)
    print(run_command)
    os.system(run_command)

    ##
    run_command=['cat ']
    for sample in sample_dict:
        run_command.append(DIR+'/peaks/'+sample+'.sig_peaks.bed ')

    run_command.append('| cut -f1,2,3 |  grep \'chr\' | coverageBed -counts -a peaks/macs_peaks.bed -b stdin > '+DIR+'/peaks/macs_peaks.coverage.txt')
    run_command=''.join(run_command)
    print(run_command)
    os.system(run_command)


def reproducible_peaks():

    min_sample=3
    input1=open(DIR+'/peaks/macs_peaks.coverage.txt','r')
    output1=open(DIR+'/peaks/reproducible_peaks.thymus.bed','w')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')

        chrom=each[0]
        pt1=each[1]
        pt2=each[2]
        coverage=int(each[3])
        if coverage >= min_sample:
            new_line=[chrom, pt1, pt2]
            output1.write('\t'.join(new_line)+'\n')

    input1.close()
    output1.close()


def merge():

    marrow_peaks='/lab-share/Hem-Sankaran-e2/Public/ajlee/projects/bcf/data/singlecell/integ_bm/1_merge/peaks/reproducible_peaks.sorted.bed'

    #
    run_command='cat '+DIR+'/peaks/reproducible_peaks.thymus.bed '+marrow_peaks+' | cut -f1,2,3 | grep \'chr\' | bedtools sort -i stdin | bedtools merge -i stdin -d 200 | sort -V -k1,1 -k2,2n - > '+DIR+'/peaks/reproducible_peaks.merged.sorted.bed'
    print(run_command)
    os.system(run_command)


#run_macs()
#macs_peaks()
#reproducible_peaks()
merge()



