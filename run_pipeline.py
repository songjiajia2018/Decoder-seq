# -*- coding: utf-8 -*-

import os
import re
import argparse
import multiprocessing as mp
from multiprocessing import Pool
import time
import numpy as np

#######input parameters#################
parser = argparse.ArgumentParser(description='''pileline for yanglab's ST-SEQ''')
parser.add_argument('--map', choices=['hg38', 'mm10', 'hm'], type = str, help = 'choose the reference',default = 'mm10')
parser.add_argument('--thread', default = 3, type = int, help = 'Specify the number of \
                    threads for qc step and dge step. The recommended number is less than 5. ')
parser.add_argument('--R1', type=str, required=True, help='input R1 fastq')
parser.add_argument('--R2', type=str, required=True, help='input R2 fastq')
parser.add_argument('--bcx', type=str, help='input barcode X', default = 'lib_barcode_X.txt')
parser.add_argument('--bcy', type=str, help='input barcode Y',default = 'lib_barcode_Y.txt')
parser.add_argument('--startx', type=int, help='position of barcode X',default = 0)
parser.add_argument('--endx', type=int, help='position of barcode X',default = 8)
parser.add_argument('--starty', type=int, help='position of barcode Y',default = 38)
parser.add_argument('--endy', type=int, help='position of barcode Y',default = 46)
parser.add_argument('--startumi', type=int, help='position of umi',default = 46)
parser.add_argument('--endumi', type=int, help='position of umi',default = 58)



args = parser.parse_args()
threads = args.thread
map_choice = args.map
fastq1 = args.R1
fastq2 = args.R2
barcode_X = args.bcx  
barcode_Y = args.bcy
start_barcode_X = args.startx 
end_barcode_X = args.endx 
start_barcode_Y = args.starty 
end_barcode_Y = args.endy
start_umi = args.startumi
end_umi = args.endumi



human_path = '/home/songjia/reference/human'
mouse_path = '/home/songjia/reference/mouse'
hm_path = '/home/songjia/newdisk/YK/reference/hg19_mm10'
#dge_filter_ref = '/home/songjia/bigdisk/cyw/PIPE/smRNA_gene_filter.txt'


if map_choice == 'hg38':
    path = human_path
    gtf = 'Homo_sapiens.GRCh38.98.gtf'
elif map_choice == 'mm10':
    path = mouse_path
    gtf = 'Mus_musculus.GRCm38.98.gtf'
elif map_choice == 'hm':
    path = hm_path
    gtf = 'hg19_mm10_transgenes_species_gene_name.gtf'




#######0_gunzip if needed###########
suffix = fastq1.split('.')[-1]
suffix_fq = fastq1.split('.')[-2]
prefix_r1 = fastq1.split('.')[0]
prefix_r2 = fastq2.split('.')[0]
begin_time = time.time()
if suffix == 'gz':
    os.system('gunzip -f {0}'.format(fastq1))
    os.system('gunzip -f {0}'.format(fastq2))
    print(' 0: input fastq R1 and R2 has been decompressed!'+ '\n')
    time_point_unzip = time.time()
    time_run_unzip = float(time_point_unzip - begin_time)/60.0
    run_time_unzip_info = ' unzip completed, run time:' + str(time_run_unzip) + ' min' + '\n'
    print(run_time_unzip_info + '\n')
    fastq1 = prefix_r1 + '.' + suffix_fq
    fastq2 = prefix_r2 + '.' + suffix_fq
else:
    fastq1 = prefix_r1 + '.' + suffix
    fastq2 = prefix_r2 + '.' + suffix


#######1_pollution#####################

print(str(fastq1)+',' +str(fastq2))
#fastq2 = str(sys.argv[1])
outfa = fastq2.split('.fq')[0] + '.fa'

if os.path.exists(prefix_r2 + '_blast_final_result.txt'):
    print('Pollution already done, passed pollution' + '\n')
else:
    print('Assessment of pollution for ' + outfa + '\n')
    begin_time_pollution = time.time()
    os.system('''head -n 20000 {0}| sed '/^@/!d;s//>/;N'> {1}'''.format(fastq2,outfa))
    print('1_pollution analysis start!' + '\n')
    os.system('python 1_blast_species.py -f {0}'.format(outfa))
    time_point_pollution = time.time()
    time_run_pollution = float(time_point_pollution - begin_time_pollution)/60.0
    run_time_pollution_info = ' pollution completed, run time:' + str(time_run_pollution) + ' min' + '\n'
    print(run_time_pollution_info + '\n')



#######2_trim_steq_before_star#####################
trim_params=[]
trim_params.append((fastq1, fastq2, barcode_X, barcode_Y,start_barcode_X,end_barcode_X, start_barcode_Y, end_barcode_Y, start_umi, end_umi))
if os.path.exists(prefix_r1 + '_trim_report.txt'):
    print('Trim already done, PASS' + '\n')
else:
    begin_time_qc = time.time()
    print('start barcode filter process,' + ' the parameters are :' +'\n')
    print(trim_params, flush=True)
    os.system('python 2_trim_steq_before_star.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(fastq1, fastq2, barcode_X, barcode_Y,start_barcode_X,end_barcode_X, start_barcode_Y, end_barcode_Y, start_umi, end_umi))
    time_point_qc = time.time()
    time_run_qc = float(time_point_qc - begin_time_qc)/60.0
    run_time_qc_info = ' barcode filter completed, run time:' + str(time_run_qc) + ' min' + '\n'
    print(run_time_qc_info + '\n')


#######3_STAR#####################
prefix_trim_r2 = prefix_r2 + '_trim_'
fastq2_trim= prefix_r2 + '_trim.fq'
if os.path.exists(prefix_r2 + '_trim_Log.final.out'):
    print('STAR already done, PASS' + '\n')
else:
    begin_time_map = time.time()
    print('start STAR process' + '\n')
    os.system('''STAR --runThreadN 20 --genomeDir {PATH}/STAR --outFileNamePrefix {prefix} --sjdbGTFfile {PATH}/{GTF} --outSAMunmapped Within --readFilesIn {fq2};'''.format(PATH=path, prefix=prefix_trim_r2, GTF=gtf, fq2=fastq2_trim))
    time_point_map = time.time()
    time_run_map = float(time_point_map - begin_time_map)/60.0
    run_time_map_info = 'STAR mapping run time:' + str(time_run_map) + ' min' + '\n'
    print(run_time_map_info + '\n')



#######4_extract barcode / umi+barcode#####################
begin_time_cluster_tag = time.time()
print('start cluster and tag '+ '\n')
prefix_trim_r1 = prefix_r1 + '_trim_'
fastq1_trim= prefix_r1 + '_trim.fq'
out_bc = prefix_trim_r1 +'barcode.fasta'
with open (fastq1_trim,'r') as ref:
    rowcount = 0
    with open(out_bc,'w') as extract:
        for line in ref:
            rowcount +=1
            if rowcount%4 ==1:
                line = line.rstrip()
                id_fasta ='>'+line+'\n'
                extract.write(id_fasta)
            elif rowcount%4 ==2:
                line = line.rstrip()
                barcodeX = line[start_barcode_X:end_barcode_X]
                barcodeY = line[start_barcode_Y:end_barcode_Y]
                out_str = barcodeX + barcodeY+'\n'
                extract.write(out_str)
            else:
                continue
print('extract barcode done')

out_umibc = prefix_trim_r1 + 'umiUBC.fasta'
with open (fastq1_trim,'r') as ref:
    rowcount = 0
    with open(out_umibc,'w') as extract:
        for line in ref:
            rowcount +=1
            if rowcount%4 ==1:
                line = line.rstrip()
                id_fasta ='>'+line+'\n'
                extract.write(id_fasta)
            elif rowcount%4 ==2:
                line = line.rstrip()
                umi = line[start_umi:end_umi]
                barcodeX = line[start_barcode_X:end_barcode_X]
                barcodeY = line[start_barcode_Y:end_barcode_Y]
                out_str = umi+barcodeX+barcodeY+'\n'
                extract.write(out_str)
            else:
                continue
print('extract umi-ubc done')

#######5_cluster barcode / umi+barcode#####################
umi_len = end_umi-start_umi
os.system('''starcode/starcode -t 10 -d 0 --seq-id -i {0}barcode.fasta -o {0}cluster_cell.txt;
	cd starcode
             python2 starcode-umi --umi-len {1} --umi-d 1 --seq-d 0 --seq-threads 10 \
             --seq-id  --seq-trim 0 ../{0}umiUBC.fasta > {0}cluster_umiubc.txt;
             mv {0}cluster_umiubc.txt ../;
             '''.format(prefix_trim_r1,umi_len))
print('starcode done!!')



#######6_tag#####################

os.system('''samtools view -b -@ 10 -o {prefix}Aligned.out.bam {prefix}Aligned.out.sam;
             TagReadWithGeneExonFunction I={prefix}Aligned.out.bam O={prefix}med.bam \
             ANNOTATIONS_FILE={PATH}/{GTF} TAG=GE;
             TagReadWithGeneFunction I={prefix}med.bam O={prefix}tagged.bam ANNOTATIONS_FILE={PATH}/{GTF};
             samtools view -h -@ 10 -o {prefix}tagged.sam {prefix}tagged.bam;
             rm {prefix}med.bam {prefix}tagged.bam 
          '''.format(PATH=path, prefix=prefix_trim_r2, GTF=gtf)
                  )
print(prefix_trim_r2 + ' done!')

print(' all done')
time_point_cluster_tag = time.time()
time_run_cls_tag = float(time_point_cluster_tag - begin_time_cluster_tag)/60.0
run_time_clust_tag_info = 'clusting and tagging completed, run time:' + str(time_run_cls_tag) + ' min' + '\n'
print(run_time_clust_tag_info + '\n')


#######7_build_expresion_matrix#####################
#sam = prefix_r2 + '_Aligned.out_tagged.sam'
umi = prefix_trim_r1 + 'cluster_umiubc.txt'
cb = prefix_trim_r1 + 'cluster_cell.txt'
tagsam = prefix_trim_r2 + 'tagged.sam'
os.system('python 3_aftercode.py {0} {1} {2} {3}'.format(fastq1_trim, umi, cb, tagsam))
end_time = time.time()
time_run = float(end_time-time_point_cluster_tag)/60.0
run_time_DGE_info = 'generating DGE completed, run time:' + str(time_run) + 'min' + '\n'
print(run_time_DGE_info + '\n')

run_time = float(end_time-begin_time)/60.0
run_time_info = 'D_pipe run time:' + str(run_time) + 'min' + '\n'

print(run_time_info + '\n')

#log_name = 'D_pipe_log'
#with open(log_name, 'w')as OUT:
#    str_write = run_time_qc_info + run_time_map_info + run_time_clust_tag_info + run_time_DGE_info + run_time_info
#    print(str_write + '\n')
#    OUT.write(str_write)

##########8.st_pipe,statistic of expression matrix#############

os.system('''cat {R1}_trim_report.txt;
	     cat {R2}_blast_final_result.txt;
	     cat {R2}_trim_Log.final.out;
	     wc -l {R1}_trim_cluster_umiubc.txt;
             conda activate r4;
             Rscript st_pipe1.R {R2}_trim_tagged_dge 1 50 1 50;
	     Rscript st_xiuqiu_3marker.R {R2}_trim_tagged_dge 2000;
          '''.format(R1=prefix_r1, R2=prefix_r2))
