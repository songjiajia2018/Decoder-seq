# -*- coding: utf-8 -*-

import sys
from collections import defaultdict ,Counter
import time

fastq1 = str(sys.argv[1]) #file used to generate seq id 
cluster_fq1_ref = str(sys.argv[2]) #file with umi-ubc cluster info
cluster_cell = str(sys.argv[3])# file with cluster info (only clustered by barcode), the file represent different cells
tagged_sam = str(sys.argv[4]) #SAM file come from STAR,tagged by TAGWITHGENEFUNCTION and TAGWITHGENEEXONFUNCTION
out_prefix = tagged_sam.split('.sam')[0]
cell_barcode_save = []# save barcode construction
DGE = out_prefix + '_dge'
report = out_prefix +'_endreport'
begin_time = time.time()
cluster_ref_dict = {}#we use this to record (reads number in fastq1: clustered umi)
fq1_ref_dict = {}#we use this to record (id:clutered umi)
#fq2_ref_dict = {}
umi_genenumber_dict = defaultdict(Counter) #use this to record the frequency each gene matched by umis
cell_umi_dict = defaultdict(dict)#use this to record cell-umi information
dge_dict = defaultdict(Counter)#data used to output digital gene expression matrix
ID_error = 0
barcode_umi_match_error = 0
umi_discarded = 0
# save umi cluster information:a group of reads share a same umi, record this in format:reads number-umi-number
with open(cluster_fq1_ref,'r') as REF:
    umi_group_count = 0
    for line in REF:
        umi_group_count += 1
        line = line.rstrip()
        fq1_seq_numstr = line.split()[2]
        fq1_samegroup_list = fq1_seq_numstr.split(',')
        for fq1_seqnumber in fq1_samegroup_list:
            cluster_ref_dict[fq1_seqnumber] = umi_group_count
#match cell and umi, because we don't want to avovid umi repeat , so we use 2-d dict in place of dict(list)
with open(cluster_cell,'r') as CELL:
    cell_group_count = 0
    for line in CELL:
        cell_group_count += 1
        line = line.rstrip()
        cell_barcode = line.split()[0]
        #using the sequence to match the barcode content and the barcode number
        cell_barcode_save.append(cell_barcode)
        fq1_cell_number = line.split()[2]
        fq1_cell_list = fq1_cell_number.split(',')
        for cell_number in fq1_cell_list:
            flag = cluster_ref_dict.get(cell_number,-1)
            if flag ==-1:
                barcode_umi_match_error +=1
            else:
                cell_umi_dict[cell_group_count][flag]=''
cell_umi_dict = dict(cell_umi_dict)
#trans cluster imput number to seq ID
with open(fastq1,'r') as ID_TRANS:
    rowcount = 0
    readscount = 0
    for line in ID_TRANS:
        rowcount +=1
        if rowcount%4 ==1:
            readscount += 1
            readscount_str = str(readscount)
            flag = cluster_ref_dict.get(readscount_str,-1)
            if flag == -1:
                continue
            else:
                line = line.rstrip()
                id_fq1 = line.split()[0].split('@')[1]
                fq1_ref_dict[id_fq1] = flag

#the most difficult step,the structure of tagged line in sam files are different.
with open(tagged_sam,'r') as GE_MATCH:
    for line in GE_MATCH:
        if line.startswith('@'):
            continue
        line = line.rstrip()
        info_list = line.split()
        ID_SAM = info_list[0]
        bitflag = int(info_list[1])#256 means secondary alignment, discard to save time
        if bitflag >=256:
            continue
#flag should not be -1, because sam id was includeed in fastq2 which is included in fastq1
        flag = fq1_ref_dict.get(ID_SAM,-1)
        if flag == -1:
            ID_error += 1
            continue
        GEXF = info_list[11]
        if GEXF.startswith('GE'):#if tagged with GE, we can get a gene count without limitation
            gene = GEXF.split(':Z:')[1]
            umi_genenumber_dict[flag][gene] += 1
        elif GEXF.startswith('XF'):#if tagged with XF, we can only get count when the aligned reigons are INTRONIC OR UTR
            if len(info_list)<=17:
                continue
            XF_info = GEXF.split(':Z:')[1]
            gf_info = info_list[16]
            gn_info = info_list[17]
            gz_info = info_list[18]
            gf_list = gf_info.split(':Z:')[1].split(',')
            gn_list = gn_info.split(':Z:')[1].split(',')
            gz_list = gz_info.split(':Z:')[1].split(',')
            possible_gene_count = 0
            temp_list = []
            for item in gf_list:
                if item =='INTRONIC' or item =='UTR':
                    possible_gene_count += 1
                    index = gf_list.index(item)
                    gfgngz = (item,gn_list[index],gz_list[index])
                    temp_list.append(gfgngz)
            if possible_gene_count < 1:# no UTR and INTRON save time
                continue
            else:
                if XF_info == 'INTRONIC' or XF_info =='UTR':#one alignment save time
                    if len(temp_list) == 1:
                        gene = temp_list[0][1]
                        umi_genenumber_dict[flag][gene] += 1
                        continue
                final_judge_list = []#multi alignment, e.g.[intron,utr,intron], we have to decide which one we used to count
                if bitflag == 0:#the final judgement depends on gz info and bitflag
                    pn_flag = '+'
                elif bitflag == 16:
                    pn_flag = '-'
                for i in temp_list:
                    if i[2] == pn_flag:
                        final_judge_list.append(i)
                if len(final_judge_list) ==1:
                    gene = final_judge_list[0][1]
                    umi_genenumber_dict[flag][gene] += 1
                # if one reads can be assigned to more then one gene, we think the assignment is ambigugous,so discard this read
                else:
                    continue
#choose the gene with the largest count as the ture assignment
umi_genenumber_dict = dict(umi_genenumber_dict)
umi_gene_dict = {}
for k,v in umi_genenumber_dict.items():
    gene_umi = max(v,key = v.get)
    umi_gene_dict[k]=gene_umi
#dge is a matrix with its columns are cells(barcodes) rows are gene counts
col_list = []
row_dict = {}
for cell,umi in cell_umi_dict.items():
    col_list.append(cell)
    for umimark in umi.keys():
        flag = umi_gene_dict.get(umimark,-1)
        if flag == -1:
            umi_discarded += 1
        else:
            dge_dict[cell][flag] += 1
            row_dict[flag] = ''

with open(DGE, 'w') as f:
    line_str = 'Gene,'
    for col in col_list:
        line_str += '{0},'.format(cell_barcode_save[col-1])
    f.write('%s\n' % line_str[:-1])
    for gene in row_dict.keys():
        line_str = '%s,' % gene
        # append every col into line
        for col in col_list:
            value_dict = dge_dict.get(col,{})
            value = value_dict.get(gene,0)
            line_str += '{0},'.format(value)
        f.write('%s\n' % line_str[:-1])  # remove ending comma and add \n
        
end_time = time.time()
run_time = float(end_time-begin_time)/60

with open(report,'w') as f:
    f.write('error occured when matching barcode with umi: ' + str(barcode_umi_match_error)+'\n')
    f.write('error occured when matching ID in SAM with reference of fastq1: ' + str(ID_error)+'\n')
    f.write('umi discarded: ' + str(umi_discarded)+'\n')
    f.write('run time: ' + str(run_time)+'\n')


