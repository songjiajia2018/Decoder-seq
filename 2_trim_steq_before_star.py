# -*- coding: utf-8 -*-

##Our experiment corrects barcodes by comparing the sequence barcode to a set of barcodes
## with the same X and Y values.If there is no exact match for a barcode in the whitelist,
## the barcode is substituted with the closest barcode from the whitelist if the edit 
## distance is not greater than 1.
import sys
import os
import re
from collections import defaultdict
import Levenshtein
import time



fastq1 = str(sys.argv[1])
fastq2 = str(sys.argv[2])
out_prefix_fq1 = fastq1.split('.')[0]
out_prefix_fq2 = fastq2.split('.')[0]
lib_barcode_path_X = str(sys.argv[3])
lib_barcode_X = []
lib_barcode_path_Y = str(sys.argv[4])
lib_barcode_Y = []
start_x = int(sys.argv[5])
end_x = int(sys.argv[6])
start_y = int(sys.argv[7])
end_y = int(sys.argv[8])
start_umi = int(sys.argv[9])
end_umi = int(sys.argv[10])

begin_time = time.time()
mismatch_x = 0
mismatch_y = 0
umi_del = 0

with open(lib_barcode_path_X, 'r') as f:
    for line in f:
        line = line.rstrip()
        barcode = str(line)
        lib_barcode_X.append(barcode)

with open(lib_barcode_path_Y, 'r') as v:
    for line in v:
        line = line.rstrip()
        barcode = str(line)
        lib_barcode_Y.append(barcode)

replace_dict = defaultdict(str)


def replace_char(string, char, index1, index2):
    return string[:index1] + char + string[index2:]


with open(fastq1, 'r') as trim_ref:
    row_count = 0
    for line in trim_ref:
        row_count += 1
        line = line.rstrip()
        line = str(line)
        if row_count%4 == 1:
            read_now = row_count
            continue
        if row_count%4 == 2:
            sim_x = end_x - start_x + 1
            sim_y = end_y - start_y + 1
            barcode_x = line[start_x:end_x]
            for tag in lib_barcode_X:
                sim = Levenshtein.distance(barcode_x, tag) #using edit distance to 
                if sim <= sim_x:
                    sim_x = sim
                    sim_x_repair = tag
            if sim_x >= 1:
                if sim_x > 1:
                    replace_dict[row_count-1] = 0
                    mismatch_x += 1
                    continue
                else:
                    replace_dict[row_count-1] = replace_dict[row_count-1] + '>'+'x:'+str(sim_x_repair)
            barcode_y = line[start_y:end_y]
            for tag in lib_barcode_Y:
                sim = Levenshtein.distance(barcode_y, tag)
                if sim <= sim_y :
                    sim_y = sim 
                    sim_y_repair = tag
            if sim_y >=1:
                if sim_y >1:
                    replace_dict[row_count-1]=0
                    mismatch_y +=1
                    continue
                else:
                    replace_dict[row_count-1] = replace_dict[row_count-1] + '>'+'y:'+str(sim_y_repair)
        if row_count%4 == 0: 
            str_umi = line[start_umi:end_umi]
            tmp=end_umi-start_umi
            for i in range(0,tmp):
                q_umi = ord(str_umi[i])-33
                if q_umi<10:
                    replace_dict[row_count-3]=0
                    umi_del +=1
                    break
reads_total_num = row_count/4
replace_dict = dict(replace_dict)
out_put_r2 = out_prefix_fq2+'_trim.fq'
out_put_r1 = out_prefix_fq1+'_trim.fq'
cr1 = 0
cr2 = 0
with open(fastq1,'r') as r1:
    r1_rowcount = 0
    with open(out_put_r1,'w') as OUTPUT:
        for line in r1:
            line = line.rstrip()
            r1_rowcount +=1
            if r1_rowcount%4 ==1:
                flag = 0
                read_now = r1_rowcount
                flag = replace_dict.get(read_now,1)
            if flag != 0:
                cr1 +=1
                if flag ==1:
                    out_str = str(line)+'\n'
                    OUTPUT.write(out_str)
                else:
                    if r1_rowcount%4 !=2:
                        out_str = str(line)+'\n'
                        OUTPUT.write(out_str)
                    else:
                        out_str = str(line)
                        temp_str = str(flag)
                        barcode_rep_list = temp_str.split('>')
                        for item in barcode_rep_list:
                            if item != '':
                                barcode_poi = str(item.split(':')[0])
                                barcode_cont = item.split(':')[1]
                                #care str.replace will replace more than one alignment if no limitation,so use new function
                                if barcode_poi == 'x':
                                    out_str = replace_char(out_str,barcode_cont,start_x,end_x)
                                if barcode_poi == 'y':
                                    out_str = replace_char(out_str,barcode_cont,start_y,end_y)
                        out_str +='\n'
                        OUTPUT.write(out_str)
            else :
                pass
#edit r2 and output
with open(fastq2,'r') as r2:
    r2_rowcount = 0
    with open(out_put_r2,'w') as OUTPUT:
        for line in r2:
            line = line.rstrip()
            r2_rowcount +=1
            if r2_rowcount%4 ==1:
                flag = 0
                read_now = r2_rowcount
                flag = replace_dict.get(read_now,1)
            if flag != 0:
                cr2 += 1
                out_str = str(line)+'\n'
                OUTPUT.write(out_str)
            else :
                pass
reads_after_1st_2nd_step = r2_rowcount/4
reads_after_1st_2nd_step_r1= r1_rowcount/4
end_time = time.time()
run_time = float(end_time-begin_time)/60

with open (out_prefix_fq1+'_trim_report.txt','w') as REPORT:
    REPORT.write('before trim: '+str(reads_total_num)+'\n')
    REPORT.write('after trim: '+str(cr1/4)+'\n')
    REPORT.write('reads discarded in step1--mismatch_X: '+str(mismatch_x)+'\n')
    REPORT.write('reads discarded in step1--mismatch_Y: '+str(mismatch_y)+'\n')
    REPORT.write('reads discarded in step2--low quality umi: '+str(umi_del)+'\n')
    REPORT.write('qc control completed,run time: '+str(run_time)+ ' min' + '\n')
    

