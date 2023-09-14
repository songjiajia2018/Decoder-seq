# Preprocessing of Decoder-seq raw data

Starting with FASTQ files after QC, this pipeline analyses the data to produce a count matrix that 
represents the predicted number of unique molecules derived from each gene at each spot. 


Usage for inference:

Input comprising pair-end fastq files and a whitelist of barcodes for X and Y. 
Utilising the STAR programme to generate genome indices for the reference file.
Our sequence library is pair-end 150bp. 
Read 2 contains cDNA. For read1, beginning at 5', 1-8 represents barcodeX, 39-46 represents barcodeY, and 47-58 represents UMIs. 
Both barcodes X and Y contain 50 unique 8bp sequences in the same order.	
The filename for the output count matrix is ALL_XQ2_sub_2_trim_tagged_dge.txt.	
The barcode and its coordinates are contained in the xy_index.txt file for further plotting and analysis.
```
python run_pipeline.py --R1 ALL_XQ2_sub_1.fq.gz \
					   --R2 ALL_XQ2_sub_2.fq.gz \
					   --map hg38 \
					   --bcx lib_barcode_X.txt \
					   --bcy lib_barcode_Y.txt \
					   --startx 0 \
					   --endx 8 \
					   --starty 38 \
					   --endy 46 \
					   --startumi 46 \
					   --endumi 58
										
Options for inference:
--map	    choose the reference from 'hg38' and 'mm10', default = 'mm10'
--thread    Specify the number of threads for qc step and dge step. The recommended number is less than 5, default = 3
--R1	    input R1 fastq
--R2	    input R2 fastq
--bcx	    input barcode white list for X, default is 'lib_barcode_X.txt'
--bcy	    input barcode white list for Y, default is 'lib_barcode_Y.txt'
--startx    position of barcode X, default is 0
--endx	    position of barcode X, default is 8
--starty    position of barcode Y, default is 38
--endy	    position of barcode Y, default is 46
--startumi  position of umi, default is 46
--endumi    position of umi, default is 58
```
			
												
This script requires the following softwares and versions:
Python 3.9.7
Nucleotide-Nucleotide BLAST 2.10.1+
STAR 2.7.3a							
starcode v1.3
samtools 1.9
Drop-seq_tools 2.3.0 

Author: Jia Song, Weizhou Qian, Di Sun
