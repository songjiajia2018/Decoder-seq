# Decoder-seq
Preprocessing of Decoder-seq raw data

Starting with FASTAQ files after QC, this pipeline analyses the data to produce a count matrix that 
represents the predicted number of unique molecules derived from each gene at each spot. 
The xy_index.txt file is the barcode and its coordinates.

run: python run_pipeline.py --R1 _1.fq \
							--R2 _2.fq \
							--map hg38 \
							--bcx lib_barcode_X.txt \
							--bcy lib_barcode_Y.txt \
							--startx 0 \
							--endx 8 \
							--starty 38 \
							--endy 46 \
							--startumi 46 \
							--endumi 58
										
												
This particular script requires the following softwares and versions:
Python 3.8.8
R 4.1.1 
Nucleotide-Nucleotide BLAST 2.10.1+
STAR 2.7.3a							
starcode v1.3
samtools 1.9
Drop-seq_tools 2.3.0 
