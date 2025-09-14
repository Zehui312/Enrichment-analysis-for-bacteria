
#!/bin/bash
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Basic Information: Bulk RNA Sequence mapping
# 2. Author: Zehui Yu
# 3. Version: 1.0
# 4. Created data: 2024/12/24
# 5. Usage: sh umi_tools_mapping_bulk.sh -i SRR15174659_1.fastq.gz -j SRR15174659_2.fastq.gz
# 6. Parameter Explanation
# -i input_r1: Input file 1
# -j input_r2: Input file 2
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#=================================================================
#++++++++++++++++++ Parse command line arguments++++++++++++++++++
#=================================================================
echo "usage: umi_tools_mapping_bulk.sh -i <input_r1> -j <input_r2> -g <gff_file> -f <fna_file>. Example: umi_tools_mapping_bulk.sh -i SRR15174659_1.fastq.gz -j SRR15174659_2.fastq.gz -g Ma_L5H_1_polished.gff3 -f Ma_L5H_1_polished.fna"

while getopts i:j:f: flag
do
    case "${flag}" in
        i) input_r1=${OPTARG};;
        j) input_r2=${OPTARG};;
        f) fna_file=${OPTARG};;
    esac
done
sample_name=$(basename ${input_r1} | cut -d '.' -f 1)
echo "sample_name: ${sample_name}"
echo "input_r1: ${input_r1}"
echo "input_r2: ${input_r2}"
echo "fna_file: ${fna_file}"
echo "==========================================================="
echo "Start mapping for sample: ${sample_name}"
echo "==========================================================="
#=================================================================
#+++++++++++++++++++++++Step 1 Mapping +++++++++++++++++++++++++++
#=================================================================

#1 QC with fastp
fastp -i ${input_r1} -o ${sample_name}.R1.qc.fq -I ${input_r2} -O ${sample_name}.R2.qc.fq -l 30 -j ${sample_name}.qc.json -h ${sample_name}.qc.html -f 5 -t 5

#1.Indexing the Reference Genome 
#bwa index ./Ma_L5H_1_polished.fna

#2 Align extracted reads with bwa
bwa mem -t 8 ${fna_file} ${sample_name}.R1.qc.fq ${sample_name}.R2.qc.fq > ${sample_name}.sam

#3 Convert BAM and sort
samtools view -b -S ${sample_name}.sam | samtools sort -o ${sample_name}.sorted.bam 

