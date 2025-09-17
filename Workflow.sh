# !/bin/bash
# Change to your own output path
output_path=/research/groups/ma1grp/home/zyu/run_enrichment
cd ${output_path}

#=================================================================
#+++++++++++++++++++++++ Requirements 1 ++++++++++++++++++++++++++
#=================================================================

git clone https://github.com/Zehui312/Enrichment-analysis-for-bacteria.git
cd Enrichment-analysis-for-bacteria
mamba env create -f enrichment.yml
conda activate BacEnrich


#=================================================================
#+++++++++++++++++++++++ Requirements 2 ++++++++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/0_data
cd ${output_path}/0_data

meta_data=${output_path}/Enrichment-analysis-for-bacteria/meta_data.txt
cat ${meta_data} | while read line; do
    sample_name=$(echo $line | cut -f2 -d " ")
    srr_id=$(echo $line | cut -f1 -d " ")
    echo "fastq-dump --gzip --split-3 -O ${srr_id} -A ${srr_id}"
done > running_download.sh

sh running_download.sh
# count=1
# while read runcode; do
#     bsub -P cb_${count} -J cb_${count} -n 2 -R "rusage[mem=8GB]" -eo cb_${count}.err -oo cb_${count}.out $runcode
#     count=$((count + 1))
# done < running_download.sh
#=================================================================
#+++++++++++++++++++++++ Step 1 Functional annotation ++++++++++++
#=================================================================
mkdir -p ${output_path}/1_functional_annotation
cd ${output_path}/1_functional_annotation

#Step 1-1 Download the eggNOG database and decompress the files
EGGNOG_DATA_DIR=${output_path}/1_functional_annotation/database/emapperdb-5.0.2 
mkdir -p $EGGNOG_DATA_DIR
wget -r -np -nH --cut-dirs=1 -c -P $EGGNOG_DATA_DIR http://eggnog6.embl.de/download/emapperdb-5.0.2/

# bsub -P download -J download -n 16 -R "rusage[mem=4GB]" -eo download.err -oo download.out "
# wget -r -np -nH --cut-dirs=1 -c -P $EGGNOG_DATA_DIR http://eggnog6.embl.de/download/emapperdb-5.0.2/"
# -r → recursive download
# -np → no parent (don’t go above /download/emapperdb-5.0.2/)
# -nH → no host directory (don’t create eggnog6.embl.de/)
# --cut-dirs=1 → removes the first directory (download/) so you only get emapperdb-5.0.2/
# -c → continue partial downloads (like --continue in lftp)
# -P $EGGNOG_DATA_DIR → download into your chosen folder

cd $EGGNOG_DATA_DIR/emapperdb-5.0.2/
gunzip *.gz
for f in *.tar; do
    echo "Extracting $f ..."
    tar -xf "$f" -C ${EGGNOG_DATA_DIR}
done

#Step 1-2 Running eggNOG
cd ${output_path}/1_functional_annotation

export EGGNOG_DATA_DIR=${output_path}/1_functional_annotation/database/emapperdb-5.0.2
faa_file=${output_path}/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.faa

#you'd better backgroud running using nohup
# nohup emapper.py -i ${faa_file} -o Ma_L5H --tax_scope Bacteria --excel &
emapper.py -i ${faa_file} -o Ma_L5H --tax_scope Bacteria --excel

# bsub -P egg -J egg -n 16 -R "rusage[mem=4GB]" -eo egg.err -oo egg.out "
# emapper.py -i ${faa_file} -o Ma_L5H --tax_scope Bacteria --excel"


#=================================================================
#+++++++++++++++++++++++Step 2 Raw reads processing ++++++++++++++
#=================================================================
# bash
mkdir -p ${output_path}/2_mapping
cd ${output_path}/2_mapping

raw_reas_path=${output_path}/0_data
fna_file=${output_path}/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.fna
metadata=${output_path}/Enrichment-analysis-for-bacteria/meta_data.txt
mapping_shell=${output_path}/Enrichment-analysis-for-bacteria/script/mapping_bulk_paired.sh

cat ${metadata} | while read line; do
    ssr_id=$(echo $line | cut -f1 -d " ")
    input_r1=${raw_reas_path}/${ssr_id}/${ssr_id}_1.fastq.gz
    input_r2=${raw_reas_path}/${ssr_id}/${ssr_id}_2.fastq.gz
    echo "sh ${mapping_shell} -i ${input_r1} -j ${input_r2} -f ${fna_file}"
done > running_map.sh

sh running_map.sh
# count=1
# while read runcode; do
#     bsub -P mapping_${count} -J mapping_${count} -n 8 -R "rusage[mem=8GB]" -eo mapping_${count}.err -oo mapping_${count}.out $runcode
#     count=$((count + 1))
# done < running_map.sh

#=================================================================
#+++++++++++++++++++++++Step 3 Gene counts generation ++++++++++++
#=================================================================
mkdir -p ${output_path}/3_featureCount
cd ${output_path}/3_featureCount

ln -s ${output_path}/2_mapping/*sorted.bam ./
gff_file=${output_path}/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished_original.gff3
all_bam_files=$(ls *sorted.bam | tr '\n' ' ') 
featureCounts -p -d 10 -D 1000 -t CDS,ncRNA,tmRNA,tRNA,regulatory_region -g ID -a ${gff_file} -o gene.count -R BAM ${all_bam_files}  -T 4

# bsub -P featureCount -J featureCount -n 8 -R "rusage[mem=16GB]" -eo featureCount.err -oo featureCount.out "
# featureCounts -p -d 10 -D 1000 -t CDS,ncRNA,tmRNA,tRNA,regulatory_region -g ID -a ${gff_file} -o gene.count -R BAM ${all_bam_files}  -T 4
# "

#=================================================================
#+++++++++++++++++++++++Step 4 Enrichment analysis +++++++++++++++
#=================================================================
mkdir -p ${output_path}/4_enrichment
cd ${output_path}/4_enrichment


enrichment_script=${output_path}/Enrichment-analysis-for-bacteria/script/GSEA_enrichment.R # path to GSEA_enrichment script ()
count_table_file=${output_path}/3_featureCount/gene.count # path to count table file (from Step 3 output) 
annotations_file=${output_path}/1_functional_annotation/Ma_L5H.emapper.annotations.xlsx # path to functional annotation file (from Step 1 output)
# annotations_file=/research/groups/ma1grp/home/zyu/my_github/enrichment_run/1_functional_annotation/Ma_L5H.emapper.annotations.xlsx
go_obo_file=${output_path}/Enrichment-analysis-for-bacteria/reference/go.obo  # path to go.obo file
meta_data=${output_path}/Enrichment-analysis-for-bacteria/meta_data.txt # path to metadata file
Rscript ${enrichment_script}  --count_table_file ${count_table_file}  --metadata_file ${meta_data}  --annotations_file ${annotations_file}  --go_obo_file ${go_obo_file}

#=================================================================
#+++++++++++++++++++++++Step 5 Visualization +++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/5_visualization
cd ${output_path}/5_visualization

Visualization_script=${output_path}/Enrichment-analysis-for-bacteria/script/Visualization.R
gsea_result_path=${output_path}/4_enrichment/GSEA_result.csv

top_num=5
Rscript ${Visualization_script} --gsea_result_path ${gsea_result_path} --topnum ${top_num}

