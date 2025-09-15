# Bulk RNA-seq enrichment analysis for bacterial.
This repository provides a bulk RNA-seq enrichment analysis pipeline for bacteria. All steps of the workflow are included in the shell script (`./Workflow.sh`).
- [Requirements](#Requirements)
  -  [1. Installing software](#1.-Installing-software-using-conda)
  -  [2. Download example fastq files](#2.-Download-example-fastq-files)
- [Pipeline description](#pipeline-description)
  -  [Step 1: Functional annotation of reference](#step-1-functional-annotation-of-reference)
  -  [Step 2: Raw reads processing](#step-2-raw-reads-processing)
  -  [Step 3: Gene counts generation](#step-3-gene-counts-generation)
  -  [Step 4: Enrichment analysis](#step-4-enrichment-analysis)
  -  [Step 5: Visualization](#step-5-visualization)
# Workflow
The schematic of our workflow is demonstrated below.

<img src="/images/Workflow.png" alt="Workflow" width="500"/>

# Requirements
## 1. Installing software using conda
The required software and packages can be easily installed using Conda. 
```
git clone https://github.com/Zehui312/Enrichment-analysis-for-bacteria.git
cd Enrichment-analysis-for-bacteria
conda env create -f enrichment.yml
conda activate BacEnrich
```
## 2. Download example fastq files
Here, we use 6 sequencing data as example, they are untreated and meropenem-treated and replicated 3 time. You can find the **meta_data.txt** in 
```
meta_data=./Enrichment-analysis-for-bacteria/meta_data.txt
cat ${meta_data} | while read line; do
    sample_name=$(echo $line | cut -f2 -d " ")
    srr_id=$(echo $line | cut -f1 -d " ")
    echo "fastq-dump --gzip --split-3 -O ${srr_id} -A ${srr_id}"
done > download_data.sh
sh download_data.sh
```
# Pipeline description
Here you can find a detail steps description of workflow.
## Step 1: Functional annotation of reference
Because that Klebsiella pneumoniae is non-model organism (not like Escherichia coli), it need to been annoted to get the Term to Gene table. Here we use  eggNOG-mapper to perform functional annotation. You can submit directly your .faa file into 
http://eggnog-mapper.embl.de/ to get Term to Gene table. Sometimes the online website is unavailable. You can try functional annotation on local.

### Step 1-1 Download the eggNOG database
You can manually download the requirement databases and decompress files
```
#Download 
export EGGNOG_DATA_DIR=./enrichment_run/1_functional_annotation/database
lftp -c "open http://eggnog6.embl.de; mirror --continue --parallel=4 /download/emapperdb-5.0.2 $EGGNOG_DATA_DIR"

#Decompress files
gunzip ${EGGNOG_DATA_DIR}/*.gz
for f in ${EGGNOG_DATA_DIR}/*.tar; do
    echo "Extracting $f ..."
    tar -xf "$f" -C ${EGGNOG_DATA_DIR}
done
```
### Step 1-2 Running eggNOG
Running eggNOG would takes severa hours, you'd better backgroud running using nohup. 
**Ma_L5H_1_polished.faa** files, see here: `./reference/Ma_L5H_1_polished.faa `
```
faa_file=./Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.faa
emapper.py -i ${faa_file} -o Ma_L5H --tax_scope Bacteria --excel
```
If you want set other parameter, you can refer [eggNOG-mapper wiki](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13#user-content-Software_Requirements).
## Step 2: Raw reads processing
```./script/mapping_bulk_paired.sh``` can achieve Raw data QC and Alignment. You can just simplely run below script. 
- **mapping_shell** : `./reference/mapping_bulk_paired.sh `
- **gff_file** and **fna_file**: `./reference/Ma_L5H_1_polished_original.gff3 ` and `./reference/Ma_L5H_1_polished.fna`
- **raw_reas_path**: The path you download example fastq files
- **metadata**: `./meta_data.txt`


```
# bash
# you need to edit script and paths to your own path
mapping_shell=./Enrichment-analysis-for-bacteria/script/mapping_bulk_paired.sh
gff_file=./Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished_original.gff3
fna_file=./Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.fna
raw_reas_path=./github_vscode/enrichment_test
metadata=./Enrichment-analysis-for-bacteria/meta_data.txt


# generate mapping commands for each sample
cat ${metadata} | while read line; do
    ssr_id=$(echo $line | cut -f1 -d " ")
    input_r1=${raw_reas_path}/${ssr_id}/${ssr_id}_1.fastq.gz
    input_r2=${raw_reas_path}/${ssr_id}/${ssr_id}_2.fastq.gz
    echo "sh ${mapping_shell} -i ${input_r1} -j ${input_r2} -g ${gff_file} -f ${fna_file}"
done > running_map.sh

# run mapping commands 
sh running_map.sh
```

## Step 3: Gene counts generation
After running Step2, the raw reads were performed QC, aligment and sort. This step use featureCounts to generate gene count table.
```
gff_file=./Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished_original.gff3
all_bam_files=$(ls *sorted.bam | tr '\n' ' ') # *sorted.bam from Step 2 output
featureCounts -p -d 10 -D 1000 -t CDS,ncRNA,tmRNA,tRNA,regulatory_region -g ID -a ${gff_file} -o feature.count -R BAM ${all_bam_files}  -T 4
```

## Step 4: Enrichment analysis
The go.obo file can download from [Gene ontology](https://geneontology.org/docs/download-ontology/)

```
#you need to edit script and paths to your own path
enrichment_script=./enrichment_run/4_enrichment/GSEA_enrichment.R # path to GSEA_enrichment script (from ./reference/GSEA_enrichment.R)
count_table_file=${output_path}/3_featureCount/gene.count # path to count table file (from Step 3 output) 
annotations_file=${output_path}/1_functional_annotation/Ma_L5H.emapper.annotations.xlsx # path to functional annotation file (from Step 1 output)
go_obo_file=./Enrichment-analysis-for-bacteria/reference/go.obo  # path to go.obo file

Rscript ${enrichment_script}  \
--count_table_file ${count_table_file}  \
--metadata_file ${meta_data}  \
--annotations_file ${annotations_file}  \
--go_obo_file ${go_obo_file}
```
## Step 5: Visualization

![GSEA_lollipop](/images/GSEA_lollipop.jpg)
