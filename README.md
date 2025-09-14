# Bulk RNA-seq enrichment analysis for bacterial.
This pipeline is a workflow for Bacteria Bulk RNA-seq enrichment analysis.

# Workflow

# Requirements
You can use conda to install required software/package.
```
conda install mamba -c conda forge
mamba env create -f enrichment.yml
conda activate BacEnrich
```

# Pipeline description
Here you can find a detail steps description of workflow.
## Step 1: Functional annotation of reference
Because that Klebsiella pneumoniae is non-model organism (not like Escherichia coli), it need to been annoted to get the Term to Gene table. Here we use  eggNOG-mapper to perform functional annotation. You can submit directly your .faa file into 
http://eggnog-mapper.embl.de/ to get Term to Gene table. Sometimes the online website is unavailable. You can try functional annotation on local.

### Step 1-1 Download 
You can manually download the requirement databases and decompress files
```
#Download the eggNOG database
export EGGNOG_DATA_DIR=/research/groups/ma1grp/home/zyu/my_github/enrichment_run/1_functional_annotation/database
lftp -c "open http://eggnog6.embl.de; mirror --continue --parallel=4 /download/emapperdb-5.0.2 $EGGNOG_DATA_DIR"
```

```
#Decompress files
gunzip ${EGGNOG_DATA_DIR}/*.gz
for f in ${EGGNOG_DATA_DIR}/*.tar; do
    echo "Extracting $f ..."
    tar -xf "$f" -C ${EGGNOG_DATA_DIR}
done
```
### Step 1-2 Running
```
faa_file=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.faa
emapper.py -i ${faa_file} -o Ma_L5H --tax_scope Bacteria --excel
```
If you want set other parameter, you can refer [eggNOG-mapper wiki](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13#user-content-Software_Requirements).
## Step 2: Raw reads processing
# bash
```
gff_file=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished_original.gff3
fna_file=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished.fna
raw_reas_path=/research/groups/ma1grp/home/zyu/my_github/github_vscode/enrichment_test
metadata=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/meta_data.txt
mapping_shell=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/script/mapping_bulk_paired.sh
cat ${metadata} | while read line; do
    ssr_id=$(echo $line | cut -f1 -d " ")
    input_r1=${raw_reas_path}/${ssr_id}/${ssr_id}_1.fastq.gz
    input_r2=${raw_reas_path}/${ssr_id}/${ssr_id}_2.fastq.gz
    echo "sh ${mapping_shell} -i ${input_r1} -j ${input_r2} -g ${gff_file} -f ${fna_file}"
done > running_map.sh
sh running_map.sh
```

## Step 3: FeatureCount
After running Step2, the raw reads were performed QC, aligment and sort. Finally generate *sort.bam file. This step use featureCounts to generate gene count table.
```
gff_file=/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/reference/Ma_L5H_1_polished_original.gff3
all_bam_files=$(ls *sorted.bam | tr '\n' ' ') 
featureCounts -p -d 10 -D 1000 -t CDS,ncRNA,tmRNA,tRNA,regulatory_region -g ID -a ${gff_file} -o feature.count -R BAM ${all_bam_files}  -T 4
```

## Step 4: Enrichment
