library(DESeq2)
library(rtracklayer)
library(dplyr)
library(ontologyIndex)
library(tidyr)
library(readxl)
library(clusterProfiler)
library(optparse)

#=================================================================
#+++++++++++++++++++++++Parameters settings ++++++++++++++++++++++
#=================================================================
option_list <- list(
    make_option(c("-c", "--count_table_file"), type = "character", default = NULL,
                            help = "Path to count table file", metavar = "character"),
    make_option(c("-m", "--metadata_file"), type = "character", default = NULL,
                            help = "Path to metadata file", metavar = "character"),
    make_option(c("-a", "--annotations_file"), type = "character", default = NULL,
                            help = "Path to eggNOG annotations file", metavar = "character"),
    make_option(c("-g", "--go_obo_file"), type = "character", default = NULL,
                            help = "Path to go.obo file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

count_table_file <- opt$count_table_file
metadata_file <- opt$metadata_file
annotations_file <- opt$annotations_file
go_obo_file <- opt$go_obo_file


# count_table_file <- "/research/groups/ma1grp/home/zyu/my_github/enrichment_run/3_featureCount/gene.count"
# metadata_file <- "/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/meta_data.txt" 

# annotations_file <- "/research/groups/ma1grp/home/zyu/my_github/enrichment_run/1_functional_annotation/Ma_L5H.emapper.annotations.xlsx" # EDIT THIS to your annotations file
# go_obo_file <- "/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/reference/go.obo" # EDIT THIS to your go.obo file

print(paste("Using count table file:", count_table_file))
print(paste("Using metadata file:", metadata_file))
print(paste("Using annotations file:", annotations_file))
print(paste("Using GO OBO file:", go_obo_file))

#=================================================================
#+++++++++++++++++++++++Step 1 Get eggNOG annotations ++++++++++++
#=================================================================
#Step 1-1 Get eggNOG annotations and process it to TERM2GENE format
eggNOG_annotations <- read_xlsx(annotations_file, skip = 2, col_names = TRUE) %>% 
    filter(GOs != "-") %>%
    dplyr::select( GOs,query) %>%
    separate_rows(GOs, sep = ",")

colnames(eggNOG_annotations) <- c("term", "gene")
dim(eggNOG_annotations)

#Step 1-2 Get GO term names from go.obo file
ontology <- get_ontology(file = go_obo_file,
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)

eggNOG_term <- eggNOG_annotations %>%
    mutate(name = ontology$name[term]) %>%
    select(c(term, name)) %>%
    distinct() %>%
    drop_na() %>%
    filter(!grepl("obsolete", name))

# Step 1-3 Filter eggNOG_annotations to keep only terms present in eggNOG_term
eggNOG_annotations <- eggNOG_annotations %>% filter(term %in% eggNOG_term$term)

dim(eggNOG_annotations)
#=================================================================
#+++++++++++++++++++++++Step 2 Count matrix input ++++++++++++++++
#=================================================================
count_table <- read.table(count_table_file, row.names = 1, header = TRUE, comment.char = "#")
colnames(count_table) <- gsub("_.*$", "", colnames(count_table)) 

meta_table <- read.table(metadata_file, sep = "\t", header = FALSE, row.names = 1)
colnames(meta_table)[1] <- "Group"

#=================================================================
#+++++++++++++++++++++++Step 3 DESeq2 analysis ++++++++++++++++++++
#=================================================================
id <- rownames(meta_table)
cts <- count_table[, id]


#Step 3-1 Check the row and col name
check_row_col <- all(rownames(meta_table) == colnames(cts))
print(paste("Row names in coldata match column names in cts:",as.character(check_row_col)))

#Step 3-2 Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = meta_table,
                              design = ~ Group)

# #Step 3-3 Filter low count genes    
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
# dds <- dds[keep,]

#Step 3-4 Run DESeq2
dds <- DESeq(dds)

#Step 3-5 Get results
group1 <- unique(meta_table$Group)[1]
group2 <- unique(meta_table$Group)[2]

res <- as.data.frame(results(dds,contrast = c("Group", group1, group2)))

write.csv(res, file = "DESeq2_results.csv")

#=================================================================
#+++++++++++++++++++++++Step 4 Enrichment analysis +++++++++++++++
#=================================================================
#Step 4-1: Transfer the DESeq2 results to GSEA input
tab_no_na <- res %>% filter(!is.na(pvalue)) # Filter out rows with NA in pvalue

gsea_list <- tab_no_na %>% dplyr::arrange(desc(log2FoldChange))

gsea_input <- gsea_list %>%
    dplyr::select(log2FoldChange) %>%
    unlist() %>%
    as.vector()

names(gsea_input) <- rownames(gsea_list)

#Step 4-2: GSEA analysis below
enrichment_gsea <- GSEA(geneList = gsea_input,
                        TERM2GENE = eggNOG_annotations,
                        TERM2NAME = eggNOG_term,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")

#Step 4-3: save the enrichment result
gsea_result <- enrichment_gsea@result
write.csv(file = "GSEA_result.csv", x = gsea_result)

