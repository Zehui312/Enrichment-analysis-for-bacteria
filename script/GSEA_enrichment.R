library(DESeq2)
library(rtracklayer)
library(dplyr)
library(pheatmap)
library(ontologyIndex)
library(tidyr)
library(readr)
library(clusterProfiler)
library(optparse)


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


count_table_file <- "/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/count_table"
metadata_file <- "/research/groups/ma1grp/home/zyu/my_github/Enrichment-analysis-for-bacteria/meta_data.txt" 

annotations_file <- "/research/groups/ma1grp/home/zyu/work_2025/scRNA_7_scRNA/learning_new_tools/eggNOG/bac/annotations_bac.tsv" # EDIT THIS to your annotations file
go_obo_file <- "/research/groups/ma1grp/home/zyu/work_2025/scRNA_7_scRNA/Bulk_RNA_seq_pipeline/develop_enrichment/go.obo" # EDIT THIS to your go.obo file

print(paste("Using count table file:", count_table_file))
print(paste("Using metadata file:", metadata_file))
print(paste("Using annotations file:", annotations_file))
print(paste("Using GO OBO file:", go_obo_file))
#=================================================================
#+++++++++++++++++++++++Step 1 Get eggNOG annotations ++++++++++++
#=================================================================

#Step 0-1 Get eggNOG annotations
eggNOG <- read_tsv(annotations_file) %>%
    dplyr::select(GOs, `#query`) %>% # Select the relevant columns
    dplyr::filter(GOs != "-") %>% # Filter out rows with no GO terms annotation
    separate_rows(GOs, sep = ",") %>% # Separate multiple GO terms into different rows
    mutate(gene = gsub("\\..*", "", `#query`)) %>% # Extract gene names by removing suffixes
    select(GOs, gene) %>% # Select relevant columns
    distinct() %>% # Remove duplicate entries
    drop_na() # Drop rows with NA values
colnames(eggNOG) <- c("term", "gene")

eggNOG <- eggNOG %>% filter(term != "NA") # Filter out rows
Go_num <- length(eggNOG$term)
locus_tab_num <- length(unique(eggNOG$gene))
print(paste("The number of unique GO terms:", Go_num))
print(paste("The number of unique genes:", locus_tab_num))

#Step 0-2 Preparing Term to name
ontology <- get_ontology(file = go_obo_file,
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)
eggNOG_term <- eggNOG %>%
    mutate(name = ontology$name[term]) %>%
    select(c(term, name)) %>%
    distinct() %>%
    drop_na() %>%
    filter(!grepl("obsolete", name))

eggNOG <- eggNOG %>% filter(term %in% eggNOG_term$term)
Go_num_filtered <- length(eggNOG$term)
locus_tab_num_filtered <- length(unique(eggNOG$gene))

print(paste("The number of unique genes after filtering:", locus_tab_num_filtered))
print(paste("The number of unique GO terms after filtering:", Go_num_filtered))


#=================================================================
#+++++++++++++++++++++++Step 2 Count matrix input ++++++++++++++++
#=================================================================
count_table <- read.table(count_table_file, row.names = 1, header = TRUE, comment.char = "#")
colnames(count_table) <- gsub("_.*$", "", colnames(count_table)) 

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
rownames(metadata) <- paste0("X", metadata$SampleName)  # Ensure row names match the count table




#=================================================================
#+++++++++++++++++++++++Step 3 DESeq2 analysis ++++++++++++++++++++
#=================================================================
id <- rownames(metadata)
cts <- count_table[, id]
coldata <- metadata[id, ]

#Step 3-1 Check the row and col name
check_row_col <- all(rownames(coldata) == colnames(cts))
print(paste("Row names in coldata match column names in cts:",as.character(check_row_col)))

#Step 3-2 Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)

#Step 3-3 Filter low count genes    
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Step 3-4 Run DESeq2
dds <- DESeq(dds)


rld <- rlog(dds, blind=FALSE)
rld_assay <- assay(rld)

#Step 3-5 Get results
group1 <- unique(coldata$Group)[1]
group2 <- unique(coldata$Group)[2]

res <- results(dds,contrast = c("Group", group1, group2))
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC <- lfcShrink(dds, coef=2, type="ashr")  # or "ashr"
resLFC$log10_padj <- -log10(resLFC$padj)

sample_name <- paste0(group1, "_vs_", group2)
res_tab <- as.data.frame(res)

write.csv(res_tab, file = paste0("DESeq2_results_", sample_name, ".csv"))
write.csv(resLFC, file = paste0("DESeq2_results_LFC_", sample_name, ".csv"))
write.csv(rld_assay, file = paste0("rlog_assay_", sample_name, ".csv"))
#=================================================================
#+++++++++++++++++++++++Step 4 Enrichment analysis +++++++++++++++
#=================================================================
#Step 4-1: Transfer the DESeq2 results to GSEA input
tab_no_na <- res_tab %>% filter(!is.na(pvalue)) # Filter out rows with NA in pvalue

gsea_list <- tab_no_na %>% dplyr::arrange(desc(log2FoldChange))

gsea_input <- gsea_list %>%
    dplyr::select(log2FoldChange) %>%
    unlist() %>%
    as.vector()

names(gsea_input) <- rownames(gsea_list)

#Step 4-2: GSEA analysis below
enrichment_gsea <- GSEA(geneList = gsea_input,
                        TERM2GENE = eggNOG,
                        TERM2NAME = eggNOG_term,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")

#Step 4-3: save the enrichment result
gsea_result <- enrichment_gsea@result
write.csv(file = paste0("GSEA_", sample_name, ".csv"),x = gsea_result)
#=================================================================
#+++++++++++++++++++++++Step 6 Lollipop plot +++++++++++++++++++++
#=================================================================
library(ggplot2)
library(cols4all)
top_num <- 10

GSEA_df <- arrange(gsea_result, desc(NES))
GSEA_topN <- rbind(GSEA_df[1:top_num,], GSEA_df[(nrow(GSEA_df)-top_num + 1):nrow(GSEA_df),])

GSEA_topN$abs_NES <- abs(GSEA_topN$NES)
GSEA_topN$log10_padj <- -log10(GSEA_topN$p.adjust)
GSEA_topN <- arrange(GSEA_topN, abs_NES)

GSEA_up_sig <-subset(GSEA_topN, NES > 0)
GSEA_down_sig <- subset(GSEA_topN, NES < 0)
GSEA_down_sig$log10_padj <- -GSEA_down_sig$log10_padj

GSEA_tab <- rbind(GSEA_up_sig, GSEA_down_sig)

level_order <- c(GSEA_down_sig$Description,GSEA_up_sig$Description)
GSEA_tab$Description <- factor(GSEA_tab$Description, levels = level_order)
GSEA_tab$tags_percent <- as.numeric(sub(".*tags=([0-9]+)%.*", "\\1", GSEA_tab$leading_edge))

p <- ggplot(GSEA_tab, aes(x = NES, y = Description)) +
    geom_col(aes(fill = log10_padj), width = 0.02) +
    geom_point(aes(size = tags_percent, color = log10_padj)) +
    scale_fill_continuous_c4a_div('sunset', mid = 0,limits = c(-5, 5),oob = scales::squish) +
    scale_color_continuous_c4a_div('sunset', mid = 0,limits = c(-5, 5),oob = scales::squish) +
    ggtitle(paste0("GSEA enrichment (TOP10): ", sample_name)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.line.x = element_line(color = 'grey60', size = 0.8)) +
          geom_vline(xintercept = 0, size = 0.8, color = 'grey60', lty = "solid") +
          scale_x_continuous(breaks = seq(-2, 2, by = 1),
                             limits = c(-3.5, 3.5),
                             labels = seq(-2, 2, by = 1)) +
          geom_text(data = GSEA_up_sig,aes(x = -0.1, y = Description, label = Description),size = 4.5,hjust = 1) +
          geom_text(data = GSEA_down_sig,aes(x = 0.1, y = Description, label = Description),size = 4.5,hjust = 0)

ggsave(filename = paste0("GSEA_lollipop_fc", sample_name, ".adjust.pdf"), dpi = 300, width = 35, height = 20, units = "cm")

