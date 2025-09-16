
library(ggplot2)
library(dplyr)
library(cols4all)


library(optparse)
#=================================================================
#+++++++++++++++++++++++Parameters settings ++++++++++++++++++++++
#=================================================================
option_list <- list(
    make_option(c("-g", "--gsea_result_path"), type = "character", default = "/research/groups/ma1grp/home/zyu/my_github/enrichment_run/4_enrichment/GSEA_result.csv",
                            help = "Path to GSEA result CSV file", metavar = "character"),
    make_option(c("-n", "--topnum"), type = "integer", default = 5,
                            help = "Number of top and bottom NES to plot", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

gsea_path <- opt$gsea_result_path
top_num <- opt$topnum

print(paste("Using GSEA result file:", gsea_path))
print(paste("Number of top and bottom NES to plot:", top_num))

#=================================================================
#+++++++++++++++++++++++Visualization ++++++++++++++++++++++++++++
#=================================================================

gsea_result <- read.csv(file = gsea_path,row.names = 1)
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
    scale_fill_continuous_c4a_div('sunset', mid = 0, limits = c(-5, 5), oob = scales::squish) +
    scale_color_continuous_c4a_div('sunset', mid = 0, limits = c(-5, 5), oob = scales::squish) +
    ggtitle(paste("GSEA enrichment (TOP", top_num, " based on NES)")) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.line.x = element_line(color = 'grey60', linewidth = 0.8)) +   # fixed
    geom_vline(xintercept = 0, linewidth = 0.8, color = 'grey60', lty = "solid") +  # fixed
    scale_x_continuous(breaks = seq(-2, 2, by = 1),
                       limits = c(-3.5, 3.5),
                       labels = seq(-2, 2, by = 1)) +
    geom_text(data = GSEA_up_sig, aes(x = -0.1, y = Description, label = Description), size = 6, hjust = 1) +
    geom_text(data = GSEA_down_sig, aes(x = 0.1, y = Description, label = Description), size = 6, hjust = 0)

ggsave(filename = "GSEA_lollipop.jpg", dpi = 300, width = 35, height = 15, units = "cm")




