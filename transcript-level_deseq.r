#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("DESeq2")
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library("pheatmap")
library("data.table")
library("ggplot2")
library("biomaRt")

perform_DESeq <- function(samples_){
  files <- file.path(dir, "results", paste("PA_", samples_$sample, sep = ""), "quant.sf")
  names(files) <- samples_$sample
  txi <- tximport(files, type="salmon", txOut = TRUE)
  ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples_, design = ~ group)
  dds <- DESeq(ddsTxi)
  return(dds)
}

plotgene <- function(x, dds_) {
  a <- plotCounts(dds_, gene = x, intgroup = "group", returnData = TRUE)
  write.csv(a, file = file.path("transcripts", paste(x, "_normalized_count.csv", sep = '')))
  
  p <- ggplot(a, aes(x=group, y=count)) + geom_boxplot()
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  p <- p + scale_y_log10()
  p <- p + xlab("Group") + ylab(paste("Normalized read count of ", x))
  ggsave(file = file.path("transcripts", paste(x, "_plot.png", sep = '')),
         plot = p,
         width = 18, 
         height = 12,
         units = "cm")
} 

setwd("E://Project/Pituitary_Adenoma/Reperform_analysis")
dir <- "E://Project/Pituitary_Adenoma/Reperform_analysis"

samples <- read.table(file.path(dir, "group.list"), header=TRUE)
rownames(samples) <- samples$sample
samples$group <- as.character(samples$group)
samples <- samples[samples$group != "WT-3",]
samples[samples$group %like% "^WT",]$group <- "WT"

dds <- perform_DESeq(samples)

samples_WT_NM <- samples[samples$group == "NM" | samples$group == "WT",]
dds_WT_NM <- perform_DESeq(samples_WT_NM)
res_WT_NM <- results(dds_WT_NM)

samples_MT_NM <- samples[samples$group == "NM" | samples$group == "MT",]
dds_MT_NM <- perform_DESeq(samples_MT_NM)
res_MT_NM <- results(dds_MT_NM)

txofESRRG <- read.csv(file.path(dir, "ESRRG_all_tx_Refseq.csv"), header = FALSE)

info_WT_NM <- data.frame(res_WT_NM[txofESRRG$V1, ]@listData[c("baseMean", "log2FoldChange", "pvalue", "padj")], 
                         row.names = res_WT_NM[txofESRRG$V1, ]@rownames)

info_MT_NM <- data.frame(res_MT_NM[txofESRRG$V1, ]@listData[c("baseMean", "log2FoldChange", "pvalue", "padj")],
                         row.names = res_MT_NM[txofESRRG$V1, ]@rownames)

info_MT_NM$log2FoldChange <- -info_MT_NM$log2FoldChange

for(tx in txofESRRG$V1) {
  plotgene(tx, dds)
}

write.csv(info_WT_NM, file = file.path("transcripts", paste("WTvsNM_diff.csv", sep = '')))
write.csv(info_MT_NM, file = file.path("transcripts", paste("MTvsNM_diff.csv", sep = '')))
