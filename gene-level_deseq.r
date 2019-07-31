library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library("pheatmap")
library("data.table")
library("ggplot2")

perform_DESeq <- function(samples_, tx2gene_){
  files <- file.path(dir, "results", paste("PA_", samples_$sample, sep = ""), "quant.sf")
  names(files) <- samples_$sample
  txi <- tximport(files, type="salmon", tx2gene=tx2gene_)
  ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples_, design = ~ group)
  dds <- DESeq(ddsTxi)
  return(dds)
}

plotgene <- function(x, dds_) {
  a <- plotCounts(dds_, gene = x, intgroup = "group", returnData = TRUE)
  write.csv(a, file=paste(x, "_normalized_count.csv", sep = ''))
  p <- ggplot(a, aes(x=group, y=count)) + geom_boxplot()
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  p <- p + scale_y_log10(breaks=c(100,200,500,1000))
  p + xlab("Group") + ylab(paste("Normalized read count of ", x))
} 

setwd("E://Project/Pituitary_Adenoma/Reperform_analysis")
dir <- "E://Project/Pituitary_Adenoma/Reperform_analysis"

samples <- read.table(file.path(dir, "group.list"), header=TRUE)
rownames(samples) <- samples$sample
samples$group <- as.character(samples$group)
samples <- samples[samples$group != "WT-3",]
samples[samples$group %like% "^WT",]$group <- "WT"

tx2gene <- read.table(file.path(dir, "GCF_000001405.25_GRCh37.p13_metadata.HGNC"))

dds <- perform_DESeq(samples, tx2gene)
plotgene("ESRRG", dds)
plotgene("PRL", dds)

samples_WT_NM <- samples[samples$group == "NM" | samples$group == "WT",]
dds_WT_NM <- perform_DESeq(samples_WT_NM, tx2gene)
res_WT_NM <- results(dds_WT_NM)

samples_MT_NM <- samples[samples$group == "NM" | samples$group == "MT",]
dds_MT_NM <- perform_DESeq(samples_MT_NM, tx2gene)
res_MT_NM <- results(dds_MT_NM)

res_WT_NM["ESRRG",]
res_MT_NM["ESRRG",]
