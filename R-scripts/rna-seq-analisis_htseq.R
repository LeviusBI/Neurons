library(tximport)
library(txdbmaker)
library(DESeq2)
library(readr)
library(GenomicFeatures)
library(GenomeInfoDbData)
library(ggplot2)
library(factoextra)
library(UpSetR)
library(dplyr)
library(stringr)

# Get parameters from Snakemake or environment variables
if (exists("snakemake")) {
  # Running from Snakemake directly
  work_dir <- snakemake@params[["work_dir"]]
  metadata_file <- snakemake@params[["metadata_file"]]
} else {
  # Running from Docker container (via environment variables)
  work_dir <- Sys.getenv("WORK_DIR")
  metadata_file <- Sys.getenv("METADATA_FILE")
}

# Read metadata
meta_data <- read.csv(metadata_file, sep = "\t", header = TRUE)
  
meta_data$File <- meta_data$sampleid
meta_data$File <- paste0(meta_data$File, ".counts")
  
meta_data <- meta_data %>% relocate(File, .after = sampleid)
 
# Get path to htseq directory
htseq_dir <- file.path(work_dir, "htseq")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data,
                                         directory = htseq_dir,
                                         design = ~ individual + diff_time)
  
rownames(ddsHTSeq) <- sapply(str_split(rownames(ddsHTSeq), "\\."), function(x) x[1])
 
keep <- rowSums(counts(ddsHTSeq) >= 10) > ncol(ddsHTSeq)*0.9
ddsHTSeq <- ddsHTSeq[keep,]


vsd <- vst(ddsHTSeq, blind = F)
vsd <- assay(vsd)

PCA_vst <- prcomp(t(vsd))
PCA_vst.points <- as.data.frame(PCA_vst$x)
PCA_vst.points <- merge(meta_data, cbind(rownames(PCA_vst.points), PCA_vst.points), by = 1)

PCA_vst_plot <- ggplot(PCA_vst.points, aes(PC1, PC2, col = diff_time, shape = individual))+
    geom_point(size = 2.75)+
    geom_text(size = 3)+
    theme_bw()+
    theme(legend.position = "right")+
    scale_color_brewer(palette = "Set1")+
    xlab("PC1")+
    ylab("PC2")

prop_vst_plot <- fviz_eig(PCA_vst, col.var="blue") + ggtitle("Proportion of variance")

ddsHTSeq <- DESeq(ddsHTSeq)

res2week <- results(ddsHTSeq, contrast = c("diff_time", "2 week neuron", "NPC"))
res4week <- results(ddsHTSeq, contrast = c("diff_time", "4 week neuron", "NPC"))
res6week <- results(ddsHTSeq, contrast = c("diff_time", "6 week neuron", "NPC"))
res8week <- results(ddsHTSeq, contrast = c("diff_time", "8 week neuron", "NPC"))

get_DEG <- function(res){
    DEG <- as.data.frame(res)
    DEG <- DEG[!is.na(DEG$padj),]
    return(DEG)
}

deg2week <- get_DEG(res2week)
deg4week <- get_DEG(res4week)
deg6week <- get_DEG(res6week)
deg8week <- get_DEG(res8week)

deg2week.filt <- deg2week[abs(deg2week$log2FoldChange) > 2 & deg2week$padj < 0.01,]
deg4week.filt <- deg4week[abs(deg4week$log2FoldChange) > 2 & deg4week$padj < 0.01,]
deg6week.filt <- deg6week[abs(deg6week$log2FoldChange) > 2 & deg6week$padj < 0.01,]
deg8week.filt <- deg8week[abs(deg8week$log2FoldChange) > 2 & deg8week$padj < 0.01,]

deg2week.filt <- deg2week.filt[order(deg2week.filt$log2FoldChange, decreasing = T),]
deg4week.filt <- deg4week.filt[order(deg4week.filt$log2FoldChange, decreasing = T),]
deg6week.filt <- deg6week.filt[order(deg6week.filt$log2FoldChange, decreasing = T),]
deg8week.filt <- deg8week.filt[order(deg8week.filt$log2FoldChange, decreasing = T),]


list.upregulation <- list(week_2 = rownames(deg2week.filt[deg2week.filt$log2FoldChange > 0,]), 
     week_4 = rownames(deg4week.filt[deg4week.filt$log2FoldChange > 0,]), 
     week_6 = rownames(deg6week.filt[deg6week.filt$log2FoldChange > 0,]), 
     week_8 = rownames(deg8week.filt[deg8week.filt$log2FoldChange > 0,])
)

upset_upregulation <- upset(fromList(list.upregulation),
      nsets = length(list.upregulation),
      nintersects = 20,
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "Количество элементов",
      sets.x.label = "Элементов в наборе",
      text.scale = c(1.5, 1.5, 1, 1, 2, 1))


list.downregulation <- list(week_2 = rownames(deg2week.filt[deg2week.filt$log2FoldChange < 0,]), 
     week_4 = rownames(deg4week.filt[deg4week.filt$log2FoldChange < 0,]), 
     week_6 = rownames(deg6week.filt[deg6week.filt$log2FoldChange < 0,]), 
     week_8 = rownames(deg8week.filt[deg8week.filt$log2FoldChange < 0,])
)

upset_downregulation <- upset(fromList(list.downregulation),
      nsets = length(list.downregulation),
      nintersects = 20,
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "Количество элементов",
      sets.x.label = "Элементов в наборе",
      text.scale = c(1.5, 1.5, 1, 1, 2, 1))

upregulation <- Reduce(intersect, list.upregulation)
downregulation <- Reduce(intersect, list.downregulation)

# Get output paths
if (exists("snakemake")) {
  up_path <- snakemake@output[["up"]]
  down_path <- snakemake@output[["down"]]
  pca_plot_path <- snakemake@output[["pca_plot"]]
  pca_variance_path <- snakemake@output[["pca_variance"]]
} else {
  # Fallback if running outside Snakemake
  output_dir <- file.path(work_dir, "dif_expression_results")
  up_path <- file.path(output_dir, "upregulation_star_htseq.txt")
  down_path <- file.path(output_dir, "downregulation_star_htseq.txt")
  pca_plot_path <- file.path(output_dir, "PCA_plot_star_htseq.png")
  pca_variance_path <- file.path(output_dir, "PCA_variance_star_htseq.png")
}

# Create output directory
dir.create(dirname(up_path), recursive = TRUE, showWarnings = FALSE)

# Save gene lists
write(upregulation, file = up_path)
write(downregulation, file = down_path)

# Save PCA plots
ggsave(filename = pca_plot_path, plot = PCA_vst_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = pca_variance_path, plot = prop_vst_plot, width = 10, height = 8, dpi = 300)



