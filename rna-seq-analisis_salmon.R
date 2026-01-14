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

DATA_DIR <- "/mnt/tank/scratch/ezaitseva/Neurodegeneration/EXTERNAL_MATURATION"

meta_data <- read.csv(file.path(DATA_DIR, "meta_data.txt"), sep = "\t", header = TRUE)

rownames(meta_data) <- meta_data$sampleid

meta_data$individual <- as.factor(meta_data$individual)
meta_data$diff_time <- as.factor(meta_data$diff_time)


Annotation_file <- "/mnt/tank/scratch/ezaitseva/Neurodegeneration/EXTERNAL_MATURATION/databases/STAR/GRCh38/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"

txdb <- txdbmaker::makeTxDbFromGFF(Annotation_file, format = "gtf", organism = "Homo sapiens")

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene) <- c("tx_id", "gene_id")

tx2gene$tx_id <- gsub("\\..*", "", tx2gene$tx_id)
tx2gene$gene_id <- gsub("\\..*", "", tx2gene$gene_id)

files <- file.path("/mnt/tank/scratch/ezaitseva/Neurodegeneration/EXTERNAL_MATURATION/salmon", meta_data$sampleid, "quant.sf")

names(files) <- meta_data$sampleid

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

ddsSalmon <- DESeqDataSetFromTximport(txi, colData = meta_data[-1], design = ~ individual + diff_time)

keep <- rowSums(counts(ddsSalmon) >= 10) > ncol(ddsSalmon)*0.9
ddsSalmon <- ddsSalmon[keep,]

vsd <- vst(ddsSalmon, blind = F)
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

ddsSalmon <- DESeq(ddsSalmon)

res2week <- results(ddsSalmon, contrast = c("diff_time", "2 week neuron", "NPC"))
res4week <- results(ddsSalmon, contrast = c("diff_time", "4 week neuron", "NPC"))
res6week <- results(ddsSalmon, contrast = c("diff_time", "6 week neuron", "NPC"))
res8week <- results(ddsSalmon, contrast = c("diff_time", "8 week neuron", "NPC"))

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

up_path   <- snakemake@output[["up"]]
down_path <- snakemake@output[["down"]]

dir.create(dirname(up_path), recursive = TRUE, showWarnings = FALSE)

write(upregulation, file = up_path)
write(downregulation, file = down_path)


