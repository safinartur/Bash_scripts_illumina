
setwd("D:/Desktop/counts/star_output")
getwd()
install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'ensembldb', 'biomaRt', 'org.Hs.eg.db'))
install.packages(c('tidyverse', 'naniar', 'writexl'))
library('DESeq2')
library('naniar')
library('ensembldb')
library('org.Hs.eg.db')
library('biomaRt')
library('tidyverse')
library('writexl')

ff <- list.files(path = ".", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(sapply(counts.files, function(x) x[ ,2]))
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff)
ff <- gsub("counts/", "", ff)
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
counts <- counts[order(row.names(counts)), ]

write_tsv(rownames_to_column(counts, "gene.id"), "star_counts.tsv")

filtered_counts <- counts[apply(counts, 1, mean) >= 10, ]


design <- tibble("Sample" = colnames(counts),
                 "group" = if_else(colnames(counts) %in% c('./SRR17948866', './SRR17948874', './SRR17948859'), "T2D", "Active"))
#___________________________________

getwd()
setwd("D:/Desktop/counts")
#Создаем таблицу
ff <- list.files(path = "/Desktop/counts", pattern = "a*.tsv$",
                 full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 1)
counts <- as.data.frame(round((sapply(counts.files, function(x) x[ , 4])), digits = 0))
ff <- gsub("/Desktop/counts/", "", ff)
ff <- gsub("[.]tsv", "", ff)
ff <- gsub("SRR179488", "", ff)
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
counts <- counts[order(row.names(counts)), ]

#Создаем вторую таблицу
sample_names <- c("58", "59", "66", "69", "70", "71", "73", "74")
group <- if_else(sample_names %in% c('59', '66', '72', '74'), "T2D", ifelse(sample_names %in% c('69', '70', '71', '73'), "Active", "Obese"))
color <- if_else(sample_names %in% c('59', '66', '72', '74'), 'red', ifelse(sample_names %in% c('67', '69', '70','71', '73'), 'blue', 'yellow'))

design <- tibble("Sample" = sample_names, "group" = group, "color" = color)
cool_names <- paste(design$group, colnames(counts), sep = "_")
colnames(counts) <- cool_names

# QC ----------------------------------------------------------------------

pdf('qc_counts.pdf')
counts %>%
  ggplot(aes(x = log10(1 + apply(., 1, mean)))) +
  geom_histogram(color="black", fill="lightblue") +
  ggtitle("Distribution of overall mean coverage") +
  theme_minimal() +
  scale_x_log10(name = "log10 coverage") +
  scale_y_continuous(name = "Frequency")

dev.off()

table(apply(counts, 1, mean) == 0)
table(apply(counts, 1, mean) >= 10)

pdf('distribution coverage.pdf')
mean_sample_cov_hist <- function(sample, sample_name, data) {
  ggplot(as.data.frame(data), aes(log(sample))) +
    geom_histogram(color = "black", fill = "lightblue") +
    ggtitle(paste("Distribution of coverage in", sample_name)) +
    scale_y_continuous(name = "Frequency") +
    scale_x_log10(name = 'log10 coverage')
}

for(i in 1:6) print(mean_sample_cov_hist(counts[, i], colnames(counts)[i], counts))
dev.off()


# filtering counts --------------------------------------------------------


filtered_counts <- counts[apply(counts, 1, mean) >= 10, ]
write_tsv(rownames_to_column(filtered_counts, 'gene.id'),
          'filtered_star_counts.tsv')


# correlations ------------------------------------------------------------


cor_logp <- cor(log(1 + counts), m = 'p')
filtered_cor_logp <- cor(log(1 + filtered_counts), m = 'p')


pdf('heatmap.pdf', height = 6, width = 6)
heatmap(1 - cor_logp, symm = T,
        distfun = function(x) {as.dist(x) })
heatmap(1 - filtered_cor_logp, symm = T,
        distfun = function(x) { as.dist(x) })
dev.off()


# PCA ---------------------------------------------------------------------


filtered_cpm <- sweep(filtered_counts, 2, apply(filtered_counts, 2, sum), '/') * 1e6
filtered_pca <- prcomp(t(filtered_cpm), scale = TRUE)

head(filtered_pca$x)

# Create a data frame for ggplot ------------------------------------------

pca_data <- as.data.frame(filtered_pca$x[, 1:2]) %>% rownames_to_column("Sample") %>%
  left_join(design)
barplot_data <- data.frame(Component = paste0('PC', 1:length(filtered_pca$sdev)),
                           Standard_Deviation = filtered_pca$sdev)

# which PC gives max variance ---------------------------------------------

pdf("clustering.pdf", height = 10, width = 10)
ggplot(barplot_data, aes(x = Component, y = Standard_Deviation)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Principal Component", y = "Standard deviation") +
  ggtitle("Variance per component in PCA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# Create custom axis labels -----------------------------------------------

filtered_var <- paste0('PC', 1:length(filtered_pca$sdev), ' ', round(filtered_pca$sdev / sum(filtered_pca$sdev) * 100, 1), '%')

pca_plt <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(group))) +
  geom_point(aes(label = rownames(pca_data)), pch = 19, size = 3) +
  labs(x = filtered_var[1], y = filtered_var[2]) +
  theme_minimal() +
  ggtitle("PCA plot") +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5) +
  labs(colour = "Condition") +
  scale_color_manual(values = c("red3", "cyan3", "purple3"))
ggsave("pca_plot.pdf", plot = pca_plt, device = "pdf", width = 10, height = 6)
#print(pca_plt)
dev.off()

# k-means -----------------------------------------------------------------
library('factoextra')

cpm_t <- t(sweep(filtered_counts, 2, apply(filtered_counts, 2, sum), '/') * 1e6)

kmeans <- kmeans(cpm_t, 2, nstart = 50)
kmeans$cluster
km_plt <- fviz_cluster(kmeans, data = cpm_t) + theme_minimal() + ggtitle("K-means clustering")
print(km_plt)

gridExtra::grid.arrange(pca_plt, km_plt, ncol = 2) -> plot_k_means

ggsave("PCA_plot+K-meands_clustering.pdf", plot = plot_k_means, device = "pdf", width = 8, height = 6)

dev.off()

# adding gene names -------------------------------------------------------

library(BiocManager)

library(BiocFileCache)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "ensembl_peptide_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(filtered_counts),
                  mart = ensembl )

add_genesym <- function(deseq_res_obj, gene_annot = genemap) {
  idx <- match(rownames(deseq_res_obj), gene_annot$ensembl_gene_id )
  deseq_res_obj$prot <- gene_annot$ensembl_peptide_id[ idx ]
  deseq_res_obj$hgnc_symbol <- gene_annot$hgnc_symbol[ idx ]
  return(deseq_res_obj)
}


# filtering and building deseq object -------------------------------------
dim(counts)
nrow(colData)

col(counts)
row(design)

design$group
counts
dds <- DESeqDataSetFromMatrix(counts, design, design = ~ group)

# calculating deseq -------------------------------------------------------

dds <- DESeq(dds)
resultsNames(dds)

# retrieving results ------------------------------------------------------

res <- results(dds) %>% as.data.frame() %>%
  add_genesym() %>% rownames_to_column("gene.id") %>%
  arrange(log2FoldChange, padj) %>% as_tibble()

# filtering significant ---------------------------------------------------

dge <- list()

dge$upreg <- dplyr::filter(res, log2FoldChange > 1.3 & pvalue < 0.05)
dge$downreg <- dplyr::filter(res, log2FoldChange < -1.3 & pvalue < 0.05)

dge$upreg <- dplyr::filter(res, log2FoldChange > 1.3 & padj  < 0.05)
dge$downreg <- dplyr::filter(res, log2FoldChange < -1.3 & padj < 0.05)


nrow(dge$upreg)

nrow(dge$downreg)



# Volcano plot ------------------------------------------------------------

#write.csv(x = signGenes, file = "./signGenes.csv", row.names = TRUE)


pval = 0.05
lfc = 1.3
res$signGenes = (abs(res$log2FoldChange) > lfc & -log10(res$padj) > -log10(pval)) 
pdf("Volcano_plot_after_DES_p=0,05.pdf")
res <- res[is.na((res$signGenes))==FALSE,]
ggplot(res, aes(x=log2FoldChange,y=-log10(padj))) +
  geom_jitter(aes(colour = signGenes), size =3) +
  geom_hline(yintercept = -log10(pval), color = "green4", size = 1) +
  geom_vline(xintercept = c(-lfc, lfc), linetype='dotted', color = "blue", size = 1) +
  ggtitle(sprintf("%s", "Volcano Plot")) +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 8)) + 
  xlim(-5,5) +
  
  theme_bw() + scale_fill_grey() 
dev.off()

# cleaning ----------------------------------------------------------------

rm(counts, counts.files, ff, dds)


library("clusterProfiler")


nrow(dge$upreg)

nrow(dge$downreg)


upregs_a <- dge$upreg
head(upregs_a)


gene <- upregs_a %>% mutate(rank = rank(log2FoldChange,  ties.method = "random")) %>%
  arrange(desc(rank)) %>% dplyr::pull(rank)

names(gene) <- upregs_a %>% arrange(desc(log2FoldChange)) %>% dplyr::pull(gene.id)

names(gene) = mapIds(org.Hs.eg.db,
                     keys=names(gene), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
head(gene)
enrich_full <- gseKEGG(geneList = gene,
                       organism = 'hsa', 
                       by = "fgsea", 
                       pvalueCutoff = 1.0)

print("Number of all pathways is the first number")
print(dim(enrich_full@result))

enrich <- gseKEGG(geneList = gene,
                  organism = 'hsa', 
                  by = "fgsea")
library(enrichplot)
print("Number of significant pathways is the first number")
print(dim(enrich@result))
head(enrich@result)

net_data <- pairwise_termsim(enrich_full, method = "JC")
net_data

emapplot(net_data, showCategory = 10) + 
  ggtitle("Network map for KEGG pathways of upregulated genes in group A") -> p1
ggsave("star_final_network.pdf", plot = p1, device = "pdf", width = 8, height = 6)



ridgeplot(enrich_full, fill = "pvalue", showCategory = 10) + 
  ggtitle("Ridge plot for KEGG pathways of upregulated genes in group A") -> p2
ggsave("star_ridge.pdf", plot = p2, device = "pdf", width = 8, height = 6)


enrichplot::dotplot(enrich_full, color = "pvalue", showCategory = 10) + 
  ggtitle("GSEA dotplot for KEGG pathways of upregulated genes in group A") -> p3
ggsave("star_dotplot.pdf", plot = p3, device = "pdf", width = 8, height = 6)


sym_enrich <- setReadable(enrich_full, 'org.Hs.eg.db', 'ENTREZID') # convert entrez to symbol

cnetplot(sym_enrich, categorySize="pvalue", foldChange = gene, showCategory = 7, layout = 'kk') +
  ggtitle("KEGG pathway-gene network of the most significant pathways for upregulated genes in group A") -> p4
ggsave("star_cnetplot.pdf", plot = p4, device = "pdf", width = 8, height = 6)


enrich_Wiki <- gseWP(geneList = gene, 
                     organism = "Homo sapiens", 
                     by = "fgsea",
                     pvalueCutoff = 1.0)


print("Number of significant pathways is the first number")
print(dim(enrich_Wiki@result))
head(enrich_Wiki@result)


net_data <- pairwise_termsim(enrich_Wiki, method = "JC")

emapplot(net_data, showCategory = 10) + 
  ggtitle("Network map for WikiPathways of upregulated genes in group A") -> p5
ggsave("wikipathways_network.pdf", plot = p5, device = "pdf", width = 8, height = 6)


ridgeplot(enrich_Wiki, showCategory = 10, fill = "pvalue") + 
  ggtitle("Ridge plot for WikiPathways of upregulated genes in group A") -> p6
ggsave("wikipathways2_ridge.pdf", plot = p6, device = "pdf", width = 8, height = 6)



enrichplot::dotplot(enrich_Wiki, showCategory = 10, color = "pvalue") + 
  ggtitle("GSEA dotplot for WikiPathways of upregulated genes in group A") -> p7
ggsave("wikipathways_dotplot.pdf", plot = p7, device = "pdf", width = 8, height = 6)


#network for genes
sym_enrich <- setReadable(enrich_Wiki, 'org.Hs.eg.db', 'ENTREZID') # convert entrez to symbol
cnetplot(sym_enrich, categorySize="pvalue", foldChange = gene, showCategory = 7, layout = 'kk') +
  ggtitle("WikiPathways pathway-gene network of the most significant pathways for upregulated genes in group A") -> p8
ggsave("wikipathways_cnetplot.pdf", plot = p8, device = "pdf", width = 8, height = 6)

