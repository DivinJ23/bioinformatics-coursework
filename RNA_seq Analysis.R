load("BMS3025_RNAseq.RData")
save.image(file = "BMS3025_RNAseq.RData")

#Installing the packages and bioconductor packages
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")

#Loading packages
library(R.utils)
library(readr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)

#Download data
download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz', destfile='GSE147507_RawReadCounts_Human.tsv.gz')

#Unzip data 
gunzip('GSE147507_RawReadCounts_Human.tsv.gz', remove=FALSE)

#Import data into R
dat <- read.delim('GSE147507_RawReadCounts_Human.tsv', sep = '\t')

#Formatting and labeling Column 1
dat <- as.data.frame(dat)
names(dat)[names(dat) == "X"] <- "Gene"
colnames(dat)

#Extracting SARS-CoV-2, RSV and IAV data from basic A549 data set
dat2 <- as.matrix(dat[,14:27])
rownames(dat2) <- dat$Gene

#Definition of metadata

sample_table <- data.frame(sample_id=colnames(dat2), 
                           condition=c('Mock_RSV', 'Mock_RSV', 'RSV', 'RSV', 
                                       'Mock_IAV', 'Mock_IAV', 'IAV', 'IAV', 
                                       'Mock_Covid', 'Mock_Covid', 'Mock_Covid', 'Covid', 'Covid', 'Covid'))
View(sample_table)

#Converting data into DESeqDataSet
dds <- DESeqDataSetFromMatrix(dat2, colData = sample_table, design=~condition)
dds

#Analysing the data with DESeq2
dds2 <- DESeq(dds)
dds2

#rlog step and PCA plot
rld <- rlog(dds2, blind=TRUE)
plotPCA(rld, intgroup='condition')

# Measuring sample distances using the rld DataSet
sampleDists <- dist(t(assay(rld))) 

#Converting to matrix and labelling data
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition)
colnames(sampleDistMatrix) <- paste(rld$condition)

#Generating heatmap of the matrix
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

save.image(file = "BMS3025_RNAseq.RData")

# Creating a reduced dataset containing only RSV and SARS-CoV-2 infections 
dat3 <- as.matrix(dat[,c(14:17, 22:27)]) 
rownames(dat3) <- dat$Gene 

# Define sample conditions for metadata
sample_table <- data.frame(sample_id=colnames(dat3), 
                           condition=c('Mock_RSV','Mock_RSV','RSV','RSV',
                                       'Mock_Covid','Mock_Covid','Mock_Covid',
                                       'Covid','Covid','Covid'))

# Constructing a DESeq2 dataset from the count matrix and metadata
dds3 <- DESeqDataSetFromMatrix(dat3, colData=sample_table, design=~condition)
dds3 <- DESeq(dds3)

# Results for RSV vs mock control
rsv_results <- results(dds3, contrast=c('condition','RSV','Mock_RSV'))

# Results for SARS-CoV-2 vs mock control
covid_results <- results(dds3, contrast=c('condition','Covid','Mock_Covid'))

# Defining significance thresholds 
cols <- ifelse(rsv_results$padj <= 0.05 & abs(rsv_results$log2FoldChange) >= 1, "red", "black")

# Plot volcano plot for RSV
plot(rsv_results$log2FoldChange, -log10(rsv_results$padj),
     col=cols, pch=16,
     xlab="Log2 Fold Change", ylab="-Log10 adjusted p-value",
     main="RSV vs Mock Volcano Plot")
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1,1), col="blue", lty=2)

# Interferon Heatmap: Transforming count data to stabilise variance
rld2 <- rlog(dds3, blind=TRUE)

# Extracting variance-stabilised expression for interferon-related genes
ifn_genes <- c("IFNA1","IFNA2","IFNB1","IFNL1","IFNL2","IFNL3") # example list
ifn_mat <- assay(rld2)[ifn_genes, ]

# Generating clustered heatmap
pheatmap(ifn_mat, 
         cluster_rows=TRUE, cluster_cols=TRUE, 
         show_rownames=TRUE, show_colnames=TRUE,
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
         main="Interferon Gene Expression: RSV vs SARS-CoV-2")



