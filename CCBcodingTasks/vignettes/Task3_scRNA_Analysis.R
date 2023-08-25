## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----load package-------------------------------------------------------------
library(CCBcodingTasks)

## ----install packages, include = FALSE----------------------------------------
#BiocManager::install('scRNAseq')
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("OSCA.intro")
#BiocManager::install("scater")
#BiocManager::install("scran")
#BiocManager::install("bluster", force=T)
#BiocManager::install("metapod")
#BiocManager::install("devtools")
#library(devtools)
#install_github("LTLA/metapod")
#install_github("MarioniLab/scran")

## ----load data----------------------------------------------------------------
# Load data set
suppressMessages(library(scRNAseq))
sce.zeisel <- ZeiselBrainData()

# Following more detailed loading instructions from the OSCA book on this specific dataset (http://bioconductor.org/books/3.17/OSCA.workflows/zeisel-mouse-brain-strt-seq.html#data-loading-1), we need to merge some redundant rows. 
suppressMessages(library(scater))
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))

# Take a look at our data summary
# We have ~19K genes across 3005 cells.
sce.zeisel 

## ----rowData------------------------------------------------------------------
head(rowData(sce.zeisel))
# Looks like we only have one column showing the featureType of the genes in this dataset.

## ----colData------------------------------------------------------------------
# List the different attributes of this experiment
names(colData(sce.zeisel))


## ----assay--------------------------------------------------------------------
dim(assay(sce.zeisel))

## ----QC-----------------------------------------------------------------------
# Quality control (using mitochondrial genes).
library(scater)

# Check which genes are mitochondrial
is.mito <- grepl("mito", rowData(sce.zeisel)$featureType)

# For each cell, calculate QC metrics and also calculates stats specifically for the subset of mt genes too
qcstats <- perCellQCMetrics(sce.zeisel, subsets=list(Mito=is.mito))

head(qcstats)

# The sum column contains the total count for each cell and the detected column contains the number of detected genes. The subsets_Mito_percent column contains the percentage of reads mapped to mitochondrial transcripts. Finally, the altexps_ERCC_percent column contains the percentage of reads mapped to ERCC transcripts.


# View distributions of our total counts, the ones that actually passed threshold, and our mito percentages
summary(qcstats$sum)
summary(qcstats$detected)
summary(qcstats$subsets_Mito_percent)

# Identify the low quality cells based on our previous QC statistics. Filter based on Mito percent and deviations from expected ERCC spike-in controls expression levels.
filtered <- quickPerCellQC(qcstats, percent_subsets=c("altexps_ERCC_percent", 
    "subsets_Mito_percent"))

# Retain only passing cells
sce.unfiltered <- sce.zeisel
sce.pass <- sce.zeisel[,!filtered$discard]
dim(sce.pass) # 2816 passing cells

# We can also pull out the reasons why certain genes did not pass 
reasons <- perCellQCFilters(qcstats, 
                            sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
summary(reasons$discard)

## ----QC mito plots, fig.height = 4, fig.width = 7-----------------------------
# First get our unfiltered dataset to compare to and add on filtering stats.
colData(sce.unfiltered) <- cbind(colData(sce.unfiltered), qcstats)
sce.unfiltered$discard <- filtered$discard

# View mitochondrial percentage.
plotColData(sce.unfiltered, y="subsets_Mito_percent",
        colour_by="discard") + ggtitle("Mito percent")

## ----QC mito plots2, fig.height = 4, fig.width = 7----------------------------

# plot across tissue types now
plotColData(sce.unfiltered, x="altexps_ERCC_percent", y="subsets_Mito_percent",
        colour_by="discard", other_fields="tissue") + facet_wrap(~tissue)


## ----lib, fig.height = 4, fig.width = 7---------------------------------------
# Library size across cells before normalization
tot <- as.data.frame(cbind(colnames(assay(sce.unfiltered)), as.integer(colSums(assay(sce.unfiltered)))))

tot <- tot[1:100,]

# Plot library for first 100 cells
p1 <- ggplot(data=tot, aes(x=V1, y=as.numeric(V2))) +
  geom_point() +theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Cell") + ylab("Library Size") 

print(p1)

## ----lib scale, fig.height = 4, fig.width = 7---------------------------------
library(scater)
lib.sf.zeisel <- librarySizeFactors(sce.pass)
summary(lib.sf.zeisel)

hist(log10(lib.sf.zeisel), xlab="Log10[Size factor]", col='grey80', main= "Histogram of library size factors") 

## ----noramlize, fig.height = 4, fig.width = 7---------------------------------
library(scran)
set.seed(100)

# Cluster similar cells based on their expression profiles, using either log-expression values or ranks.
clusters <- quickCluster(sce.pass)

# Get scaled factors and log the normalized counts based on our clustering prior
sce.zeisel <- computeSumFactors(sce.pass, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)
summary(sizeFactors(sce.zeisel))

plot(librarySizeFactors(sce.zeisel), sizeFactors(sce.zeisel), pch=16,
    xlab="Library size factors", ylab="Deconvolution factors",
    log="xy")


## ----HVGs, fig.height = 4, fig.width = 5--------------------------------------
# Model the variance of the log-expression profiles for each gene. Separate technical and biological components based on a mean-variance trend.
dec <- modelGeneVar(sce.zeisel)

# Grab the top 10% of genes with the largest biological components based on this mean-variance relationship
hvg <- getTopHVGs(dec, prop=0.1)
head(hvg)
length(hvg) #946 genes chosen

# Plot
dec.df <- dec
dec.df$Label <- "Other genes"
dec.df$Label[rownames(dec.df) %in% hvg] <- "HGVs"

dec.df <- as.data.frame(dec.df)

suppressMessages(library(dplyr))

# Plot
dec.df %>% ggplot(aes(bio, fill=Label)) +
    geom_histogram(aes(y=..count../sum(..count..)),
                 alpha=0.8,position='identity',binwidth=0.5) + theme_minimal()+
  xlab("Biological component of the variance") + ylab("Proportion of genes") 


## ---- fig.height = 4, fig.width = 7-------------------------------------------
# Performing PCA only on the chosen HVGs.
library(scater)
set.seed(1234)

# Run a PCA and specify the genes to use for downstream functions via an extra argument subset.row=
sce.pca <- runPCA(sce.zeisel, ncomponents=25, subset_row=hvg)
reducedDimNames(sce.pca)
ncol(reducedDim(sce.pca, "PCA")) # reduced down to 25 PCs

# Plot
p.pca <- plotReducedDim(sce.pca, dimred="PCA", colour_by="level1class")
print(p.pca)


## ----clust, fig.height = 4, fig.width = 7-------------------------------------
# Use the reduced dimensionality from our previous step and do clustering of cells.
library(bluster)
set.seed(1234)
sce.zeisel <- runTSNE(sce.pca, dimred="PCA")
ncol(reducedDim(sce.zeisel, "PCA"))

snn.gr <- buildSNNGraph(sce.zeisel, use.dimred="PCA")
colLabels(sce.zeisel) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

# View the number of clusters and how many cells are within each.
table(colLabels(sce.zeisel))

# Plot
tsne1 <- plotTSNE(sce.zeisel, colour_by="level1class", text_by="label")
tsne2 <-plotTSNE(sce.zeisel, colour_by="label", text_by="level1class")

tsne1
tsne2


## ----UMAP, fig.height = 4, fig.width = 7--------------------------------------
set.seed(1234)
sce <- runUMAP(sce.zeisel, dimred = 'PCA')
plotUMAP(sce, colour_by="label", text_by="level1class")

## ----heatmap, fig.height = 4, fig.width = 7-----------------------------------
# Let's look at all the upregulated genes using pairwise wilcoxon rank sum tests.
markers <- findMarkers(sce, test.type="wilcox", direction="up", lfc=1)

# Look genes that appear more differentially expressed in cluster 4 compared to others
marker.set <- markers[["4"]]
head(marker.set[,1:8], 10) 

# Cnp, Enpp2, Klhl2 seem to the top DE expressed in cluster 4. 
top.markers <- rownames(marker.set)[marker.set$Top <= 10]

# Create a heatmap 
plotHeatmap(sce.zeisel, features=top.markers, order_columns_by="label")


## ---- fig.height = 4, fig.width = 7-------------------------------------------
plotExpression(sce.zeisel, features=c("Cnp", "Enpp2",
    "Klhl2", "Mag"), x="label", colour_by="label")

## -----------------------------------------------------------------------------
sessionInfo()

