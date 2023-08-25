## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 100), tidy = TRUE
)


## ----load package-------------------------------------------------------------
library(CCBcodingTasks)

## ----loadPackages, message=FALSE, warning=FALSE-------------------------------
#BiocManager::install("airway")
#BiocManager::install(SummarizedExperiment)

library(SummarizedExperiment)
library(airway)
library(dplyr)
library(reshape2)
library(tibble)

# load data and view it
data(airway)
se <- airway

# RangedSummarizedExperiment object of read counts in genes for an RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone

## ----View Data----------------------------------------------------------------
# View experiment and sample information.
colData(se)

# View the gene information.
rowData(se) # 64102 gene names


## ----Filter by Read Count 8---------------------------------------------------
# Access the read count matrix of the experiment - For the airway dataset there is only one assay.
Count <- assay(se)
head(Count)

# Identify rows where read count is greater than or equal to 8
row_indices <- apply(Count, 1, function(row) all(row >= 8))

# How many genes remain after our filtering?
length(which(row_indices==T)) #14644

# Filter matrix based on identified rows
filtered_mat <- Count[row_indices,]

# Check the new matrix to make sure it only includes 14644 genes that passed.
dim(filtered_mat) 
head(filtered_mat)


## -----------------------------------------------------------------------------
# Create a function to divide read counts by the library size (sum of each column) multiplied by 10^6 to obtain counts per million reads mapped

scaleCounts <- function(mat) {
  col_sums <- colSums(mat)
  normalized_mat <- (mat / col_sums) * 10^6
  return(normalized_mat)
}

# Use the function on our filtered matrix above
normalized_mat <- scaleCounts(filtered_mat)

# Check that the read counts were changed correctly
head(normalized_mat)

## -----------------------------------------------------------------------------
# Grab the first 100 rows (genes)
mat100 <- normalized_mat[1:100,]

# Check the dimensions
dim(mat100)
head(mat100)


## -----------------------------------------------------------------------------
# Update the our SummarizedExperiment object to include genes in the subsetted assay matrix

# Find row indices in original dataset that are not present in the filtered
non_matching_indices <- which(!(rownames(se) %in% rownames(mat100)))

# Filter out the rows 
new_se <- se[-non_matching_indices, ]
new_se #top 100 genes that are normalized by library size.

# Now replace with our noramlized matrix.
assay(new_se) <- mat100
head(assay(new_se))

## -----------------------------------------------------------------------------
# Split SummarizedExperiment by treatment groups
trt <- new_se[, new_se$dex == "trt"] # treated
untrt <- new_se[, new_se$dex == "untrt"] # untreated

## ----t.tests------------------------------------------------------------------
# First, make a new column in rowData to store your P-values.
rowData(new_se)$pval <- NA

# loop through each gene, and update the pval column with the results of the t.test comparing the two groups.
for(i in 1:nrow(assay(new_se))){
  
  # loop through each gene in our SummarizedExperiment object. Extract current gene.
  c.gene <- rownames(new_se[i]) 
  
  # For the current gene, perform a t.test between the sample groups, adjusting for multiple correction.
  result <- t.test(assay(trt[c.gene,]), assay(untrt[c.gene,]))
  
  # Add the p.value to our rowData in our SummarizedExperiment object.
  rowData(new_se)$pval[i] <- as.numeric(result$p.value) 

}

# View your updated row information with the P-values for each gene.
head(rowData(new_se))


## -----------------------------------------------------------------------------
# Grab the top 5 signficant p values
Topvals <- sort(rowData(new_se)$pval, decreasing=F)[1:5]
topidx <- which(rowData(new_se)$pval %in% Topvals)

# Subset to just the top5 genes
top5 <- new_se[topidx]

# Extract treatment groups for just the top 5 genes
trt5 <- top5[, top5$dex == "trt"] # treated
untrt5 <- top5[, top5$dex == "untrt"] # untreated
# For each gene make a boxplot comparing the two groups

# Convert to dataframe for ease of plotting
trt5 <- as.data.frame(assay(trt5))
trt5$Condition <- "Treated"
untrt5 <- as.data.frame(assay(untrt5))
untrt5$Condition <- "Untreated"

# Function to convert count data into a combined dataframe in long format for plotting
edit_df <- function(df) {
  
  df %>% 
  rownames_to_column(var = "GeneID") %>% 
  melt(id.vars = c("Condition", "GeneID"), variable.name = "Sample", value.name = "ReadCount") 
  
}

# Combine the data frames into a list and apply the function. 
combo.list <- list(trt5, untrt5)
combo.list <- lapply(combo.list, edit_df)

# Combined the edited dataframes together into one
combo5 <- data.frame(do.call(rbind, combo.list))
head(combo5)




## ---- fig1, fig.height = 4, fig.width = 7-------------------------------------
library(ggplot2)

# Reorder so untreated is shown first
combo5$Condition <- factor(combo5$Condition, levels = c("Untreated", "Treated"))
combo5 <- combo5[order(combo5$Condition), ]

# Basic box plot of scaled read counts
p1 <-ggplot(combo5, aes(x= Condition, y= ReadCount, fill=GeneID)) + 
  geom_boxplot(alpha=0.85) + 
  theme_light() +
  ylab("Scaled Read Count (per million reads)") +
  theme(legend.title = element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14, vjust=-0.5), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14), axis.text.y=element_text(size=14))

print(p1)

## -----------------------------------------------------------------------------
# Make our assay read counts into a data frame
new_se.df <- as.data.frame(assay(new_se))

# Grab sample IDs and their indices belonging to each group.
treatSamps <- colnames(trt5)[1:4]
t.idx <- which(colnames(new_se.df) %in% treatSamps)

untreatSamps <- colnames(untrt5)[1:4]
u.idx <- which(colnames(new_se.df) %in% untreatSamps)

# Calculate log2 fold change for each gene in our 100 gene subset. Start by getting the mean read count of each gene per group. Then divide treatment mean with untreated mean to get a ratio of gene expression change. Finally, log2 that ratio.

# Create a dataframe that we will populate with for loop.
Log2FC.df <- data.frame(Gene=rowData(new_se), Log2FC=NA)


for (i in 1:nrow(new_se.df)){
  
  c.df <- new_se.df[i,] # For each gene at a time
  tmean <- rowMeans(c.df[,t.idx], na.rm=TRUE) # Caluclate treated mean
  untmean <- rowMeans(c.df[,u.idx], na.rm=TRUE) # Caluclate untreated mean
  Log2ratExp <- log2(tmean/untmean) # get Log2 ratio
  gene <- rownames(new_se.df[i,])
  Log2FC.df$Log2FC[i] <- Log2ratExp # Populate 100 gene dataframe
  
}


## ----fig2, fig.height = 4, fig.width = 5--------------------------------------
# Create a volcano plot
volcano.p <- ggplot(Log2FC.df, aes(x = Log2FC, y = -log10(pval))) +
  geom_point(aes(color = ifelse(rownames(Log2FC.df) %in% rownames(top5), "red", "black")), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Log Fold Change") +
  ylab("-log10(p-value)") +
  theme_light() +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size=14, vjust=-0.5), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14), axis.text.y=element_text(size=14))


# Display the volcano plot
print(volcano.p)

## -----------------------------------------------------------------------------
sessionInfo()

