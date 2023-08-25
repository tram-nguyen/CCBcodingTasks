## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 100), tidy = TRUE
)

## ----setup--------------------------------------------------------------------
library(CCBcodingTasks)

## ----install, include=F-------------------------------------------------------
#BiocManager::install("GenomicRanges")
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm9") Includes full genome sequences for Mus musculus (Mouse) as provided by UCSC (mm9, Jul. 2007) and stored in Biostrings objects.
# BiocManager::install("BSgenome")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("GenomicFeatures")

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load packages
library(GenomicRanges)
library(plyranges)

# Load mouse
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome)

# Load in packages for the human genome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

## -----------------------------------------------------------------------------
# Load the Mmusculus mm9 genome
mm9_genome <- BSgenome.Mmusculus.UCSC.mm9

# Calculate chromosome lengths using GenomicRanges
chrom_lengths <- seqlengths(mm9_genome)

# Print the chromosome lengths
chrom_lengths


## ----subset Chr1--------------------------------------------------------------
# Load the TxDb object for hg38 human genome annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# View the object to ensure we loaded it correctly
txdb

# To check which chromosomes are active
head(seqlevels(txdb))

# Grab just Chromosome 1
seqlevels(txdb) <- "chr1"

# Extract transcript isoforms for each gene on Chr1
txn_by_genes <- transcriptsBy(txdb, "gene")

# View the isoforms for the first gene in the object. 
print(txn_by_genes[[1]]) 

## ----maxExon------------------------------------------------------------------
# We also want to grab the corresponding exons for all transcripts
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)

# Make an empty data frame with GeneID and the maximum Exon count for that gene per transcript.
MaxExon.df <- data.frame(GeneID = names(txn_by_genes), MaxExonCountPerTxn = NA)

# For loop to extract the maximum number of exons per transcript for each gene on Chr1.
for (gene in 1:length(txn_by_genes)){
  
    # Extract unique transcript names for each gene at a time
    txNames <- txn_by_genes[[gene]]$tx_name
    
    # Extract exons for just that gene
    c.df <- ex_by_tx[names(ex_by_tx) %in% txNames] # This show exons as separate GRanges objects for each transcript in a GRangesList.
    
    # Get the maximum exon count for transcripts 
    maxCount <- max(elementNROWS(c.df))
    
    # Populate your dataframe
    MaxExon.df[gene, 2] <- maxCount
}

head(MaxExon.df)

## ----old load, include=F------------------------------------------------------
#prom <- read.table("./data/hg19proms.bed", header = FALSE, sep = "\t", colClasses = c("character", "integer", "integer", "character","character","character")) # Promoters

#CpG <- read.table("./data/haowulab.org_software_makeCGI_model-based-cpg-islands-hg19.txt", header = TRUE, sep = "\t", colClasses = c("character", "integer", "integer", "integer","integer","integer", "numeric", "numeric")) # CpG

## ---- warning=F---------------------------------------------------------------
library(readr)
data(haowulabCpG, package = "CCBcodingTasks")
CpG <- as.data.frame(haowulabCpG)
head(CpG)

data(hg19proms.bed, package = "CCBcodingTasks")
prom <- as.data.frame(hg19proms.bed)
head(prom)

## -----------------------------------------------------------------------------
# Create GRanges objects
promGRanges <- GRanges(seqnames = prom$V1, ranges = IRanges(start = prom$V2, end = prom$V3), names = paste0("Prom", 1:nrow(prom)))

CpGRanges <- GRanges(seqnames = CpG$chr, ranges = IRanges(start = CpG$start, end = CpG$end), names = paste0("CpG", 1:nrow(CpG)))

## -----------------------------------------------------------------------------
Ch21.prom <- promGRanges %>% 
  filter(seqnames == 'chr21')

Ch21.CpG <- CpGRanges %>% 
  filter(seqnames == 'chr21')

## -----------------------------------------------------------------------------
overlaps <- findOverlaps(Ch21.prom, Ch21.CpG)

# Some promoters will overlap with CpG islands more than once. We want to just count these once.
length(unique(queryHits(overlaps)))

## -----------------------------------------------------------------------------
# Alternative, we can do:
num_overlapping_regions <- subsetByOverlaps(Ch21.prom, Ch21.CpG) 
length(num_overlapping_regions)

# Calculate percentage
# Get the total number of promoters on Chr21
tot <- length(Ch21.prom)
percentProm <- (length(num_overlapping_regions) / tot)*100
print(percentProm)

## -----------------------------------------------------------------------------
# For all promoters that overlap at least once, count how many times it overlaps with CpG island.
overlap_counts <- countOverlaps(Ch21.prom, Ch21.CpG) 
print(max(overlap_counts))

## -----------------------------------------------------------------------------
sessionInfo()

