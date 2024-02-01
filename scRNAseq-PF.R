###################################
# Single-cell RNA-seq Analysis    #
# Pulmonary Fibrosis  GSE135893   #
#                                 #
# M. Volkan Atalay                #
# Rengul Cetin Atalay             #
# Robert Hamanaka                 # 
# Gokhan Mutlu                    #
#                                 #
# January 2024                    #
###################################



## ----LIBRARIES----------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)


## ----READ COUNT MATRIX--------------------------------------------------------------------------------
# Read `matrix.mtx`
counts <- readMM("../Data/Source/GSE135893_matrix.mtx")
genes <- read_tsv("../Data/Source/GSE135893_genes.tsv", col_names = FALSE, show_col_types = FALSE)$X1
cell_ids <- read_tsv("../Data/Source/GSE135893_barcodes.tsv", col_names = FALSE, show_col_types = FALSE)$X1


## ----GENE NAMES and CELL IDs--------------------------------------------------------------------------
# gene names in the rows and cell identifiers in the columns
rownames(counts) <- genes
colnames(counts) <- cell_ids


## ----SUMMARY COUNTS-----------------------------------------------------------------------------------
# summary of counts for cells
counts_per_cell <- Matrix::colSums(counts)
head(counts_per_cell)


## ----DISTRIBUTION COUNTS PER CELL---------------------------------------------------------------------
# distribution of the counts per cell
cpc_hist.plot <- hist(log10(counts_per_cell+1),main='counts per cell',col='wheat', breaks = 100)


## -----------------------------------------------------------------------------------------------------
# summary of counts for genes
counts_per_gene <- Matrix::rowSums(counts)
head(counts_per_gene)


## ----DISTRIBUTION COUNTS PER GENE---------------------------------------------------------------------
# distribution of the counts per gene
cpg_hist.plot <- hist(log10(counts_per_gene+1),main='counts per gene',col='wheat', breaks = 200)


## ----GENES IN EACH CELL and CELLS PER GENE------------------------------------------------------------
# number of genes with a count more than 0 for a cell; count gene only if it has non-zero reads mapped.
genes_per_cell <- Matrix::colSums(counts>0) 

# number of cells with count more than 0 for a gene; only count cells where the gene is expressed
cells_per_gene <- Matrix::rowSums(counts>0) 


## ----DISTRIBUTION GENES PER CELL----------------------------------------------------------------------
gpcell_hist.plot <- hist(log10(genes_per_cell+1),main='genes per cell',col='wheat', breaks = 100)


## ----DISTRIBUTION CELLS PER GENE----------------------------------------------------------------------
cpgene_hist.plot <- hist(log10(cells_per_gene+1),main='cell per gene',col='wheat', breaks = 100)


## ----CREATE ild_seurat--------------------------------------------------------------------------------
# create a Seurat object with a minimum criterion
ild_seurat <- CreateSeuratObject(counts= counts, min.features=100)


## ----INFORMATION ild_seurat---------------------------------------------------------------------------
ild_seurat


## ----INFORMATION ild_seurat METADATA------------------------------------------------------------------
# Explore ild_seurat object
head(ild_seurat@meta.data)


## ----NUMBER of GENES per UMI--------------------------------------------------------------------------
# Add the number of genes per UMI for each cell to metadata
ild_seurat$log10GenesPerUMI <- log10(ild_seurat$nFeature_RNA) / log10(ild_seurat$nCount_RNA)


## ----MITOCHONDRIAL GENES------------------------------------------------------------------------------
grep("^MT-", rownames(ild_seurat@assays$RNA$counts), value = TRUE)


## ----MITOCHONDRIAL RATIO------------------------------------------------------------------------------
ild_seurat$mitoRatio <- PercentageFeatureSet(object = ild_seurat, pattern = "^MT-")
ild_seurat$mitoRatio <- ild_seurat@meta.data$mitoRatio / 100

head(ild_seurat$mitoRatio)


## -----------------------------------------------------------------------------------------------------
levels(ild_seurat@active.ident) <- rep("source",43)


## -----------------------------------------------------------------------------------------------------
VlnPlot(ild_seurat, features=c("nCount_RNA","mitoRatio"), raster = FALSE, alpha = 0.01)


## -----------------------------------------------------------------------------------------------------
ild_seurat@meta.data    %>% 
  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  	geom_point() + 
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() 


## ----COPY METADATA and RENAME COLUMNS-----------------------------------------------------------------
# Create metadata dataframe
metadata <- ild_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(Library_ID = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)


## -----------------------------------------------------------------------------------------------------
df_description = read.csv(file = "../Data/Source/aba1972_table_s2.csv")
df_LibDiag <- select(df_description, Library_ID, Diagnosis)

ident <- metadata$Library_ID 
list_ident = c()

ind <- 1
for (el in ident) {
  list_ident[ind] <- df_LibDiag[df_LibDiag$Library_ID == el,]$Diagnosis
  ind <- ind+1
}

factor_ident <- factor(list_ident)
metadata$diagnosis <- factor_ident 



## -----------------------------------------------------------------------------------------------------
ident <- metadata$diagnosis
list_ident = c()
pf_samples <- c("IPF", "ILD", "NSIP", "cHP", "Sarcoidosis")

list_ident <- sapply(ident, function(ita) ifelse(ita %in% pf_samples,"PF","Control")) 

factor_ident <- factor(list_ident)
metadata$condition <- factor_ident


## -----------------------------------------------------------------------------------------------------
df_description = read.csv(file = "../Data/Source/aba1972_table_s2.csv")
df_LibDiag <- select(df_description, Library_ID, Sample_Name)

ident <- metadata$Library_ID 
list_ident = c()

ind <- 1
for (el in ident) {
  list_ident[ind] <- df_LibDiag[df_LibDiag$Library_ID == el,]$Sample_Name
  ind <- ind+1
}

factor_ident <- factor(list_ident)
metadata$Sample_Name <- factor_ident 



## -----------------------------------------------------------------------------------------------------
# Add metadata back to the Seurat object
ild_seurat@meta.data <- metadata


## -----------------------------------------------------------------------------------------------------
# Visualize the number of cell counts per diagnosis
metadata %>% 
  	ggplot(aes(x=diagnosis, fill=diagnosis)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")


## -----------------------------------------------------------------------------------------------------
# Visualize the number of cell counts per condition
metadata %>% 
  	ggplot(aes(x=condition, fill=condition)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")


## -----------------------------------------------------------------------------------------------------
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  	ggplot(aes(color=condition, x=nUMI, fill= condition)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)



## -----------------------------------------------------------------------------------------------------
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=condition, x=nGene, fill= condition)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 250)


## -----------------------------------------------------------------------------------------------------
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = condition, fill=condition)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)


## -----------------------------------------------------------------------------------------------------
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=condition, x=mitoRatio, fill=condition)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.20)


## -----------------------------------------------------------------------------------------------------
# Visualize the correlation between genes detected and the number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio, group = condition)) +   
  	geom_point() + 
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  	stat_smooth(method = lm, aes(group = 1)) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) + 
  	facet_wrap(~condition)


## -----------------------------------------------------------------------------------------------------
# Filter out low-quality cells using selected thresholds - these will change with the experiment
filtered_seurat <- subset(x = ild_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))


## -----------------------------------------------------------------------------------------------------
# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


## -----------------------------------------------------------------------------------------------------
filtered_seurat


## -----------------------------------------------------------------------------------------------------
# Visualize the correlation between genes detected and the number of UMIs and determine whether the strong presence of cells with low numbers of genes/UMIs
filtered_seurat@meta.data %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio, group = condition)) + 
  	geom_point() + 
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~condition)


## -----------------------------------------------------------------------------------------------------
# Create RData object to load at any time
save(filtered_seurat, file="../data/generated/filtered_seurat-20240122.RData")


