# scRNA-seq of single-cell suspensions generated from 20 PF and 10 nonfibrotic control lungs 

##### M. Volkan Atalay, Rengul Cetin Atalay, Robert Hamanaka, Gokhan Mutlu

##### January 2024



### Introduction

This is an analysis of the publicly available single-cell RNA-sequencing data from Habermann et al. Science Advances 2020, using R package `Seurat`. The single-cell RNA-sequencing dataset was formed by the single-cell suspensions generated from 20 pulmonary fibrosis (PF) and 10 nonfibrotic control lungs. The article and dataset can be accessed at: https://doi.org/10.1126/sciadv.aba1972 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893, respectively.

This is the first part of the analysis and it contains pre-processing the data and its pre-analysis. The steps include

- importing packages,
- reading the data and extracting its basic statistics,
- creating the Seurat object and adding information to metadata,
- adding more information to metadata,
- visualizing various ratios and distributions,
- filtering low-quality cells,
- saving the Seurat object to load at any time.

Remark this is a combined dataset of 30 samples (n=30); 10 control and 20 PF samples. There are 43 library identifiers which means that for some of the samples, there is more than one library. The samples are also classified in terms of six diagnosis types: chronic hypersensitivity pneumonitis (cHP; n = 3), nonspecific interstitial pneumonia (NSIP; n = 2), sarcoidosis (n = 2), unclassifiable ILD (ILD; n = 1), and nonfibrotic controls (n = 10). The samples and cells can also be labeled as Control and PF. The statistics show that the data is very heterogeneous.



### Load packages

```
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(Matrix)

```

### Read in Data and Basic Statistics

The data can be downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893. You can download three files which are at the bottom of the webpage: GSE135893_barcodes.tsv.gz, GSE135893_genes.tsv.gz, and GSE135893_matrix.mtx.gz. These files contain a set of matrices that correspond to the cell identifiers (cellular barcodes), gene names, and expression values for each cell.

Expression values can be read by `readMM()` while `read_tsv()` can be used to access the other two files.

```
# Read in `matrix.mtx`
counts <- readMM("../Data/Source/GSE135893_matrix.mtx")
genes <- read_tsv("../Data/Source/GSE135893_genes.tsv", col_names = FALSE, show_col_types = FALSE)$X1
cell_ids <- read_tsv("../Data/Source/GSE135893_barcodes.tsv", col_names = FALSE, show_col_types = FALSE)$X1
```

The expression matrix should have gene names in the rows and cell identifiers in the columns. 

```
rownames(counts) <- genes
colnames(counts) <- cell_ids
```

Let’s have the summary of counts for cells. `colSums()` is a function that works on each column in a matrix and returns the sum of the elements in the corresponding column as a vector element.

```
# Summary of counts for cells
counts_per_cell <- Matrix::colSums(counts)
head(counts_per_cell)
## F01157_AAACCTGAGCATCATC F01157_AAACCTGGTGAACCTT F01157_AAACCTGTCATCGCTC 
##                   18200                    1072                     600 
## F01157_AAACCTGTCATGCTCC F01157_AAACCTGTCCACGTTC F01157_AAACGGGAGCTGCCCA 
##                   13747                     563                     548
```

A general view of the distribution of the counts per cell.

```
# distribution of the counts per cell
cpc_hist.plot <- hist(log10(counts_per_cell+1),main='counts per cell',col='wheat', breaks = 100)
```

![1-counts-per-cell](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/d639bf54-e572-4180-8bc4-4d46c03284e6)



Similarly, let’s have the summary of counts for genes. `rowSums()` is a function that works on each row in a matrix and returns the sum of the elements in the corresponding row as a vector element.

```
# Summary of counts for genes
counts_per_gene <- Matrix::rowSums(counts)
head(counts_per_gene)
##  RP11-34P13.3       FAM138A         OR4F5  RP11-34P13.7  RP11-34P13.8 
##            11             1             0           102             4 
## RP11-34P13.14 
##            11
```

A general view of the distribution of the counts per gene.

```
# distribution of the counts per gene
cpg_hist.plot <- hist(log10(counts_per_gene+1),main='counts per gene',col='wheat', breaks = 200)
```

![2-counts-per-gene](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/77c0203d-80de-4f80-9b41-5f02648ca25e)


Let’s also check the number of genes in each cell and in how many cells a particular gene is found.

```
# number of genes with a count more than 0 for a cell; count gene only if it has non-zero reads mapped.
genes_per_cell <- Matrix::colSums(counts>0) 

# number of cells with count more than 0 for a gene; only count cells where the gene is expressed
cells_per_gene <- Matrix::rowSums(counts>0) 
gpcell_hist.plot <- hist(log10(genes_per_cell+1),main='genes per cell',col='wheat', breaks = 100)
```

![3-genes-per-cell](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/07a58e26-774b-4e48-b67f-c7e3ae72425f)


```
cpgene_hist.plot <- hist(log10(cells_per_gene+1),main='cell per gene',col='wheat', breaks = 100)
```


![4-cells-per-gene](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/f3a84fce-f292-446d-84eb-41b49971c235)



### Create Seurat object and information to metadata

Create a Seurat object with a minimum criterion: if the number of genes in a cell is less than 100, then discard the cell.

```
# create a Seurat object with a minimum criterion
ild_seurat <- CreateSeuratObject(counts= counts, min.features=100)
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
## Warning: Data is of class dgTMatrix. Coercing to dgCMatrix.
```

Information regarding the Seurat object

```
ild_seurat
## An object of class Seurat 
## 33694 features across 216528 samples within 1 assay 
## Active assay: RNA (33694 features, 0 variable features)
##  1 layer present: counts
# Explore ild_seurat object
head(ild_seurat@meta.data)
##                         orig.ident nCount_RNA nFeature_RNA
## F01157_AAACCTGAGCATCATC     F01157      18200         4259
## F01157_AAACCTGGTGAACCTT     F01157       1072          587
## F01157_AAACCTGTCATCGCTC     F01157        600          421
## F01157_AAACCTGTCATGCTCC     F01157      13747         3492
## F01157_AAACCTGTCCACGTTC     F01157        563          396
## F01157_AAACGGGAGCTGCCCA     F01157        548          393
```

Seurat automatically generates very useful two metrics:

`nCount_RNA` is the total number of UMIs (or reads) 
`nFeature_RNA` is the number of observed genes

Let’s add the number of genes per UMI for each cell and the percentage of transcripts that map to mitochondrial genes to metadata.

```
# Add the number of genes per UMI for each cell to metadata
ild_seurat$log10GenesPerUMI <- log10(ild_seurat$nFeature_RNA) / log10(ild_seurat$nCount_RNA)
```

Large numbers of reads coming from mitochondria may fill in single cell dataset. Such cells are often sick cells undergoing apoptosis and they should be eliminated. Gene names can be accessed by `rownames` of the `@assays$RNA$counts` slot of the Seurat object and the mitochondrial genes can be identified by their names starting with “MT-”.

```
grep("^MT-", rownames(ild_seurat@assays$RNA$counts), value = TRUE)
##  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3" 
##  [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"
ild_seurat$mitoRatio <- PercentageFeatureSet(object = ild_seurat, pattern = "^MT-")
ild_seurat$mitoRatio <- ild_seurat@meta.data$mitoRatio / 100

head(ild_seurat$mitoRatio)
## F01157_AAACCTGAGCATCATC F01157_AAACCTGGTGAACCTT F01157_AAACCTGTCATCGCTC 
##              0.04697802              0.09608209              0.04500000 
## F01157_AAACCTGTCATGCTCC F01157_AAACCTGTCCACGTTC F01157_AAACGGGAGCTGCCCA 
##              0.04044519              0.03730018              0.04379562
```

There are 43 levels in `active.ident` (corresponding to library identifiers). Let’s make them a single label to observe all data (not splitted by 43 groups).

```
levels(ild_seurat@active.ident) <- rep("source",43)
```

Now, we can observe the number of counts and mitocondrial ratio for all of the dataset.

```
VlnPlot(ild_seurat, features=c("nCount_RNA","mitoRatio"), raster = FALSE, alpha = 0.01)
```

![5-nCount-mitoRatio](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/10c76e57-97b4-4866-9c5a-630d3bf19db1)





Visualize the correlation between the number of transcripts and the number of genes colored by the mitoRatio

```
ild_seurat@meta.data    %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
    geom_point() + 
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() 
```


![6-corr-transcripts-genes](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/7eecfd50-5e5f-419a-9482-7185d279f246)



### Add More Information to metadata

Let’s copy the metadata of the Seurat object to a separate dataframe and perform all operations (for the time being) on this dataframe.

```
# Create metadata dataframe
metadata <- ild_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(Library_ID = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
```

It would be useful to add a new entry to metadata as Diagnosis. Mark the cells according to Control, cHP, ILD, NSIP, IPF, and Sarcoidosis.

```
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
```

Similarly, add a new entry to metadata as condition. Mark the cells according to Control or PF (one of these cHP, ILD, NSIP, IPF, Sarcoidosis).

```
ident <- metadata$diagnosis
list_ident = c()
pf_samples <- c("IPF", "ILD", "NSIP", "cHP", "Sarcoidosis")

list_ident <- sapply(ident, function(ita) ifelse(ita %in% pf_samples,"PF","Control")) 

factor_ident <- factor(list_ident)
metadata$condition <- factor_ident
```

Finally, let’s add a new entry to metadata as Sample_Name.

```
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
```

Add metadata back to the Seurat object.

```
# Add metadata back to the Seurat object
ild_seurat@meta.data <- metadata
```



### Visualize various ratios and distributions

```
Visualize the number of cell counts per diagnosis.

# Visualize the number of cell counts per diagnosis
metadata %>% 
    ggplot(aes(x=diagnosis, fill=diagnosis)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")


    
```

![7-histogram-NCells-diagnosis](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/7bbfc12d-214e-4840-b4ad-dbe7a35a676d)


Visualize the number of cell counts per condition.

```
# Visualize the number of cell counts per condition
metadata %>% 
    ggplot(aes(x=condition, fill=condition)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
```

![8-histogram-NCells-condition](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/dbb1cb72-ad80-4807-a299-331e08348d4d)


Visualize the number UMIs/transcripts per cell.

```
# Visualize the number UMIs/transcripts per cell
metadata %>% 
    ggplot(aes(color=condition, x=nUMI, fill= condition)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
```


![9-density-UMIs-per-cell](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/0d93fb62-630e-4ff3-87ee-47bbcae42f20)



Visualize the distribution of genes detected per cell via histogram.

```
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
    ggplot(aes(color=condition, x=nGene, fill= condition)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 250)
```

![10-density-genes-per-cell](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/524c8c42-d77c-4b45-a2d6-808ad0d789e3)


Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score).

```
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = condition, fill=condition)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
```

![11-density-complexity](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/a6472c5e-affd-44ef-a35b-1d72754b104c)


Visualize the distribution of mitochondrial gene expression detected per cell.

```
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
    ggplot(aes(color=condition, x=mitoRatio, fill=condition)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.20)
```

![12-density-miyoRatio-per-cell](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/8de437f7-2319-4882-9f9f-d0691ecab7c9)




Visualize the correlation between genes detected and the number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs.

```
# Visualize the correlation between genes detected and the number of UMIs and determine whether the strong presence of cells with low numbers of genes/UMIs
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
```

![13-correlation-genes-UMIs-before](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/5ba9d2f1-862d-4ce7-8705-2cdb6293f050)



### Filter low-quality cells

Filter out cells with less than 250 genes and more than 20% mitoRatio

Filter out low-quality cells using selected thresholds.

```
# Filter out low-quality cells using selected thresholds - these will change with the experiment
filtered_seurat <- subset(x = ild_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
```

Some more filtering.

```
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
filtered_seurat
## An object of class Seurat 
## 25153 features across 190689 samples within 1 assay 
## Active assay: RNA (25153 features, 0 variable features)
##  1 layer present: counts
```

Sanity check for the filtered data.

```
# Visualize the correlation between genes detected and the number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
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
```

![14-correlation-genes-UMIs-after](https://github.com/mvolkanatalay/scRNAseq-PulmonaryFibrosis/assets/134640971/ca79d6b7-cef8-4bc5-9c37-f6d346b443eb)




### Save Seurat object

Save the Seurat object to load at any time.

```
# Create RData object to load at any time
save(filtered_seurat, file="../data/generated/filtered_seurat-20240122.RData")
```



#### More Code

https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Setup

