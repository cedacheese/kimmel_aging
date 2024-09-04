install.packages("biomaRt")

BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("fastmap", force = TRUE)

library(biomaRt)
library(Seurat)
library(SeuratData)
library(Matrix)
library(patchwork)
library(dplyr)

tabula_senis <- readRDS("/Users/ceda_/Downloads/275d18e7-dbe1-47e6-93e9-642e895d0af7.rds")

unique(tabula_senis$age)
# [1] 18m 21m 24m 30m 1m  3m 
# Levels: 1m 3m 18m 21m 24m 30m

head(rownames(tabula_senis))

# Initialize biomaRt and select dataset
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Example for mouse genes

# Extract rownames from the Seurat object
ensembl_ids <- rownames(tabula_senis)

# Retrieve gene symbols for ENSEMBL IDs
gene_mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),  # Change to "hgnc_symbol" for human genes
  values = ensembl_ids,
  mart = ensembl
)

# Check the first few rows of the gene mapping
head(gene_mapping)

# Ensure your Seurat object and the gene mapping have the same row names
ensembl_ids <- rownames(tabula_senis)
gene_mapping <- gene_mapping %>% filter(ensembl_gene_id %in% ensembl_ids)

# Create a named vector for mapping
mapping_vector <- setNames(gene_mapping$mgi_symbol, gene_mapping$ensembl_gene_id)

# Update the row names of the Seurat object
new_rownames <- mapping_vector[ensembl_ids]

# Replace NA values with original ENSEMBL IDs (optional)
new_rownames[is.na(new_rownames)] <- ensembl_ids[is.na(new_rownames)]

# Set the new row names to the Seurat object
rownames(tabula_senis) <- new_rownames

# Check the updated row names
head(rownames(your_seurat_object))


DimPlot(tabula_senis, reduction = "umap")


Idents(tabula_senis) <- tabula_senis@meta.data$cell_type

muscle_stem_cells <- subset(tabula_senis, idents = "skeletal muscle satellite cell")

Idents(muscle_stem_cells) <- muscle_stem_cells@meta.data$age

satellite_cells <- subset(muscle_stem_cells, idents = c("3m", "24m"))


VlnPlot(satellite_cells, features = "ENSMUSG00000032231") #Anxa2 Ensembl gene ID

VlnPlot(satellite_cells, features = "ENSMUSG00000017009") #Sdc4 Ensembl gene ID

Differential_Markers <- FindMarkers(satellite_cells,
                                        ident.1 = "24m",
                                        ident.2 = "3m")

Differential_Markers["ENSMUSG00000032231",] #Anxa2 stats for 24 month vs 3 month
#                         p_val avg_log2FC pct.1 pct.2 p_val_adj
# ENSMUSG00000032231 0.07647971 -0.1135917 0.673 0.816         1

Differential_Markers["ENSMUSG00000017009",] #Sdc4 stats for 24 month vs 3 month
