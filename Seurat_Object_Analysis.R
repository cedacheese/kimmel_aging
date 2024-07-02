library(Seurat)
library(SeuratData)
library(Matrix)
library(patchwork)
library(dplyr)

setwd("/Users/ceda_/Documents/PhD/Coding/R/kimmel_aging")

mtx_file <- "/Users/ceda_/Documents/PhD/Coding/R/kimmel_aging/GSE143476_matrix.mtx"
barcodes_file <- "/Users/ceda_/Documents/PhD/Coding/R/kimmel_aging/GSE143476_barcodes.tsv"
genes_file <- "/Users/ceda_/Documents/PhD/Coding/R/kimmel_aging/GSE143476_genes.tsv"
cell_metadata_file <- "/Users/ceda_/Documents/PhD/Coding/R/kimmel_aging/GSE143476_cell_metadata.csv"

#Read the expression matrix file
expression_matrix <- readMM(mtx_file)

# Check the dimensions of the expression matrix
dim(expression_matrix)
# 27998 rows (genes) and 22699 columns (cells)

#Read the barcodes and gene names
barcodes <- readLines(barcodes_file)

# Creates a table of the 27998 unique genes in this dataset, split between their ENSEMBL ID and Gene ID
genes <- read.table(genes_file, header = FALSE, sep = "\t")

# Returns 22699 length (cells), which matches with the expression matrix dimensions
length(barcodes)

# Returns 27998 length (genes), which matches with the expression matrix dimensions
length(genes$V2)

# Read the cell metadata
cell_metadata <- read.csv(cell_metadata_file)

# Returns (among other data) the class (data.frame) and 21555 observations across 18 variables. Each observation corresponds to a cell barcode, which suggests that this file contains the metadata for 21,555 cells from the dataset.
# However, the original expression matrix has 22,699 cells, so now we have to filter out the barcodes which are not present in the metadata file.
str(cell_metadata)

# Filter barcodes
common_barcodes <- barcodes[barcodes %in% cell_metadata$X]

# Filter the expression matrix to include only these barcodes
barcode_indices <- which(barcodes %in% common_barcodes)
filtered_expression_matrix <- expression_matrix[, barcode_indices]

# Create new object of filtered barcodes for downstream
filtered_barcodes <- barcodes[barcode_indices]

# All non-duplicated genes from the 
unique_genes <- !duplicated(genes$V2)

# Modifies the object, changing dimensions from 27,998 genes to 27,933 genes, removing 65 duplicated genes
genes <- genes[unique_genes, ]

# Modifies the filtered expression matrix to remove duplicated genes
filtered_expression_matrix <- filtered_expression_matrix[unique_genes, ]

# Create the Seurat object
kimmel_object <- CreateSeuratObject(counts = filtered_expression_matrix)

#Add gene and cell names
rownames(kimmel_object) <- genes$V2
colnames(kimmel_object) <- filtered_barcodes

# Add cell metadata to the Seurat object
kimmel_object <- AddMetaData(kimmel_object, metadata = cell_metadata)

# Perform standard Seurat preprocessing and analysis steps
kimmel_object <- NormalizeData(kimmel_object)
kimmel_object <- FindVariableFeatures(kimmel_object)
kimmel_object <- ScaleData(kimmel_object)
kimmel_object <- RunPCA(kimmel_object)
kimmel_object <- FindNeighbors(kimmel_object)
kimmel_object <- FindClusters(kimmel_object)
kimmel_object <- RunUMAP(kimmel_object, dims = 1:10)

# Visualize the results
DimPlot(kimmel_object, reduction = "umap")

# Set identity class of the kimmel object to activation
Idents(kimmel_object) <- kimmel_object$Activation

# Subset the quiescent MuSCs from kimmel object
quiescent <- subset(kimmel_object, idents = "Q")

# Set identity class of the kimmel object to age
Idents(quiescent) <- quiescent$Age

# DEGs analysis of the quiescent subset between Young and Aged (Aged is the referential point, e.g., negative fold change means transcription of the gene goes down in aged)
differential_quiescent <- FindMarkers(quiescent,
                                        ident.1 = "Aged",
                                        ident.2 = "Young")


# Creates character vectors of the cell barcodes for all aged/young cells
aged_cells <- WhichCells(quiescent, ident = "Aged")
young_cells <- WhichCells(quiescent, ident = "Young")

# Extract expression data for Sdc4 from the quiescent object
sdc4_expression <- FetchData(quiescent, vars = "Sdc4")

# Get the expression levels for aged and young cells using the character vector of cell barcodes
sdc4_aged <- sdc4_expression[aged_cells, "Sdc4"]
sdc4_young <- sdc4_expression[young_cells, "Sdc4"]

# Perform Wilcoxon rank-sum test
sdc4_wilcox_test <- wilcox.test(sdc4_aged, sdc4_young)

sdc4_wilcox_test$p.value

# View results
wilcox_test

library(ggplot2)

anxa2_expression <- FetchData(quiescent, vars = "Anxa2")
anxa2_aged <- anxa2_expression[aged_cells, "Anxa2"]
anxa2_young <- anxa2_expression[young_cells, "Anxa2"]

# Create a data frame for plotting Anxa2
anxa2_data <- data.frame(
    Expression = c(anxa2_young, anxa2_aged),
    Group = factor(c(rep("Young", length(anxa2_young)), rep("Aged", length(anxa2_aged))),
                        levels = c("Young","Aged"))
)

anxa2_data$Group <- as.factor(anxa2_data$Group)
anxa2_data <- factor(x = anxa2_data$Group, levels = c("Young","Aged"))

anxa2_wilcox_test <- wilcox.test(anxa2_aged, anxa2_young)

anxa2_wilcox_test

Quiescent_Anxa2 <- ggplot(anxa2_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0,4)) +
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "red") +
    labs(title = "Expression of Anxa2: Young vs. Aged",
        x = "Group",
        y = "Normalized Expression") +
    theme(panel.background = element_rect(fill = "white", color = NA),
        panel.ontop = FALSE,
        panel.grid.major = element_blank(),
        axis.line.y = element_line(color = "black",
                                    linewidth = 0.5),
        axis.line.x = element_line(color = "black",
                                    linewidth = 0.5),
        plot.title = element_text(hjust = 0.5,
                                    size = 25,
                                    face = "bold"),
        axis.title.y = element_text(size = 25,
                                    face = "italic"),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20)) +
        NoLegend()

pdf(file = "Anxa2_Young_Aged.pdf",
    width = 10,
    height = 30)
plot(Quiescent_Anxa2)
dev.off()

# Create a data frame for plotting
sdc4_data <- data.frame(
  Expression = c(sdc4_young, sdc4_aged),
  Group = factor(c(rep("Young", length(sdc4_young)), rep("Aged", length(sdc4_aged))), 
                        levels = c("Young","Aged"))
)

sdc4_data$Group <- as.factor(sdc4_data$Group)
sdc4_data <- factor(x = sdc4_data$Group, levels = c("Young","Aged"))

# Plot the data
Quiescent_Sdc4 <- ggplot(sdc4_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "red") +
  labs(title = "Expression of Sdc4: Young vs. Aged", 
       x = "Condition", 
       y = "Normalized Expression") +
 theme(panel.background = element_rect(fill = "white", color = NA),
        panel.ontop = FALSE,
        panel.grid.major = element_blank(),
        axis.line.y = element_line(color = "black",
                                    linewidth = 0.5),
        axis.line.x = element_line(color = "black",
                                    linewidth = 0.5),
        plot.title = element_text(hjust = 0.5,
                                    size = 25,
                                    face = "bold"),
        axis.title.y = element_text(size = 25,
                                    face = "italic"),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20)) +
        NoLegend()

# Save the plot as pdf
pdf(file = "Sdc4_Young_Aged.pdf",
    width = 10,
    height = 30)
plot(Quiescent_Sdc4)
dev.off()

# Creates data table of differential markers for the MuSCs in the dataset, with the referential comparison being the aged samples
differential_markers <- FindMarkers(kimmel_object,
                                        ident.1 = "Aged",
                                        ident.2 = "Young")

# Returns the statistical values of gene transcribing Syndecan-4 from the differential analysis between Aged and Young
sdc4_stats <- differential_markers["Sdc4",]
#             p_val avg_log2FC pct.1 pct.2   p_val_adj
# Sdc4 3.045018e-47 -0.3253447 0.874 0.918 8.50565e-43

# Drop in percentage of Sdc4+ MuSCs in Age
(sdc4_stats$pct.2 - sdc4_stats$pct.1)*100
# 4.4

FeaturePlot(kimmel_object, features = "Anxa2")

# Returns the statistical values of gene transcribing Annexin-A2 from the differential analysis between Aged and Young
anxa2_stats <- differential_markers["Anxa2",]
#              p_val avg_log2FC pct.1 pct.2    p_val_adj
# Anxa2 2.255446e-22 -0.3242899 0.433 0.486 6.300137e-18

(anxa2_stats$pct.2 - anxa2_stats$pct.1)*100
# 5.3

# Setting age column in meta.data as factor to manipulate order of variables
kimmel_object@meta.data$Age <- as.factor(kimmel_object@meta.data$Age)

# 
kimmel_object@meta.data$Age <- factor(x = kimmel_object@meta.data$Age, levels = c("Young", "Aged"))

Idents(kimmel_object) <- kimmel_object@meta.data$Age

VlnPlot(kimmel_object, features = c("Anxa2","Sdc4"), pt.size = 0, y.max = 6)

plot
plot + theme
seurat_cell_names <- colnames(kimmel_object)
 
metadata_cell_names <- rownames(cell_metadata)

all(seurat_cell_names %in% metadata_cell_names)



duplicate_genes <- genes$V2[duplicated(genes$V2)]
print(duplicate_genes)
