library(Seurat)
library(SeuratData)
library(Matrix)
library(patchwork)
library(dplyr)
library(readr)

setwd("/Users/ceda_/Documents/PhD/Coding/R/shcherbina_aging")

df <- readr::read_tsv("GSE121589_bulkRNAseq.TPM.tsv")

colnames(df)
rownames(df)
head(df, 10:10)
df[,GeneName == "Anxa2"]

rownames(df) <- df$GeneName

nrow(df)

smallestGroupSize <- 4
keep <- rowSums(counts(df) >= 10) >= smallestGroupSize
df <- df[keep,]
nrow(df)

??counts

colData(df)
