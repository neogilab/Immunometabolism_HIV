#!/usr/bin/env Rscript

install.packages("BiocManager")

BiocManager::install(c("SNFtool", "netDx", "ggplots2", "limma", "DESeq2","ComplexHeatmap","Seurat","Hmisc","piano","ggalluvial","caret", "randomForest", "Boruta"))




library("SNFtool")
library("netDx")
library("ggplots2")
library("limma")
library("DESeq2")
library("ComplexHeatmap")
library("Seurat")
library("Hmisc")
library("piano")
library("ggalluvial")
library("caret")
library("randomForest")
library("Boruta")
