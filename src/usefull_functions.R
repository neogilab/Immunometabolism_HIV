#' Run limma and extract table with logFC and pvalue for all proteins
#'
#' @param frame dataframe containing proteim abundances 
#' @param group vector containing in order conditions of observations
#' @param name_ana name of analysis
#' @param comp comparison done
#' @param c1 condition control
#' @param c2 condition case
#' @return dataframe containing uniprot ID LogFC and statistics
extract_limma <- function(fit, design, name_ana, c1, c2, viral_proteins){
  print(names(frame))
  print(c(c1,c2))
  ## limma
  ## pvalue = 0.05
  ## adjustment Benjamini-Hotschberg
  # construct of groups for model design
  # Fit linear model : estimate fold change and standard errors
  cont.matrix <- makeContrasts(paste0(c2, "-", c1), levels=design)
  print(cont.matrix)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  # Execute eBayes to obtain attenuated t statistics(smoothing standard errors)
  fit.cont <- eBayes(fit.cont)
  top_table <- topTable(fit.cont, adjust="fdr", n = Inf, sort.by = "P")
  top_table$Accession <- rownames(top_table)
  top_table_virus <- top_table[top_table$Accession %in% viral_proteins, ]
  name <- paste0(c1, "_vs_", c2, "_DGE_results_virus")
  save_file_excel(top_table_virus, name, paste0("results/virus_levels/"))
  top_table_2 <- merge(annotationInfo, top_table, by = "Accession")
  top_table_2 <- top_table_2[!duplicated(top_table_2$Gene_name), ]
  print(top_table_2)
  name <- paste0(c1, "_vs_", c2, "_DGE_results")
  print(name)
  save_file_excel(top_table_2, name, paste0("results/", name_ana, "/DGE/"))
  top_table_3 <- top_table_2
  name <- paste0(c1, "_vs_", c2, "_DGE_results_filt")
  top_table_3 <-top_table_3[top_table_3$adj.P.Val < 0.05, ]
  print(c2)
  if(nrow(top_table_3) != 0){
    top_table_3$Significance <- ifelse(top_table_3$logFC < 0, paste0("downregulated in ", c2), paste0("upregulated in ", c2))
    top_table_3$Significance_2 <- ifelse(top_table_3$logFC < 0, paste0("downregulated"), paste0("upregulated"))
    save_file_excel(top_table_3, name, paste0("results/", name_ana, "/DGE/"))
  }
  return(top_table_3)
}

#' Save figure as PDF in results
#'
#' @param figure variable containing figure
#' @param name name of the figure
#' @param height height of the figure
#' @param width width of the figure
save_figure <- function(figure, name, height, width){
  path_fig <- paste0("results/", name_ana, "/figures/", name_ana, "_", name,".pdf")
  dev.copy(pdf, path_fig, width = width, height = width)
}

#' Save figure as PDF in results
#'
#' @param ggplot variable containing figure
#' @param name name of the figure
#' @param height height of the figure
#' @param width width of the figure
save_ggplot <- function(name, height, width){
  path_fig <- paste0("results/figures/", name_ana, "_", name,".pdf")
  ggsave(path_fig, width = width, height = width)
  dev.off()
}

#' Save file as csv
#'
#' @param file variable containing file
#' @param name name of the file
#' @param destination destination folder for file
save_file_csv <- function(file, name, destination){
  path_file <- paste0(destination, name_ana, "_", name,".csv")
  print(path_file)
  write.csv(file, file = path_file)
}


#' Save file as txt
#'
#' @param file variable containing file
#' @param name name of the file
#' @param destination destination folder for file
save_file <- function(file, name, destination){
  path_file <- paste0(destination, name_ana, "_", name,".txt")
  print(path_file)
  write.table(file, file =path_file, sep = "\t", row.names = TRUE, 
              col.names = NA)
}

#' Save file as txt
#'
#' @param file variable containing file
#' @param name name of the file
#' @param destination destination folder for file
save_file_excel <- function(file, name, destination){
  path_file <- paste0(destination, name_ana, "_", name,".xlsx")
  print(path_file)
  write.xlsx(file, file =path_file)
}

#' Set up data_matrix for NormalizerDE
#'
#' @param filtered_file file after filtration ready for normalization
#' @param name_ana name of analysis
make_data_matrix <- function(filtered_file,  name_ana){
  rownames(filtered_file) <- c()
  path_data <- paste0("processing/normalization/", name_ana, "_data_matrix.tsv")
  print(path_data)
  write.table(filtered_file, file = path_data, sep = "\t",quote = FALSE, row.names=FALSE)
  return(path_data)
}

#' Set up design_matrix for NormalizerDE
#'
#' @param filtered_file file after filtration ready for normalization
#' @param name_analysis name of analysis
make_design_matrix <- function(filt_file, name_analysis, group){
  filt_file$Accession <- NULL
  design_table <-data.frame(sample = colnames(filt_file), group = group)
  path_design <- paste0("processing/normalization/", name_analysis, "_design_matrix.tsv")
  write.table(design_table, file = path_design,
              sep = "\t",quote = FALSE, row.names=FALSE)
  return(path_design)
}

#' Run limma and extract table with logFC and pvalue for all proteins
#'
#' @param frame dataframe containing proteim abundances 
#' @param group vector containing in order conditions of observations
#' @param name_ana name of analysis
#' @param comp comparison done
#' @param c1 condition control
#' @param c2 condition case
#' @return dataframe containing uniprot ID LogFC and statistics
run_limma <- function(frame, group, name_ana, comp, c1, c2){
  print(names(frame))
  print(c(c1,c2))
  ## limma
  ## pvalue = 0.05
  ## adjustment Benjamini-Hotschberg
  # construct of groups for model design
  groups <- as.factor(group)
  design <- model.matrix( ~ 0 + groups)
  colnames(design) <- levels(groups) # change colnames to KO and WT
  
  # Fit linear model : estimate fold change and standard errors
  fit <- lmFit(frame, design)
  # Fit linear model : estimate fold change and standard errors
  cont.matrix <- makeContrasts(paste0(c2, "-", c1),levels=design)
  print(cont.matrix)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  # Execute eBayes to obtain attenuated t statistics(smoothing standard errors)
  fit.cont <- eBayes(fit.cont)
  top_table <- topTable(fit.cont, adjust="fdr", n = Inf, sort.by = "P")
  #top_table <-top_table[top_table$adj.P.Val < 0.05, ]
  #top_table$Accession <- rownames(top_table)
  #data <-merge(top_table, annotationInfo, by = "Accession")
  #data <- select(data, Gene_name, logFC, adj.P.Val)
  return(top_table)
}

#' Run limma and extract table with logFC and pvalue for all proteins
#'
#' @param frame dataframe containing proteim abundances 
#' @param group vector containing in order conditions of observations
#' @param name_ana name of analysis
#' @param comp comparison done
#' @param c1 condition control
#' @param c2 condition case
#' @return dataframe containing uniprot ID LogFC and statistics
run_limma_2 <- function(frame, group, name_ana, comp, ca, cb){
  print(names(frame))
  print(c(ca,cb))
  ca = ca
  cb = cb
  ## limma
  ## pvalue = 0.05
  ## adjustment Benjamini-Hotschberg
  # construct of groups for model design
  groups <- as.factor(group)
  design <- model.matrix( ~ 0 + groups)
  colnames(design) <- levels(groups) # change colnames to KO and WT
  
  # Fit linear model : estimate fold change and standard errors
  fit <- lmFit(frame, design)
  # Fit linear model : estimate fold change and standard errors
  cont.matrix <- makeContrasts(paste0(cb, "-", ca),levels=design)
  print(cont.matrix)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  # Execute eBayes to obtain attenuated t statistics(smoothing standard errors)
  fit.cont <- eBayes(fit.cont)
  top_table <- topTable(fit.cont, adjust="fdr", n = Inf, sort.by = "P")
  top_table <-top_table[top_table$adj.P.Val < 0.1, ]
  top_table$Accession <- rownames(top_table)
  data <-merge(top_table, annotationInfo, by = "Accession")
  data <- select(data, Gene_name, logFC, adj.P.Val)
  return(data)
}

#' Print heatmap of genes differentially expressed
#'
#' @param r_norm dataframe containing proteim abundances of conditions to compare
#' @param DF_top_table_2 Filt DGE results table containing uniprot ID
#' @param name_ana name of analysis
#' @param comp comparison done
#' @param annotationInfo table containing uniprot accession and gene name
#' @param label_cutoff name FDR cutoff value
clustering_DGE <- function(r_norm, DF_top_table_2, name_ana, comp, 
                           annotationInfo){
  # get subset of genes differentially expressed 
  subset <- r_norm[rownames(r_norm) %in% DF_top_table_2$Accession,]
  ## save proteins abundance of significant proteins
  subset$Accession <- rownames(subset)
  subset <-merge(subset, annotationInfo, by = "Accession")
  subset <- subset[!duplicated(subset$Gene_name), ]
  rownames(subset) <- subset$Gene_name
  subset$Gene_name <- NULL
  subset$Accession <- NULL
  hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
  path_fig <- paste0("results/", name_ana, "/figures/",  name_ana, "_", comp, 
                     "_clustering.pdf")
  pdf(path_fig)
  heatmap.2(as.matrix(subset), scale="row", Rowv = TRUE, Colv = TRUE, trace="none",
            col = hmcol, keysize = 0.9, dendrogram = 'both')
  dev.off()
  
}

