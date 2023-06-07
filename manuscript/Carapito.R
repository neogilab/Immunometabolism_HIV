### FUNCTIONS
calBICdiff <- function(net,data) {
  arcList = bnlearn::arc.strength(net,data,criterion='bic-g')
  return (arcList)
}
#####################################################################################################
extractGraph <- function(bn.hc_dataSets=NULL, label_string="")
{
  #' Extract arc, number of descendant, and number of ancestor from final DAG
  #' @param bn.hc_dataSets 
  #' @param label_string base name 
  
  # extract number of childnodes and rank nodes 
  childlist <- lapply(names(bn.hc_dataSets$nodes), function(x) bnlearn::descendants(bn.hc_dataSets, x) )
  names(childlist) <- names(bn.hc_dataSets$nodes)
  childLs_sortIndex <- sort.int(unlist(lapply(childlist, length)), decreasing = TRUE, index.return=TRUE)
  childLs_sorted <- childlist[childLs_sortIndex$ix]
  
  # extract number of ancestors and rank nodes
  ancestorlist <- lapply(names(bn.hc_dataSets$nodes), function(x) bnlearn::ancestors(bn.hc_dataSets, x) )
  names(ancestorlist) <- names(bn.hc_dataSets$nodes)
  ancestorLs_sortIndex <- sort.int(unlist(lapply(ancestorlist, length)), decreasing = TRUE, index.return=TRUE)
  ancestorLs_sorted <- ancestorlist[ancestorLs_sortIndex$ix]
  
  # extract count of descendants and ancestors; merge two 
  descendantsCount <- sapply(childLs_sorted, function(x) length(x) )
  descendantsCountDF <-data.frame(geneSymbol=names(descendantsCount), descendantCount=descendantsCount,stringsAsFactors = F ) 
  
  ancestorsCount <- sapply(ancestorLs_sorted, function(x) length(x) )
  ancestorsCountDF <-data.frame(geneSymbol=names(ancestorsCount), ancestorCount=ancestorsCount,stringsAsFactors = F ) 
  totalCountDF <- merge(descendantsCountDF, ancestorsCountDF, by ="geneSymbol")
  
  # arc list
  arcList <- as.data.frame(bnlearn::arcs(bn.hc_dataSets), stringsAsFactors=F)
  # write.table(arcList,file=file.path(output_folder,'tables', paste0("arcs_", label_string, ".csv")), sep = ",", row.names = F)
  
  # merge count DF and arc List
  finalDF <- arcList
  finalDF$fromDescendantCount <- descendantsCountDF$descendantCount[match(finalDF$from, descendantsCountDF$geneSymbol)]
  finalDF$fromAncestorCount <- ancestorsCountDF$ancestorCount[match(finalDF$from, ancestorsCountDF$geneSymbol)]
  finalDF$toDescendantCount <- descendantsCountDF$descendantCount[match(finalDF$to, descendantsCountDF$geneSymbol)]
  finalDF$toAncestorCount <- ancestorsCountDF$ancestorCount[match(finalDF$to, ancestorsCountDF$geneSymbol)]
  
  write.table(finalDF,file=file.path(output_folder,'tables', paste0("arc2count_", label_string, ".txt")), quote=F, sep = "\t", row.names = F)
  
  return (finalDF)
}