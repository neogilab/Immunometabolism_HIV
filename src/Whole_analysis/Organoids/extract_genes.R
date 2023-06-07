
## Custom functions used in the analysis should go into this chunk.
## They will be listed in their own section of the appendix.
#' Function to create UpSetR compatible table from taxonomy table
#'
#' @param data taxonomy data 
#' @param cname name of the column with taxon of interest
#' @param ename name of the column with entity names
#'
#' @return UpSetR compatible data.frame
prepareUpSet<-function(data,cname,ename='GeneID'){
  data<-as.data.frame(data)
  myFromList<-function (input) 
  {
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    rownames(data)<-elements
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    return(data)
  }
  if(any(is.na(data[,cname]))){
    data[is.na(data[,cname]),cname]<-'Unspecified'
  }
  sets<-unique(data[,cname])
  l<-lapply(sets, function(.x){ 
    as.character(data[data[,cname]==.x,ename])
  })
  names(l)<-sets
  pt<-myFromList(l)
}
