


#' This function get OncoTarget genes
#' @export 
GetOncoTarget<-function(){
  
require(drugbank)
require(n1database)
  
  
data(drugbase)

pos <- sapply(names(dbase),function(x) dbase[[x]][['oncology']] & !dbase[[x]][['topical']])

dbase2 <- dbase[pos]

targets <- unique(unlist(sapply(names(dbase2),function(x) names(dbase2[[x]][['targets']]),USE.NAMES = FALSE)))

names<-NameConvertor()

tmp<-match(targets, names$entrez)

OncoTarget<-names$hgnc_symbol[tmp]

OncoTarget<-OncoTarget[is.na(OncoTarget)==FALSE]

targets<-OncoTarget

return(targets)

}



