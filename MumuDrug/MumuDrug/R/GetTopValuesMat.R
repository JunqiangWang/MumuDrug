



#' This function get top master regulators 
#' @param  top_n the top n values needed
#' @param  order.cluster the clusters ordered
#' @export
GetTopValuesMat<-function(mat,
                          top_n=10,
                          remove.duplicates=TRUE,
                          order.cluster=NULL
                          ){


require(dplyr)
  
mat<-mat %>% as.matrix()  
  
  
  
top.n.markers<-vector()

cluster.names<-colnames(mat)

for (i in cluster.names){
  
  print(i)
  
  tmp<-mat[,i]
  
  tmp<-sort(tmp, decreasing=TRUE)
  
  tmp<-sort(tmp, decreasing=TRUE)[1:top_n]  %>% as.data.frame()
  
  tmp<-cbind(tmp, i)
  
  colnames(tmp)<-c("Values", "clusters")
  
  tmp$label<-rownames(tmp)
  
  top.n.markers<-rbind(top.n.markers, tmp)
  
}



# The duplicate will be removed

if(remove.duplicates==TRUE){

top.n.markers %>%
  arrange(-Values) %>%
  distinct(label, .keep_all = TRUE) %>%
  arrange(clusters, -Values)->top.n.markers
  
}




tmp.markers<-vector()

for (i in cluster.names){
  
  tmp<-top.n.markers %>%
    dplyr::filter(clusters==i)
  
  tmp.markers<-rbind(tmp.markers, tmp)
  
}

top.n.markers<-tmp.markers

tmp<-match(top.n.markers$label, rownames(mat))

mat<-mat[tmp, ]

return(mat)

}












