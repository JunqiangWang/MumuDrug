

#' This function generate the pseudo signature from the drug perturbed signature
#' @param  vpmat The protein activity signature of the perturbed drugs
#' @return  The pseudo signatures of two combined drugs
#' @author Junqiang Wang
#' @export


PseudoSignatureCombinedDrug<-function(vpmat, method="maximalEffects"){


  require(stringr)


  vpmat<-vpmat

drugs<-colnames(vpmat)

drugs<-as.vector(str_split_fixed(drugs, pattern="_", 2)[,1])



pseudo.signature<-vector()

pseudo.signature.names<-vector()

message("process the ith column:")

for (i in 1:ncol(vpmat)){

  print(i)

for (j in 1:ncol(vpmat)){

tmp<-vpmat[, c(i,j)]

tmp<-apply(tmp, 1, FUN = function(x){

if (abs(x[1])>abs(x[2])){x<-x[1]}
else{
x<-x[2]}

})

pseudo.signature<-cbind(pseudo.signature, tmp)

tmp.name<-paste0(drugs[i], ":",  drugs[j])

pseudo.signature.names<-c(pseudo.signature.names, tmp.name)

}

}

colnames(pseudo.signature)<-pseudo.signature.names

return(pseudo.signature)

}

















