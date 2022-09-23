


#' This function generate the protein reliability score
#' @param  vpmat the vpmat function
#' @export
proteinReliabilityScore <- function(vpmat) {
  rsc <- apply(vpmat, 1, median, na.rm=TRUE)
  rsc <- sort(rsc[is.finite(rsc)])
  rscfit <- lm(y~x, data=list(x=1:length(rsc), y=rsc), weights=1-seq(-1, 1, length=length(rsc))^2)
  rsc <- residuals(rscfit)
  rsc[1:round(length(rsc)/2)] <- -rsc[1:round(length(rsc)/2)]
  rsc[rsc<0] <- 0
  protweight <- 1-sigT(rsc, 20, .2)
  protweight <- protweight[match(rownames(vpmat), names(protweight))]
  names(protweight) <- rownames(vpmat)
  protweight[is.na(protweight)] <- 1
  return(protweight)
}


#' This function convet the value using the sigmoid function
#' @param  vpmat the vpmat function
#' @export
sigT <- function (x, slope = 20, inflection = 0.5){
  1 - 1/(1 + exp(slope * (x - inflection)))
}



#
#' DarwinOncoMatch
#'
#' This function performs oncoMatch analysis
#'
#' @param vpmat Viper matrix for the phenotypes
#' @param vpsig Viper matrix for the signatures or perturbation object: list of 2 components vpmat (viper matrix) and protweight (named vector of weights for the rows)
#' @param nn Number indicating the number of MR to consider for the enrichment analysis
#' @param modsig Optional protein activity signature for the model. If provided, only MRs conserved in the model are used for the analysis
#' @return Matrix of Normalized Enrichment Scores
#' @aliases oncoTreat
#' @export
OncoMatch <- function(vpmat, vpsig, nn=50, adaptive=FALSE, modsig=NULL, filter=TRUE) {
  if (length(vpsig)<4) {
    prweight <- vpsig[["protweight"]]
    vpsig <- vpsig[["vpmat"]]
  }
  else {
    prweight <- rep(1, nrow(vpsig))
    names(prweight) <- rownames(vpsig)
  }

 if (filter==TRUE){
  vpmat <- filterRowMatrix(vpmat, rownames(vpmat) %in% rownames(vpsig))
 }

  if (is.null(modsig)) {
    reg <- apply(vpmat, 2, function(x, prweight, nn) {
      x <- x[is.finite(x)]
      pos <- as.vector(filterColMatrix(rbind(order(x, decreasing=TRUE), order(x)), 1:floor(length(x)/2)))
      tfmode <- sign(x)[pos]
      lklh <- prweight[match(names(tfmode), names(prweight))]
      tfmode <- tfmode[!is.na(lklh)]
      lklh <- lklh[!is.na(lklh)]

      if(adaptive==FALSE){
        pos<-nn
      }else{
      pos <- which(cumsum(lklh^2)>=nn)[1]
      }

      list(tfmode=tfmode[1:pos], likelihood=lklh[1:pos])
    }, prweight=prweight, nn=nn)
  }
  else {
    reg <- apply(vpmat, 2, function(x, nn, prweight, modsig) {
      x <- x[is.finite(x)]
      f <- 2
      set <- NULL
      lh <- 0
      while(sum(lh^2)<nn*2 & f < 5) {
        tfmode <- x[order(x)[c(1:nn*f, (length(x)-nn*f+1):length(x))]]
        le <- ledge(modsig, list(ledge=list(tfmode=sign(tfmode))))[[1]]
        set1 <- tfmode[tfmode<0 & names(tfmode) %in% le]
        set2 <- tfmode[tfmode>0 & names(tfmode) %in% le]
        set <- c(set1[order(set1)[1:min(nn, length(set1))]], set2[order(set2, decreasing=TRUE)[1:min(nn, length(set2))]])
        lh <- prweight[match(names(set), names(prweight))]
        lh[is.na(lh)] <- 0
        f <- f+.1
      }
      list(tfmode=sign(set), likelihood=lh)
    }, nn=round(nn/2), prweight=prweight, modsig=modsig)
  }
  t(aREA(vpsig, reg, minsize=2)$nes)
}







#
#' DarwinOncoMatch
#'
#' This function performs oncoMatch analysis
#' The changes. The choices can be done
#'
#'
#' @param  top.mrs logic vector, whether only the top MRs are used.
#' @param vpmat Viper matrix for the phenotypes
#' @param vpsig Viper matrix for the signatures or perturbation object: list of 2 components vpmat (viper matrix) and protweight (named vector of weights for the rows)
#' @param nn Number indicating the number of MR to consider for the enrichment analysis
#' @param modsig Optional protein activity signature for the model. If provided, only MRs conserved in the model are used for the analysis
#' @return Matrix of Normalized Enrichment Scores
#' @aliases oncoTreat
#' @export
OncoMatchExtended <- function(vpmat, vpsig, nn=50, adaptive=FALSE, modsig=NULL, top.mrs=FALSE) {
  if (length(vpsig)<4) {
    prweight <- vpsig[["protweight"]]
    vpsig <- vpsig[["vpmat"]]
  }
  else {
    prweight <- rep(1, nrow(vpsig))
    names(prweight) <- rownames(vpsig)
  }

 # vpmat <- filterRowMatrix(vpmat, rownames(vpmat) %in% rownames(vpsig))

  if (is.null(modsig)) {
    reg <- apply(vpmat, 2, function(x, prweight, nn) {
      x <- x[is.finite(x)]

if(top.mrs==FALSE){

      pos <- as.vector(filterColMatrix(rbind(order(x, decreasing=TRUE), order(x)), 1:floor(length(x)/2)))

}else if (top.mrs==TRUE){

  pos <- order(abs(x), decreasing=TRUE)

}




      tfmode <- sign(x)[pos]
      lklh <- prweight[match(names(tfmode), names(prweight))]
      tfmode <- tfmode[!is.na(lklh)]
      lklh <- lklh[!is.na(lklh)]

      if(adaptive==FALSE){
        pos<-nn
      }else{
        pos <- which(cumsum(lklh^2)>=nn)[1]
      }

      list(tfmode=tfmode[1:pos], likelihood=lklh[1:pos])
    }, prweight=prweight, nn=nn)
  }
  else {
    reg <- apply(vpmat, 2, function(x, nn, prweight, modsig) {
      x <- x[is.finite(x)]
      f <- 2
      set <- NULL
      lh <- 0
      while(sum(lh^2)<nn*2 & f < 5) {
        tfmode <- x[order(x)[c(1:nn*f, (length(x)-nn*f+1):length(x))]]
        le <- ledge(modsig, list(ledge=list(tfmode=sign(tfmode))))[[1]]
        set1 <- tfmode[tfmode<0 & names(tfmode) %in% le]
        set2 <- tfmode[tfmode>0 & names(tfmode) %in% le]
        set <- c(set1[order(set1)[1:min(nn, length(set1))]], set2[order(set2, decreasing=TRUE)[1:min(nn, length(set2))]])
        lh <- prweight[match(names(set), names(prweight))]
        lh[is.na(lh)] <- 0
        f <- f+.1
      }
      list(tfmode=sign(set), likelihood=lh)
    }, nn=round(nn/2), prweight=prweight, modsig=modsig)
  }
  t(aREA(vpsig, reg, minsize=2)$nes)
}




#' DarwinOncoMatch
#'
#' This function performs oncoMatch analysis
#'
#' @param vpmat Viper matrix for the phenotypes
#' @param vpsig Viper matrix for the signatures or perturbation object: list of 2 components vpmat (viper matrix) and protweight (named vector of weights for the rows)
#' @param nn Number indicating the number of MR to consider for the enrichment analysis
#' @param balanced Logical, whether to use equal number of activated and inactivated MRs
#' @param modsig Optional protein activity signature for the model. If provided, only MRs conserved in the model are used for the analysis
#' @param reverse Logical, whether the enrichment of the drug MRs should e computed on the phenotype signature instead of the phenotype MRs on the drug signature
#' @return Matrix of Normalized Enrichment Scores
#' @export
oncoMatch_Darwin <- function(vpmat, vpsig, nn=50, balanced=TRUE, modsig=NULL, reverse=FALSE) {
  if (reverse) {
    if (length(vpsig)<4) {
      tmp <- vpmat
      vpmat <- vpsig[["vpmat"]]
      vpsig[["vpmat"]] <- tmp
    }
    else {
      tmp <- vpmat
      vpmat <- vpsig
      vpsig <- tmp
    }
    rm(tmp)
  }
  if (length(vpsig)<4) {
    prweight <- vpsig[["protweight"]]
    vpsig <- vpsig[["vpmat"]]
  }
  else {
    prweight <- rep(1, nrow(vpsig))
    names(prweight) <- rownames(vpsig)
  }
  vpmat <- filterRowMatrix(vpmat, rownames(vpmat) %in% rownames(vpsig))
  if (is.null(modsig)) {
    reg <- apply(vpmat, 2, function(x, prweight, nn, balanced) {
      x <- x[is.finite(x)]
      if (balanced) {
        pos <- as.vector(filterColMatrix(rbind(order(x, decreasing=TRUE), order(x)), 1:floor(length(x)/2)))
      }
      else {
        pos <- order(abs(x), decreasing=TRUE)
      }
      tfmode <- sign(x)[pos]
      lklh <- prweight[match(names(tfmode), names(prweight))]
      tfmode <- tfmode[!is.na(lklh)]
      lklh <- lklh[!is.na(lklh)]
      pos <- which(cumsum(lklh^2)>=nn)[1]
      list(tfmode=tfmode[1:pos], likelihood=lklh[1:pos])
    }, prweight=prweight, nn=nn, balanced=balanced)
  }
  else {
    reg <- apply(vpmat, 2, function(x, nn, prweight, modsig, balanced) {
      x <- x[is.finite(x)]
      f <- 2
      set <- NULL
      lh <- 0
      while(sum(lh^2)<nn*2 & f < 5 & 2*nn*f < length(x)) {
        if (balanced) {
          tfmode <- x[order(x)[c(1:nn*f, (length(x)-nn*f+1):length(x))]]
        }
        else {
          tfmode <- x[order(abs(x), decreasing=TRUE)[1:nn*2*f]]
        }
        le <- ledge(modsig, list(ledge=list(tfmode=sign(tfmode))))[[1]]
        set1 <- tfmode[tfmode<0 & names(tfmode) %in% le]
        set2 <- tfmode[tfmode>0 & names(tfmode) %in% le]
        set <- c(set1[order(set1)[1:min(nn, length(set1))]], set2[order(set2, decreasing=TRUE)[1:min(nn, length(set2))]])
        lh <- prweight[match(names(set), names(prweight))]
        lh[is.na(lh)] <- 0
        f <- f+.1
      }
      list(tfmode=sign(set), likelihood=lh)
    }, nn=round(nn/2), prweight=prweight, modsig=modsig, balanced=balanced)
  }
  tmp <- aREA(vpsig, reg, minsize=2)$nes
  if (reverse) return(tmp)
  t(tmp)
}







#' This function integrate NES matrices
#'
#' @param nes.list a list of nes matrices
#' @param nes.weight a weight vector for the matrices
#' @return integrated NES
#'@author Junqiang Wang
#'@export
IntegrateOncoTreatNESs<-function(nes.list, nes.weight){

nes.list<-lapply(nes.list, FUN=function(x){x<-x[order(rownames(x)), ]})

# filter rows

nes.stouffer<-nes.list[[1]]

drugs.list<-lapply(nes.list, rownames)

nes.combined <- do.call(cbind,nes.list)

for (i in 1:nrow(nes.list[[1]])){
  for (j in 1:ncol(nes.list[[1]])){

   tmp<- which(colnames(nes.combined)  %in% colnames(nes.list[[1]])[j]==TRUE)
   tmp<-nes.combined[i, tmp]

   tmp<-sum(tmp*nes.weight, na.rm=TRUE)/sqrt(sum((nes.weight^2)))

   nes.stouffer[i,j]<-tmp

  }

}

return(nes.stouffer)

}












#' This function plot the oncoTreat/oncoMatch result
#'@param oncoTreat the OncoTreat/OncoMatch output; a NES matrix
#'@param nn Number of drugs for each phenotype
#'@return A Heatmap with p-values
#' @export
PlotOncoTreat<-function(oncoTreat,
                        top.n=10,
                        round.n=0,
                        break.legend=10,
                        geom.text.size=4,

                        labs.x="",
                        remove.concentration.at.10=TRUE
                        ){

  require(RColorBrewer)
  require(ggplot2)
  require(reshape2)

  #
  nes<-as.matrix(oncoTreat)




  # remove concentration at

  if(remove.concentration.at.10==TRUE){

    nes<-nes[!grepl("_10_", rownames(nes)),]

  }





  pval<-pnorm(nes, lower.tail = TRUE)


  #adjust the p-value by FDR

  unlist(pval)

  pval.fdr<-matrix(as.numeric(p.adjust(unlist(pval), method="fdr")), nrow=nrow(pval))

  rownames(pval.fdr)<-rownames(pval)

  colnames(pval.fdr)<-colnames(pval)

  #

  pval.fdr<--log10(pval.fdr)


  pval.fdr<-as.data.frame(pval.fdr)

  tmp<-apply(pval.fdr,2, as.numeric)
  rownames(tmp)<-rownames(pval.fdr)

  pval.fdr<-tmp

  #pval.fdr<-pval.fdr[,-1]


  # chose the top 10 and do heatmap

  dat<-as.data.frame(pval.fdr)

  dat$drugs<-rownames(dat)

  dat<-melt(dat)


  dat %>%
    arrange(desc(value)) %>%
    arrange(variable) %>%
    group_by(variable) %>%
    dplyr::filter(value>2)%>%
    top_n(n = top.n, wt=value) %>%
     as.data.frame() %>%
    distinct(drugs, .keep_all= TRUE)  -> all.drugs


  drugs<-all.drugs$drugs

  drugs<-pval.fdr[match(drugs, rownames(pval.fdr)),]


  #--do heatmap no

  require(RColorBrewer)
  require(ggplot2)
  require(reshape2)

  # convert to long format
  df <- melt(drugs)
  summary(df)

  df$log10FDR<-df$value


  p1 <- ggplot(df, aes(x = Var2, y = Var1, fill = log10FDR)) +
    geom_raster() +
    geom_text(aes(label = round(value, 0)), size=geom.text.size)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ggtitle("OncoTreat")+
    theme(plot.title = element_text(hjust = 0.5))
  #theme_minimal(base_size = 16)


  limit <- c(0, max(df$log10FDR))

  d1<-1/max(df$log10FDR)

  d2<-2/max(df$log10FDR)
  d3<-3/max(df$log10FDR)

  d4<-4/max(df$log10FDR)
  d5<-5/max(df$log10FDR)


  b0<-seq(0, 200, by=break.legend)
  b1<-rep(floor(max(df$log10FDR)), length(b0))

  tmp<-which(b0>b1)[1]

  #tmp<-max(1, tmp-2)

  #breaks=c(b0[1:tmp], floor(max(df$log10FDR)))

  breaks=b0

  p2<-p1 + scale_fill_distiller(palette = 1, values=c(0,d1,d2,d3, d4,d5, 1), breaks=breaks, limit = limit, direction=1) +labs(fill="-log10(FDR)") + labs(y="Drugs")+labs(x=labs.x) +scale_colour_brewer(palette="BuPu", direction=-1)


   #scale_fill_distiller(palette = 1, values=c(0,d1,d2,d3, d4,d5, 1), breaks=breaks, limit = limit, direction=1) +labs(fill="-log10(FDR)") + labs(y="Drugs")+labs(x="Cellular States") +scale_colour_brewer(palette="BuPu", direction=-1)


  p2


}


#' This function plot the oncoTreat/oncoMatch result
#'@param oncoTreat the OncoTreat/OncoMatch output; a NES matrix
#'@param nn Number of drugs for each phenotype
#'@return A Heatmap with p-values
#' @export
TxtOncoTreat<-function(oncomatch.out, nn=10, round.n=0){

  require(RColorBrewer)
  require(ggplot2)
  require(reshape2)

  #
  nes<-as.matrix(oncomatch.out)

  pval<-pnorm(nes, lower.tail = TRUE)

  #adjust the p-value by FDR

  unlist(pval)

  pval.fdr<-matrix(as.numeric(p.adjust(unlist(pval), method="fdr")), nrow=nrow(pval))

  rownames(pval.fdr)<-rownames(pval)

  colnames(pval.fdr)<-colnames(pval)


  #pval.fdr<--log10(pval.fdr)

  pval.fdr<-as.data.frame(pval.fdr)

  tmp<-apply(pval.fdr,2, as.numeric)
  rownames(tmp)<-rownames(pval.fdr)

  pval.fdr<-tmp

  #pval.fdr<-pval.fdr[,-1]
  # chose the top 10 and do heatmap

  dat<-as.data.frame(pval.fdr)

  dat$drugs<-rownames(dat)

  dat<-melt(dat)


  dat %>%
    arrange(desc(-value)) %>%
    arrange(variable) %>%
    group_by(variable) %>%
   # dplyr::filter(value>1)%>%
    top_n(n = nn, wt= -value) -> all.drugs


  drugs<-all.drugs$drugs

  drugs<-pval.fdr[match(drugs, rownames(pval.fdr)),]

return(drugs)

}









#' This function plot the oncoTreat/oncoMatch result
#'@param oncoTreat the OncoTreat/OncoMatch output; a NES matrix
#'@param nn Number of drugs for each phenotype
#'@return A Heatmap with p-values
#' @export
PlotOncoMatch<-function(oncoMatch,
                        top.n=50,
                        break.legend=5,
                        lower.tail=FALSE,
                        labs.x="",
                        labs.y="",
                        geom.text.size=5,
                        plot.all=c(FALSE, TRUE),
                        return.plot=c(FALSE, TRUE),
                        return.dat=c(FALSE, TRUE)
                        ){

  require(RColorBrewer)

  require(ggplot2)

  require(reshape2)

  nes<-as.matrix(oncoMatch)

  pval<-pnorm(nes, lower.tail = lower.tail)

  #adjust the p-value by FDR

  unlist(pval)

  pval.fdr<-matrix(as.numeric(p.adjust(unlist(pval), method="fdr")), nrow=nrow(pval))

  rownames(pval.fdr)<-rownames(pval)

  colnames(pval.fdr)<-colnames(pval)

  #

  pval.fdr<--log10(pval.fdr)


  pval.fdr<-as.data.frame(pval.fdr)

  tmp<-apply(pval.fdr,2, as.numeric)
  rownames(tmp)<-rownames(pval.fdr)

  pval.fdr<-tmp


  # chose the top 10 and do heatmap

  dat<-as.data.frame(pval.fdr)

  dat$drugs<-rownames(dat)

  dat<-melt(dat)


  dat %>%
    arrange(desc(value)) %>%
    arrange(variable) %>%
    group_by(variable) %>%
   # dplyr::filter(value>1)%>%
    top_n(n = top.n, wt=value)%>%
    as.data.frame() %>%
  distinct(drugs, .keep_all= TRUE)  -> all.drugs


  if(plot.all==TRUE){

    all.drugs<-dat
  }


  drugs<-all.drugs$drugs

  drugs<-pval.fdr[match(drugs, rownames(pval.fdr)),]

  if(return.dat==TRUE) {return(drugs)}

  #--do heatmap no

  require(RColorBrewer)
  require(ggplot2)
  require(reshape2)


  # convert to long format
  df <- melt(drugs)
  summary(df)

  df$log10FDR<-df$value




  p1 <- ggplot(df, aes(x = Var2, y = Var1, fill = log10FDR)) +
    geom_raster() +
    geom_text(aes(label = round(value, 0)), size=geom.text.size)+
    theme(axis.text.x = element_text(angle = 30, hjust=1))+
    ggtitle("OncoMatch")+
    theme(plot.title = element_text(hjust = 0.5))

  #theme_minimal(base_size = 16)


  limit <- c(0, max(df$log10FDR))

  d1<-1/max(df$log10FDR)

  d2<-2/max(df$log10FDR)
  d3<-3/max(df$log10FDR)

  d4<-4/max(df$log10FDR)
  d5<-5/max(df$log10FDR)


  b0<-seq(0, 200, by=break.legend)

  b1<-rep(floor(max(df$log10FDR)), length(b0))

  tmp<-which(b0>b1)[1]

  tmp<-max(1, tmp-2)

  #breaks=c(b0[1:tmp], floor(max(df$log10FDR)))

 breaks<-b0

  p2<-p1 +

    scale_fill_distiller(palette = 1, values=c(0,d1,d2,d3, d4,d5, 1), breaks=breaks, limit = limit, direction=1) +labs(fill="-log10(FDR)") + labs(y=labs.y)+labs(x=labs.x) +scale_colour_brewer(palette="BuPu", direction=-1)


  if(return.plot==TRUE){return(p2)}

   p2


}





#' This function plot the oncoTreat/oncoMatch result
#'@param oncoTreat the OncoTreat/OncoMatch output; a NES matrix
#'@param nn Number of drugs for each phenotype
#'@return A Heatmap with p-values
#' @export
GetOncoTreat<-function(oncoTreat, nn=20){

  require(RColorBrewer)
  require(ggplot2)
  require(reshape2)
require(dplyr)

  nes<-oncoTreat

  dat<-nes

dat<-as.data.frame(dat)

dat$drugs<-rownames(dat)

  dat<-melt(dat)


  dat %>%
    arrange(desc(-value)) %>%
    arrange(variable) %>%
    group_by(variable) %>%
    top_n(n = -nn, wt=value) -> all.drugs


  drugs<-all.drugs$drugs

  drugs<-nes[match(drugs, rownames(nes)),]




}





