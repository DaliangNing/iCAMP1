quali.perf<-function(exp.v=NULL,PNmexp=NULL,obs.v=NULL,PNm=NULL,ab.v,pron=c("HoS","HeS","HD","DL","DR"))
{
  if(length(unique(c(length(exp.v),nrow(PNmexp),length(obs.v),nrow(PNm),length(ab.v),0)))!=2){stop("Different length!")}
  TFm<-matrix(NA,nrow=length(ab.v),ncol=length(pron))
  if(is.null(PNmexp))
  {
    PNmexp<-TFm
    PNmexp[]<-"N"
    PNmexp[cbind(1:nrow(PNmexp),match(exp.v,pron))]="P"
  }else{PNmexp=as.matrix(PNmexp)}
  #else{PNmexp=PNmexp[,match(pron,colnames(PNmexp))]}
  
  if(is.null(PNm))
  {
    PNm<-TFm
    PNm[]<-"N"
    PNm[cbind(1:nrow(PNm),match(obs.v,pron))]="P"
  }else(PNm=as.matrix(PNm))
  #else{PNm=PNm[,match(pron,colnames(PNm))]}
  
  TFm[which(PNm==PNmexp,arr.ind = TRUE)]="T"
  TFm[which(PNm!=PNmexp,arr.ind = TRUE)]="F"
  TPNm=matrix(paste0(TFm,PNm),nrow=nrow(TFm),ncol=ncol(TFm))
  id.TP=which(rowSums(TPNm=="TP")>0)
  id.FN=which(rowSums(TPNm=="FN")>0)
  if(sum(id.FN %in% id.TP)>0)
  {
    TPNm[id.FN[id.FN %in% id.TP],]=gsub("FN","TN",TPNm[id.FN[id.FN %in% id.TP],])
  }
  TPNabs<-function(type="TP",...)
  {
    outx=(sapply(1:ncol(TPNm),function(i){sum(ab.v[which(TPNm[,i]==type)])}))/sum(ab.v)
    c(outx,mean(outx))
  }
  TPab=TPNabs("TP")
  TNab=TPNabs("TN")
  FPab=TPNabs("FP")
  FNab=TPNabs("FN")
  ACC=(TPab+TNab)/(TPab+TNab+FPab+FNab)
  PRC=TPab/(TPab+FPab)
  SPC=TNab/(TNab+FPab)
  SST=TPab/(TPab+FNab)
  #out=rbind(TPab,TNab,FPab,FNab,ACC,PRC,SPC,SST)
  out=rbind(ACC,PRC,SPC,SST)
  colnames(out)=c(pron,"Total")
  out
}