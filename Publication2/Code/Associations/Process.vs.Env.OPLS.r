##
datf="file:///E:/Dropbox/2_1-100wells_Daliang/Writing/Paper1_MSI.crct/2-CommunityAssembly/iCAMP/Envlog2PIcrl.ProcEnvGeoScaled.Dist.csv"
prefix="Envlog2PIcrl"

##
datf="file:///E:/Dropbox/2_1-100wells_Daliang/Writing/Paper1_MSI.crct/2-CommunityAssembly/iCAMP/EnvlogPIcrl.ProcEnvGeoScaled.Dist.csv"
prefix="EnvlogPIcrl"

rm(list=ls())
wd="/Code/Associations"
datf="/Code/Associations/EnvlogPIcrl.ProcEnvGeoScaled.Dist.csv"
prefix="EnvlogPIcrl"
rand=1000 # randomization time


###################
wd=iwd(wd)
dat=lazyopen(datf)
dim(dat)
dat[1:3,1:10]
ym=as.matrix(dat[,3:7,drop=FALSE])
xm=as.matrix(dat[,8:ncol(dat),drop=FALSE])

library(ropls)
modsum=list()
vip=list()
for(i in 1:5)
{
  message("--- i=",i,"  ",date())
  colnames(ym)[i]
  t0=Sys.time()
  tpls <- opls(x=xm,y=ym[,i,drop=FALSE],fig.pdfC='none',predI=NA,orthoI=NA,permI = 1000)
  (t1=format(Sys.time()-t0))
  save(tpls,file = paste0(prefix,".OPLS.",colnames(ym)[i],".rda"))
  
  coefs=coef(tpls)
  head(coefs)
  coefs[order(abs(coefs[,1]),decreasing = TRUE),,drop=FALSE][1:10,,drop=FALSE]
  
  #samp.scores=getScoreMN(tpls)
  wt=getWeightMN(tpls)
  wt[order(abs(wt[,1]),decreasing = TRUE),,drop=FALSE][1:10,,drop=FALSE]
  
  modsum[[i]]=getSummaryDF(tpls)
  modsum[[i]]
  
  vip[[i]]=getVipVn(tpls)
  sort(vip[[i]],decreasing = TRUE)[1:10]
}
modsumm=Reduce(rbind,modsum)
vipm=Reduce(rbind,vip)
sumout=t(cbind(modsumm,vipm))
colnames(sumout)=colnames(ym)
head(sumout)
save.file(sumout,prefix = prefix,filename = "OPLSsummary")

##############################
datm1=col3.dist(dat[,1:3,drop=FALSE],to.dist = FALSE)
dim(datm1)
perm <- vegan:::getPermuteMatrix(rand, nrow(datm1))

samps=sort(unique(as.vector(as.matrix(dat[,1:2]))));length(samps)
oldid1=dat[,1]
oldid2=dat[,2]
oldid12=paste0(oldid1,"_",oldid2)

idrd=lapply(1:nrow(perm),
            function(i)
            {
              permi=perm[i,]
              newid1=samps[permi[match(oldid1,samps)]]
              newid2=samps[permi[match(oldid2,samps)]]
              newid12=paste0(newid1,"_",newid2)
              newid12b=paste0(newid2,"_",newid1)
              idrand=match(newid12,oldid12)
              idrdb=match(newid12b,oldid12)
              idrand[which(is.na(idrand))]=idrdb[which(is.na(idrand))]
              if(sum(is.na(idrand))>0){stop("something wrong!")}
              idrand
            })

oplseach<-function(k,xm,ym)
{
  na.count=integer(0)
  for(i in 1:5)
  {
    message("---k=",k," i=",i,"  ",date())
    tpls <- try(opls(x=xm,y=ym[,i,drop=FALSE],fig.pdfC='none',predI=NA,orthoI=NA))
    if(class(tpls)=="try-error")
    {
      modsum[[i]]=NA
      vip[[i]]=NA
      na.count=c(na.count,i)
    }else{
      modsum[[i]]=getSummaryDF(tpls)
      vip[[i]]=getVipVn(tpls)
    }
  }
  if(length(na.count)==5)
  {
    sumout=NA
  }else{
    nonaid=setdiff((1:5),na.count)[1]
    for(i in na.count)
    {
      modsum[[i]]=modsum[[nonaid]]
      modsum[[i]][]=NA
      vip[[i]]=vip[[nonaid]]
      vip[[i]][]=NA
    }
    modsumm=Reduce(rbind,modsum)
    vipm=Reduce(rbind,vip)
    sumout=t(cbind(modsumm,vipm))
    colnames(sumout)=colnames(ym)
  }
  sumout
}

arch.point=seq(from=1, to=length(idrd), by=20)
oplsr=list()
for(k in 1:length(idrd))
{
  ymr=ym[idrd[[k]],,drop=FALSE]
  rownames(ymr)=rownames(ym)
  oplsr[[k]]=oplseach(k,xm,ymr)
  if(k %in% arch.point)
  {
    message("---- Now archiving k=",k,"  ",date())
    save(oplsr,file = paste0(prefix,".OPLS.permute.rda"))
  }
}
save(oplsr,file = paste0(prefix,".OPLS.permute.rda"))

oplsra=array(unlist(oplsr),dim = c(nrow(oplsr[[1]]),ncol(oplsr[[2]]),length(oplsr)))
opls.obs=array(sumout,dim = c(nrow(oplsr[[1]]),ncol(oplsr[[2]]),length(oplsr)))
p1m=apply(oplsra>opls.obs,c(1,2),sum,na.rm=TRUE)
rownames(p1m)=rownames(sumout)
colnames(p1m)=colnames(sumout)
na1m=apply(is.na(oplsra),c(1,2),sum,na.rm=TRUE)
head(na1m)
rownames(na1m)=rownames(sumout)
colnames(na1m)=paste0("NAcount.",colnames(sumout))

save.file(cbind(p1m,na1m),prefix=prefix,filename = paste0(".OPLS.permtest.pcount",rand))

