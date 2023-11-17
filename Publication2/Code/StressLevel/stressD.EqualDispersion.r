rm(list=ls())
wd="/Code/StressLevel"
com.f="/Data/100WSc.OTUtable.csv"
env.f="/Data/100WSc.Env.csv"
ct.f="Reference.taxa.ra0.001fq0.5.rda"
SIraw.f="SI.20220922.csv"

source("/Code/handytool.r")
wd=iwd(wd)
envin=lazyopen(env.f)
SIraw=lazyopen(SIraw.f)
comm=t(lazyopen(com.f))
sum(!(rownames(comm) %in% rownames(SIraw)))

colnames(envin)
colnames(SIraw)

simx=SIraw[match(rownames(comm), rownames(SIraw)),,drop=FALSE]
simx[,1:4]    
wdst=iwd(SIraw.f)

dim(simx)
sum(simx<0) # should be zero
id0=unique(which(simx==0,arr.ind = TRUE)[,2])
colnames(simx)[id0]
i=id0[length(id0)]
simx[which(simx[,i]==0),i,drop=FALSE]


lnx<-function(v){out=log(v);if(sum(v==0)>0){out[v==0]=min(density(log(v[v>0]))$x)};out}
simxln=apply(simx,2,lnx)

silnd=dist(simxln)
silndm=as.matrix(silnd)

################################
hc.siln=hclust(silnd)
plot(hc.siln)
ct.siln=cutree(hc.siln,k=9)
table(ct.siln)

treat=data.frame(group=paste0("S",ct.siln))
rownames(treat)=names(ct.siln)
(trt.lev=sort(unique(treat[,1])))

head(treat)


################
# exploring each cluster
##############
samps=list()
i=1 # test i=1 to i=9
if(i==5 | i==9)
{
  sampi=rownames(treat)[which(treat[,1]==trt.lev[i])]
}else if(i==6){
  sampi=rownames(treat)[which(treat[,1]==trt.lev[6] & (!(rownames(treat) %in% c("GW101.11.13.12.02","DP16D.67.11.27.12.02"))))]
}else{
  sampi=rownames(treat)[sample(which(treat[,1]==trt.lev[i]),8)]
}
length(sampi)
samps[[i]]=sampi
mean(as.dist(silndm[samps[[i]],samps[[i]]]))
sd(as.dist(silndm[sampi,sampi]))

################
# after exploring, use 5 levels
##############

trt.grps=list(S01=1,S02=c(2,3),S03=c(5,6),S04=c(4,7),S05=c(8,9))
table(ct.siln)

samps=list()
mds=list()

################
i=5 # test i=1 to i=5
#####
mdi=15
mdsi=list()
k=1
while((mdi>12.4 | mdi<11) & k<200)
{
  if(i==1)
  {
    sampi=rownames(treat)[sample(which(treat[,1]==trt.lev[i]),8)]
  }else{
    idi1=which(treat[,1]==trt.lev[trt.grps[[i]][1]])
    idi2=which(treat[,1]==trt.lev[trt.grps[[i]][2]])
    idi3=c(sample(idi1,1),sample(idi2,1))
    idi4=c(idi1,idi2)
    idi5=idi4[!(idi4 %in% idi3)]
    sampi=rownames(treat)[c(idi3,sample(idi5,6))]
  }
  length(sampi)
  samps[[i]]=sampi
  mdi=mean(as.dist(silndm[samps[[i]],samps[[i]]]))
  mdsi[[k]]=mdi
  k=k+1
  message("mdi=",mdi," k=",k,". ",date())
}
(mds[[i]]=unlist(mdsi))
mean(as.dist(silndm[samps[[i]],samps[[i]]]))
################

sapply(samps,length)

#####################
treat.old=treat
comm.old=comm
simx.old=simx
#####################
timecode=format(Sys.time(),format = "%Y%m%d%H%M%S")
treat=data.frame(group=unlist(lapply(1:length(samps),function(i){outi=rep(paste0("D",i),length(samps[[i]]));names(outi)=samps[[i]];outi})))
save.file(treat,filename = paste0("grouping.stressD.",timecode),folder = wd)
silnds=as.dist(silndm[rownames(treat),rownames(treat)])
(trt.lev=unique(treat[,1]))
library(vegan)
btdisp=betadisper(d=silnds, group=treat[,1], type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
permutest(btdisp, pairwise = FALSE,permutations = 999)
permutest(btdisp, pairwise = TRUE,permutations = 999)

########################
comm=comm.old[match(rownames(treat),rownames(comm.old)),,drop=FALSE]
comm=comm[,colSums(comm)>0,drop=FALSE]
simx=simx.old[match(rownames(treat),rownames(simx.old)),,drop=FALSE]
stn=rowSums(simx>1)
msi=apply(simx,1,max)
min(msi)
msi.rk=order(msi,decreasing = FALSE)

###########################

msi=apply(simx,1,max)
min(msi)
si1n=sapply(1:nrow(simx),function(i){sum(simx[i,]>1)})
si2n=sapply(1:nrow(simx),function(i){sum(simx[i,]>2)})
names(si1n)<-names(si2n)<-rownames(simx)


ctsp=lazyopen(ct.f)
comct=comm[,which(colnames(comm) %in% ctsp),drop=FALSE]
ctra=(rowSums(comct)/rowSums(comm))
ctab=ctra*(envin$AODC[match(rownames(comm),rownames(envin))])
ctri=rowSums(comct>0)


idck=match.name(rn.list = list(treat),v.list = list(msi,si1n,si2n,ctra,ctri,ctab))
# should match

msi.trt.median=sapply(1:length(trt.lev),function(i){median(log(msi[which(names(msi) %in% rownames(treat)[treat$group==trt.lev[i]])]))})
msi.trt.mean=sapply(1:length(trt.lev),function(i){mean(log(msi[which(names(msi) %in% rownames(treat)[treat$group==trt.lev[i]])]))})

si1n.trt.median=sapply(1:length(trt.lev),function(i){median(si1n[which(names(si1n) %in% rownames(treat)[treat$group==trt.lev[i]])])})
si1n.trt.mean=sapply(1:length(trt.lev),function(i){mean(si1n[which(names(si1n) %in% rownames(treat)[treat$group==trt.lev[i]])])})

ctri.trt.median=sapply(1:length(trt.lev),function(i){median(ctri[which(names(ctri) %in% rownames(treat)[treat$group==trt.lev[i]])])})
ctri.trt.mean=sapply(1:length(trt.lev),function(i){mean(ctri[which(names(ctri) %in% rownames(treat)[treat$group==trt.lev[i]])])})

output=data.frame(MSI.median=msi.trt.median,MSI.mean=msi.trt.mean,
           STnum.median=si1n.trt.median,STnum.mean=si1n.trt.mean,
           CTrich.median=ctri.trt.median,CTrich.mean=ctri.trt.mean)
rownames(output)=trt.lev
output
save.file(output,filename = paste0("grouping.stressD.statistic.",timecode),folder = wd)

datmsi=data.frame(lnMSI=log(msi),STn1=si1n,STn2=si2n,CTra=ctra,CTrich=ctri,CTabln=log(ctab),group=factor(treat$group,levels = trt.lev[order(msi.trt.mean,decreasing = FALSE)]))

library(ggplot2)
(p1 <- ggplot(datmsi, aes(x=group, y=lnMSI)) + geom_boxplot())
(p2 <- ggplot(datmsi, aes(x=group, y=si1n)) + geom_boxplot())
(p3 <- ggplot(datmsi, aes(x=group, y=si2n)) + geom_boxplot())
(p4 <- ggplot(datmsi, aes(x=group, y=CTrich)) + geom_boxplot())

treat

