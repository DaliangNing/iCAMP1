rm(list=ls())
wd="/Code/StressLevel"
SIraw.f="SI.20220922.csv"
com.f="/Data/100WSc.OTUtable.csv"
env.f="/Data/100WSc.Env.csv"
ct.f="Reference.taxa.ra0.001fq0.5.rda"

source("/Code/handytool.r")
wd=iwd(wd)
envin=lazyopen(env.f)
SIraw=lazyopen(SIraw.f)
comm=t(lazyopen(com.f))
sum(!(rownames(comm) %in% rownames(SIraw)))

simx=SIraw[match(rownames(comm), rownames(SIraw)),,drop=FALSE]
simx[,1:4]    
wdst=iwd(SIraw.f)

dim(simx)
sum(simx<0) # should be zero

head(simx)
stn=rowSums(simx>1)
table(stn)
treat=data.frame(group=rep("B1",nrow(comm)),stringsAsFactors = FALSE)
rownames(treat)=rownames(comm)

treat[which(rownames(treat) %in% names(stn)[which(stn==0)]),1]="B1"
treat[which(rownames(treat) %in% names(stn)[which(stn==1)]),1]="B2"
treat[which(rownames(treat) %in% names(stn)[which(stn==2)]),1]="B3"
treat[which(rownames(treat) %in% names(stn)[which(stn %in% c(3,4))]),1]="B4"
treat[which(rownames(treat) %in% names(stn)[which(stn %in% 5:9)]),1]="B5"
treat[which(rownames(treat) %in% names(stn)[which(stn>=11)]),1]="B6"
table(treat[,1])

save.file(treat,filename = "grouping.stressB.stressornumber",folder = wd)

treat.dat=data.frame(treat,STN=stn[match(rownames(treat),names(stn))])
save.file(treat.dat,filename = "grouping.stressB.stressornumber.data",folder = wd)



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

trt.lev=sort(unique(treat[,1]))

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
save.file(output,filename = "grouping.stressB.statistic",folder = wd)

datmsi=data.frame(lnMSI=log(msi),STn1=si1n,STn2=si2n,CTra=ctra,CTrich=ctri,CTabln=log(ctab),group=factor(treat$group,levels = trt.lev[order(msi.trt.mean,decreasing = FALSE)]))

library(ggplot2)
(p1 <- ggplot(datmsi, aes(x=group, y=lnMSI)) + geom_boxplot())
(p2 <- ggplot(datmsi, aes(x=group, y=si1n)) + geom_boxplot())
(p3 <- ggplot(datmsi, aes(x=group, y=si2n)) + geom_boxplot())
(p4 <- ggplot(datmsi, aes(x=group, y=CTrich)) + geom_boxplot())

treat

