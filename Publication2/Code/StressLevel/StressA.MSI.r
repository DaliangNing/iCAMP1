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

simx=SIraw[match(rownames(comm), rownames(SIraw)),,drop=FALSE]
simx[1:4,1:4]    
#save.file(simx,filename = "SI.matched.20220922")
dim(simx)
sum(simx<0) # should be zero

head(simx)
stn=rowSums(simx>1)

msi=apply(simx,1,max)
min(msi)
msi.rk=order(msi,decreasing = FALSE)

treat=data.frame(group=rep(NA,nrow(comm)),stringsAsFactors = FALSE)
rownames(treat)=rownames(comm)

treat[which(rownames(treat) %in% names(msi)[which(msi.rk<=13)]),1]="A1";sum(treat$group=="A1",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>13 & msi.rk<=(13*2))]),1]="A2";sum(treat$group=="A2",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>(13*2) & msi.rk<=(13*3))]),1]="A3";sum(treat$group=="A3",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>(13*3) & msi.rk<=(13*4))]),1]="A4";sum(treat$group=="A4",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>(13*4) & msi.rk<=(13*5))]),1]="A5";sum(treat$group=="A5",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>(13*5) & msi.rk<=(13*6))]),1]="A6";sum(treat$group=="A6",na.rm = TRUE)
treat[which(rownames(treat) %in% names(msi)[which(msi.rk>(13*6))]),1]="A7";sum(treat$group=="A7",na.rm = TRUE)
table(treat[,1])
save.file(treat,filename = "grouping.stressA.msi",folder = wd)

treat.dat=data.frame(treat,MSI=msi[match(rownames(treat),names(msi))])
save.file(treat.dat,filename = "grouping.stressA.msi.data",folder = wd)

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
save.file(output,filename = "grouping.stress01.statistic",folder = wd)

datmsi=data.frame(lnMSI=log(msi),STn1=si1n,STn2=si2n,CTra=ctra,CTrich=ctri,CTabln=log(ctab),group=factor(treat$group,levels = trt.lev[order(msi.trt.mean,decreasing = FALSE)]))

library(ggplot2)
(p1 <- ggplot(datmsi, aes(x=group, y=lnMSI)) + geom_boxplot())
(p2 <- ggplot(datmsi, aes(x=group, y=si1n)) + geom_boxplot())
(p3 <- ggplot(datmsi, aes(x=group, y=si2n)) + geom_boxplot())
(p4 <- ggplot(datmsi, aes(x=group, y=CTrich)) + geom_boxplot())

treat

