rm(list=ls())
wd="/Code/StressLevel"
com.f="/Data/100WSc.OTUtable.csv"
env.f="/Data/100WSc.Env.csv"
SIraw.f="SI.20220922.csv"
ct.f="Reference.taxa.ra0.001fq0.5.rda"

source("/Code/handytool.r")
wd=iwd(wd)
envin=lazyopen(env.f)
SIraw=lazyopen(SIraw.f)
comm=t(lazyopen(com.f))
sum(!(rownames(comm) %in% rownames(SIraw)))

simx=SIraw[match(rownames(comm), rownames(SIraw)),,drop=FALSE]
simx[1:4,1:4]    

dim(simx)
sum(simx<0) # should be zero

head(simx)
stn=rowSums(simx>1)

ctsp=lazyopen(ct.f)
save.file(data.frame(ctsp),filename = "RefTaxa")
comct=comm[,which(colnames(comm) %in% ctsp),drop=FALSE]
ctra=(rowSums(comct)/rowSums(comm))
ctab=ctra*(envin$AODC[match(rownames(comm),rownames(envin))])
ctri=rowSums(comct>0)

table(ctri)
treat=data.frame(group=rep(NA,nrow(comm)),stringsAsFactors = FALSE)
rownames(treat)=rownames(comm)
SRTR=(max(ctri)-ctri)/max(ctri)

treat[which(rownames(treat) %in% names(ctri)[which(ctri>0 & ctri<=33)]),1]="C1";sum(treat$group=="C1",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>33 & ctri<=41)]),1]="C2";sum(treat$group=="C2",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>41 & ctri<=50)]),1]="C3";sum(treat$group=="C3",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>50 & ctri<=54)]),1]="C4";sum(treat$group=="C4",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>54 & ctri<=59)]),1]="C5";sum(treat$group=="C5",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>59 & ctri<=64)]),1]="C6";sum(treat$group=="C6",na.rm = TRUE)
treat[which(rownames(treat) %in% names(ctri)[which(ctri>64)]),1]="C7";sum(treat$group=="C7",na.rm = TRUE)
table(treat[,1])

save.file(treat,filename = "grouping.stressC.SRTR",folder = wd)

srrti=(max(ctri)-ctri)/max(ctri)
treat.dat=data.frame(treat,RTR=ctri[match(rownames(treat),names(ctri))],SRTR=srrti[match(rownames(treat),names(srrti))])
save.file(treat.dat,filename = "grouping.stressC.SRTR.data",folder = wd)


msi=apply(simx,1,max)
min(msi)
si1n=sapply(1:nrow(simx),function(i){sum(simx[i,]>1)})
si2n=sapply(1:nrow(simx),function(i){sum(simx[i,]>2)})
names(si1n)<-names(si2n)<-rownames(simx)



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
save.file(output,filename = "grouping.stressC.statistic",folder = wd)

datmsi=data.frame(lnMSI=log(msi),STn1=si1n,STn2=si2n,CTra=ctra,CTrich=ctri,CTabln=log(ctab),group=factor(treat$group,levels = trt.lev[order(msi.trt.mean,decreasing = FALSE)]))

library(ggplot2)
(p1 <- ggplot(datmsi, aes(x=group, y=lnMSI)) + geom_boxplot())
(p2 <- ggplot(datmsi, aes(x=group, y=si1n)) + geom_boxplot())
(p3 <- ggplot(datmsi, aes(x=group, y=si2n)) + geom_boxplot())
(p4 <- ggplot(datmsi, aes(x=group, y=CTrich)) + geom_boxplot())

treat

