
rm(list=ls(all=TRUE))

###########################
# 1 # file/folder paths/names and package loading
# you may need to change them to the paths on your computer before testing the code.

wd="./Data/EmpiricalData" # the folder to save results
wd.in=wd # the folder saving input files
tool.code.wd="./Code/Tools" # the folder saving tools.r
com.file="OKCW5Y.OTU.csv"  # OTU table
tree.file="OKCW5Y.Tree.nwk" # phylogenetic tree
treat.file="OKCW5Y.Treatment.csv" # treatment information

######
source(paste0(tool.code.wd,"/tools.r"))
library(ape)
library(vegan)
library(bigmemory)
library(NST)
library(iCAMP)

#############################
# 2 # define settings
prefix="OKCW5Y"
nworker=4 # parallel thread number # i used 20 or more on our server. otherwise, it may take too long.

x=1 # test numbering
rand.time=1000 # i changed this to test the impact of different randomization times.
sig.index="SES.RC" # null model significance testing index for iCAMP 
sig.index2="Confidence" # alternative option of null model significance testing index for iCAMP, for comparison

#############################
# 3 # load data

iwd(wd.in)
comm=t(lazyopen(com.file))
tree=lazyopen(tree.file)
treat=lazyopen(treat.file)

sampc=match.name(rn.list=list(comm=comm,treat=treat))
comm=sampc$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
treat=sampc$treat

spc=match.name(cn.list=list(comm=comm),tree.list=list(tree=tree))
comm=spc$comm
tree=spc$tree

#############################
# 4 # phylogenetic distance matrix

pd.wd=paste0(wd,"/pdbig")
if(!dir.exists(pd.wd)){dir.create(pd.wd)}
iwd(pd.wd)
if(file.exists("pd.bin"))
{
  pdbig=list()
  pdbig$tip.label=lazyopen("pd.taxon.name.csv")[,1]
  pdbig$pd.wd=pd.wd
  pdbig$pd.file="pd.desc"
}else{
  pdbig=iCAMP::pdist.big(tree = tree, wd = pd.wd, nworker = nworker)
}


#############################
# 5 # explore community assembly mechanisms with different methods
prefixx=paste0(prefix,"T",x,"R",rand.time)
wd.outx=paste0(wd,"/output.T",x,"R",rand.time)
if(!dir.exists(wd.outx)){dir.create(wd.outx)}
iwd(wd.outx)

# 5.1 # iCAMP
# 5.1.1 # iCAMP based on betaNRI and RC
icout=iCAMP::icamp.big(comm=comm, tree=tree, pd.desc = pdbig$pd.file, pd.spname = pdbig$tip.label,
                       pd.wd = pdbig$pd.wd, rand = rand.time, 
                       prefix = paste0(prefixx,".iCAMP.",sig.index), ds = 0.2, pd.cut = NA, 
                       sp.check = TRUE, phylo.rand.scale = "within.bin",
                       taxa.rand.scale = "across.all", phylo.metric = "bMNTD",
                       sig.index=sig.index, bin.size.limit = 12,
                       nworker = nworker, memory.G = 50, rtree.save = FALSE,
                       detail.save = TRUE, qp.save = TRUE, detail.null=TRUE,
                       ignore.zero = TRUE, output.wd = wd.outx, correct.special = TRUE,
                       unit.sum = NULL, special.method = "depend",
                       ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                       omit.option="no",meta.ab=NULL)
save.file(icout$bNRIiRCa, prefix=prefixx, filename = "iCAMP.bNRIiRCa.sum")
t1=Sys.time()
(t.icp=format(t1-t0))

# 5.1.2 # iCAMP based on Confidence
icout2=iCAMP::change.sigindex(icamp.output = icout, sig.index = sig.index2, detail.save = TRUE)
save(icout2,file = paste0(prefix,".iCAMP.",sig.index2,".detail.rda"))

# 5.1.3 # comparing different sig.index
head(treat)
treat.use=treat[,"yeartreatment",drop=FALSE]
icbin=iCAMP::icamp.bins(icamp.detail = icout,treat = treat.use,between.group = TRUE)
icbin2=iCAMP::icamp.bins(icamp.detail = icout2,treat = treat.use,between.group = TRUE)
Pt1=icbin$Pt
Pt2=icbin2$Pt
if(sum(Pt1$Group!=Pt2$Group)>0){stop("something wrong in Pt.")}
colnames(Pt1)[4:ncol(Pt1)]=paste0("Confidence.",colnames(Pt1)[4:ncol(Pt1)])
colnames(Pt2)[4:ncol(Pt2)]=paste0("SESRC.",colnames(Pt2)[4:ncol(Pt2)])
Pt12=data.frame(Dataset=prefix,Group=Pt1$Group,
                Pt1[,4:ncol(Pt1),drop=FALSE],
                Pt2[,4:ncol(Pt2),drop=FALSE],stringsAsFactors = FALSE)
save.file(Pt12,prefix = prefixx,filename = "iCAMP.Pt.Conf.SESRC")

Ptktrs<-function(icbin,idprefix)
{
  Ptk=icbin$Ptk
  bw=icbin$Binwt
  bwv=as.numeric(as.vector(as.matrix(bw[,4:ncol(bw)])))
  grpbinnb=paste0(rep(bw$Group,times=ncol(bw)-3),"_",rep(colnames(bw)[4:ncol(bw)],each=nrow(bw)))
  
  procn=c("HeS","HoS","DL","HD","DR")
  Ptk1=Ptk[which(Ptk$Index==procn[1]),]
  grpn=rep(Ptk1$Group,times=ncol(Ptk1)-4)
  binn=rep(colnames(Ptk1)[5:ncol(Ptk1)],each=nrow(Ptk1))
  grpbinn=paste0(grpn,"_",binn)
  if(sum(grpbinn!=grpbinnb)>0){stop("wrong in grpbinnb")}
  
  Ptktm=sapply(1:length(procn),
               function(i)
               {
                 Ptki=Ptk[which(Ptk$Index==procn[i]),]
                 grpbinni=paste0(rep(Ptki$Group,times=ncol(Ptki)-4),"_",rep(colnames(Ptki)[5:ncol(Ptki)],each=nrow(Ptki)))
                 if(sum(grpbinn!=grpbinni)>0){stop("wrong in grpbinn")}
                 as.numeric(as.vector(as.matrix(Ptki[,5:ncol(Ptki)])))
               })
  colnames(Ptktm)=paste0(idprefix,".",procn)
  data.frame(Group=grpn,BinID=binn,BinWt=bwv,Ptktm,stringsAsFactors = FALSE)
}
Ptk1=Ptktrs(icbin = icbin,idprefix = "Confidence")
Ptk2=Ptktrs(icbin = icbin2,idprefix = "SESRC")
if(sum(Ptk1[,1:3]!=Ptk2[,1:3])>0){stop("wrong in Ptk12")}
Ptk12=data.frame(Dataset=prefix,Ptk1,Ptk2[,4:ncol(Ptk2)],stringsAsFactors = FALSE)
save.file(Ptk12,prefix = prefixx,filename = "iCAMP.Ptk.Conf.SESRC")

# 5.2 # QPEN
qpout=iCAMP::qpen(comm = comm, pd = pdbig$pd.file, 
                  pd.big.wd = pdbig$pd.wd, pd.big.spname = pdbig$tip.label,
                  tree = tree, bNTI = NULL, RC = NULL, ab.weight = TRUE, 
                  exclude.conspecifics = FALSE, rand.time = rand.time, 
                  sig.bNTI = 1.96, sig.rc = 0.95, nworker = nworker,
                  memory.G = 50, project = paste0(prefixx,".QPEN"), wd = wd.outx,
                  output.detail = TRUE, save.bNTIRC = TRUE)
save(qpout,file = paste0(wd.outx,"/",prefixx,".QPEN.detail.rda"))

# 5.3 # Neutral taxa percentage
prefixi=paste0(prefix,".NP")
treat.use=treat[,"yeartreatment",drop=FALSE]
metagroup=treat[,"year",drop=FALSE]
npout=iCAMP::snm.comm(comm,treat=treat.use,meta.group=metagroup,
                      rand=1000,taxon=NULL,alpha=0.05,two.tail=TRUE,output.detail=TRUE)
save.file(snm.comp$stats,prefix = prefixi, filename = "NeutralModel.statistics")
save.file(snm.comp$plot.detail,prefix=prefixi,filename = "NeutralModel.Details")
save.file(snm.comp$ratio.summary,prefix=prefixi,filename = "NeutralModel.Ratios")
save.file(snm.comp$pvalues,prefix = prefixi,filename = "NeutralModel.comparison")

# 5.4 # tNST and pNST were done in a previous study (Guo et al 2018)
# tNST
prefixj=paste0(prefix,".tNST")
treat.use=treat[,"yeartreatment",drop=FALSE]
tnstout=NST::tNST(comm=comm, group=treat.use,
                  dist.method="bray",abundance.weighted=TRUE,
                  rand=1000,output.rand=FALSE,nworker=nworker,
                  null.model="PF",between.group=TRUE,
                  SES=TRUE,RC=TRUE)
save(tnstout,file = paste0(prefixj,".detail.rda"))
save.file(tnstout$index.grp,prefix=prefixi,filename = "group.summary")
save.file(tnstout$index.pair.grp,prefix=prefixi,filename = "pairwise")

# pNST
prefixj=paste0(prefix,".pNST")
treat.use=treat[,"yeartreatment",drop=FALSE]
pnstout=NST::pNST(comm=comm, pd.desc=pdbig$pd.file,
                  pd.wd=pdbig$pd.wd, pd.spname=pdbig$tip.label,
                  group=treat.use, abundance.weighted=TRUE,
                  rand=1000,output.rand=FALSE, taxo.null.model="PF",
                  phylo.shuffle=TRUE, nworker=nworker,
                  between.group=TRUE, SES=TRUE,RC=TRUE)
save(pnstout,file = paste0(prefixj,".detail.rda"))
save.file(pnstout$index.grp,prefix=prefixi,filename = "group.summary")
save.file(pnstout$index.pair.grp,prefix=prefixi,filename = "pairwise")

# End #