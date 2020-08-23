# retest 2018.10.18
# retest 2018.10.19
# revise 2020.4.24
# 2020.5.2 change to server version
# 2020.8.3 change to use Confidence rather than bNRI or RC.
# GitHub backup 2020.8.20

rm(list=ls(all=TRUE))

###########################
# 1 # file/folder paths and package loading
# you may need to change them to the paths on your computer before testing the code.
code.wd="./Code/Optimization" # a path to the folder you saved the code.
tool.code.wd="./Code/Tools" # the folder saving tools.r
wdm="./Data/SimulatedData/NoComp" # to test with the data without competition
with.competition=FALSE

######
source(paste0(tool.code.wd,"/tools.r"))
library(ape)
library(vegan)
library(bigmemory)
library(iCAMP)

#############################
# 2 # define settings
# give a precode for a certain algorithm/parameter setting. Run different settings and compare the performance summary to optimize the iCAMP settings.
# for example
# precode b #
## phylo.rand.scale="within.bin";taxa.rand.scale="across.all";phylo.metric="bMPD"
## ds=0.2;bin.size=24;sig.index="Confidence"
## correct.special=TRUE;unitsum=TRUE;ses.cut=1.96;rc.cut=0.95;conf.cut=0.975

precode="b" 
phylo.rand.scale="within.bin";taxa.rand.scale="across.all";phylo.metric="bMPD"
ds=0.2;bin.size=24;sig.index="Confidence"
correct.special=TRUE;unitsum=TRUE;ses.cut=1.96;rc.cut=0.95;conf.cut=0.975

######
nworker=4 # parallel thread number
dirpattern="d0" # set some pattern if you want to test the performance with a certian type of simulated situations

#########################
# 3 # basic evalution for each situation

sim.test<-function(sim.file,wd,...)
{
  
  wd=iwd(wd)
  prefix=paste0(precode,sub(".comm.tree.pd.*.rda$","",substring(sim.file,regexpr("\\/[^\\/]*$",sim.file)+1)))
  print(prefix)

  simd=lazyopen(sim.file)
  comm=simd$comm
  pd=simd$pd
  tree=simd$tree
  ABP=simd$ABP
  
  pd.wd=paste0(wd,"/",prefix,".pdbig")
  if(!dir.exists(pd.wd)){dir.create(pd.wd)}
  iwd(pd.wd)
  if(file.exists("pd.bin"))
  {
    pd.big=list()
    pd.big$tip.label=lazyopen("pd.taxon.name.csv")[,1]
    pd.big$pd.wd=pd.wd
    pd.big$pd.file="pd.desc"
    pd.big$pd.name.file="pd.taxon.name.csv"
  }else{
    pdbig = bigmemory::big.matrix(nrow = nrow(pd), ncol = nrow(pd), 
                                  type = "double", backingfile = "pd.bin", descriptorfile = "pd.desc", 
                                  shared = TRUE)
    pdbig[] = as.matrix(pd)
    write.csv(data.frame(taxon.name = rownames(pd)), file = "pd.taxon.name.csv")
    pd.big=list()
    pd.big$tip.label=rownames(pd)
    pd.big$pd.wd=pd.wd
    pd.big$pd.file="pd.desc"
    pd.big$pd.name.file="pd.taxon.name.csv"
  }
  
  if(unitsum){unit.sum=rowSums(comm)}else{unit.sum=NULL}
  
  library(iCAMP)
  icamp.res=icamp.big(comm=comm,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, 
                      pd.wd=pd.big$pd.wd,rand=10,tree=tree,
                      prefix=prefix,ds=ds,pd.cut=NA,sp.check=TRUE,
                      phylo.rand.scale=phylo.rand.scale,taxa.rand.scale=taxa.rand.scale,
                      phylo.metric=phylo.metric,sig.index=sig.index,
                      bin.size.limit=bin.size,nworker=nworker,
                      rtree.save=FALSE,detail.save=TRUE,qp.save=TRUE,detail.null=FALSE,
                      ignore.zero=TRUE,output.wd=wd,
                      correct.special=correct.special,unit.sum=unit.sum,
                      special.method="depend",ses.cut = ses.cut,rc.cut = rc.cut,
                      conf.cut=conf.cut,omit.option="no",meta.ab=NULL)
  source(paste0(code.wd,"/evaluate.sim.r"))
  if(!with.competition)
  {
    iceva=evaluate.sim(icamp.res.detail=icamp.res$detail,ABP=ABP,comm=comm)
  }else{
    iceva=evaluate.sim(icamp.res.detail=icamp.res$detail,ABP=ABP,comm=comm, sel=TRUE)
  }
  save(iceva,file = paste0(wd,"/",prefix,".iCAMP.evaluate.rda"))
  prefix
}

###############
simtestx<-function(prefix,...)
{
  wd=paste0(wdm,"/",prefix)
  wd=iwd(wd)
  sim.files=sort(list.files(pattern = ".comm.tree.pd.*.rda"))
  for(i in 1:length(sim.files))
  {
    sim.test(sim.file=sim.files[i],wd=wd)
  }
  sim.files
}
iwd(wdm)
prefixss=list.dirs(path=wdm,full.names = FALSE,recursive = FALSE)
prefixs=prefixss[grep(pattern = dirpattern,prefixss)]
prefixs
length(prefixs)
########
for(ps in 1:length(prefixs))
{
  simtestx(prefixs[ps])
}

#########################
# 4 # summary the evalution from different situations
iwd(wdm)
prefixs=sort(list.dirs(path=wdm,full.names = FALSE,recursive = FALSE))

pron=c("HoS","HeS","HD","DL","DR")
pron2=c("Sel","HD","DL","DR")
source(paste0(code.wd,"/quali.perf.r"))
source(paste0(code.wd,"/quant.perf.r"))
perf=list()
STevl=list()

for(x in 1:length(prefixs))
{
  prefixi=prefixs[x]  
  message("---Now x=",x," prefix=",prefixi,". ",date())
  Puvml<-Pgkml<-Puvkml<-list()
  
  senai=substr(prefixi,1,2)
  abszi=substr(prefixi,regexpr("b",prefixi)+1,regexpr("d",prefixi)-1)
  if(with.competition)
  {
    dsi=substr(prefixi,regexpr("d",prefixi)+1,regexpr("C",prefixi)-1)
    compri=substring(prefixi,regexpr("C",prefixi)+1)
  }else{
    dsi=substring(prefixi,regexpr("d",prefixi)+1)
    compri=0
  }
  
  wdi=iwd(paste0(wdm,"/",prefixi))
  evsi.fs=list.files(path = wdi,pattern = ".evaluate.rda$")
  
  for(j in 1:length(evsi.fs))
  {
    fsij=evsi.fs[j]
    situj=substr(fsij,regexpr(".S[1-15]",fsij)+1,regexpr(".iCAMP",fsij)-1)
    evsij=lazyopen(fsij)
    Puvml[[j]]=data.frame(Scenario=senai,Abinsize=abszi,ds=dsi,situation=situj,evsij$Puv,stringsAsFactors = FALSE)
    Pgkml[[j]]=data.frame(Scenario=senai,Abinsize=abszi,ds=dsi,situation=situj,evsij$Pgk,stringsAsFactors = FALSE)
    Puvkml[[j]]=data.frame(Scenario=senai,Abinsize=abszi,ds=dsi,situation=situj,evsij$Puvk,stringsAsFactors = FALSE)
  }
  Puvm=Reduce(rbind,Puvml)
  Pgkm=Reduce(rbind,Pgkml)
  Pgkm=Pgkm[which(Pgkm$RA.bin>0),,drop=FALSE]
  Puvkm=Reduce(rbind,Puvkml)
  Puvkm=Puvkm[which(Puvkm$RelativeAbundance>0),,drop=FALSE]
  
  STexp=rowSums(Puvm[,paste0(c("HD","DL","DR"),".PtExpect"),drop=FALSE])
  STest=rowSums(Puvm[,paste0(c("HD","DL","DR"),".PtObs"),drop=FALSE])
  STuv=data.frame(Approach="iCAMP",Scenario=senai,ABinSize=abszi,Act.ds=dsi,CompetitionRatio=compri,Situation=Puvm$situation,
                  Level="community", Plot=paste0(Puvm$name1,Puvm$name2),STexp=STexp,STest=STest,stringsAsFactors = FALSE)
  stev=ieggr::ccc(x=STexp,y=STest,out.vector = TRUE)[1:2]
  STexpk=rowSums(Pgkm[,paste0(c("HD","DL","DR"),".PtExpect"),drop=FALSE])
  STestk=rowSums(Pgkm[,paste0(c("HD","DL","DR"),".PtObs"),drop=FALSE])
  STuvk=data.frame(Approach="iCAMP",Scenario=senai,ABinSize=abszi,Act.ds=dsi,CompetitionRatio=compri,Situation=Pgkm$situation,
                   Level="bin", Plot=Pgkm$plots,STexp=STexpk,STest=STestk,stringsAsFactors = FALSE)
  stevk=ieggr::ccc(x=STexpk,y=STestk,out.vector = TRUE)[1:2]
  stevkw=ieggr::ccc(x=STexpk*Pgkm$RA.bin,y=STestk*Pgkm$RA.bin,out.vector = TRUE)[1:2]
  STevl[[x]]=data.frame(STuvk[1:3,c(2:5)],Level=c("community","bin","bin.wt"),rbind(stev,stevk,stevkw),stringsAsFactors = FALSE)
  
  
  if(with.competition)
  {
    Puv.perf=rbind(quant.perf(expm=Puvm[,paste0(pron2,".PtExpect2")],
                              obsm=Puvm[,paste0(pron2,".PtObs2")],pron = pron2),
                   quali.perf(PNmexp = Puvm[,paste0(pron2,".WtExpect2")],
                              PNm = Puvm[,paste0(pron2,".WtObs2")],
                              ab.v = rep(1,nrow(Puvm)),pron = pron2))
    
    Pk.perf=rbind(quant.perf(expm=Pgkm[,c(paste0(pron2[1],".PtExpect2"),paste0(pron2[2:length(pron2)],".PtExpect"))],
                             obsm=Pgkm[,c(paste0(pron2[1],".PtObs2"),paste0(pron2[2:length(pron2)],".PtObs"))],pron = pron2,
                             RA=Pgkm$RA.bin),
                  quali.perf(PNmexp = Puvkm[,c(paste0(pron2[1],".WtExpect2"),paste0(pron2[2:length(pron2)],".WtExpect"))],
                             obs.v = Puvkm$ObsPuvk2,
                             ab.v = Puvkm$RelativeAbundance,pron = pron2))
    perfi=data.frame(Scenario=senai,Abinsize=abszi,ds=dsi,CompetitionRatio=compri,level=c(rep("comm",nrow(Puv.perf)),rep("bin",nrow(Pk.perf))),
                     index=c(rownames(Puv.perf),rownames(Pk.perf)),rbind(Puv.perf,Pk.perf),stringsAsFactors = FALSE)
  }else{
    Puv.perf=rbind(quant.perf(expm=Puvm[,paste0(pron,".PtExpect")],
                              obsm=Puvm[,paste0(pron,".PtObs")],pron = pron),
                   quali.perf(PNmexp = Puvm[,paste0(pron,".WtExpect")],
                              PNm = Puvm[,paste0(pron,".WtObs")],
                              ab.v = rep(1,nrow(Puvm)),pron = pron))
    Pk.perf=rbind(quant.perf(expm=Pgkm[,paste0(pron,".PtExpect")],
                             obsm=Pgkm[,paste0(pron,".PtObs")],pron = pron,
                             RA=Pgkm$RA.bin),
                  quali.perf(PNmexp = Puvkm[,paste0(pron,".WtExpect")],
                             obs.v = Puvkm$ObsPuvk,
                             ab.v = Puvkm$RelativeAbundance,pron = pron))
    perfi=data.frame(Scenario=senai,Abinsize=abszi,ds=dsi,CompetitionRatio=compri,level=c(rep("comm",nrow(Puv.perf)),rep("bin",nrow(Pk.perf))),
                     index=c(rownames(Puv.perf),rownames(Pk.perf)),rbind(Puv.perf,Pk.perf),stringsAsFactors = FALSE)
  }
  rownames(perfi)=c()
  perf[[x]]=perfi
}
perfm=Reduce(rbind,perf)
head(perfm)
dim(perfm)
time.code=format(Sys.time(),"%y%m%d%H%M")
save.file(perfm,filename = paste0("iCAMP.Perform.Summary.",time.code),folder = wdm)

# End #