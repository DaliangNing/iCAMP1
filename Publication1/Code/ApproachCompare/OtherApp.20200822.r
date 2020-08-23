# retest 2018.10.18
# retest 2018.10.19
# revise 2020.4.24
# 2020.5.2 change to server version
# 2020.5.2 other approaches
# for GitHub version, 2020.8.22

rm(list=ls(all=TRUE))

###########################
# 1 # file/folder paths
# you may need to change them to the paths on your computer before testing the code.
code.wd="./Code/Optimization" # the folder saving quant.perf.r and quali.perf.r
wdmm="./Data/SimulatedData" # the folder saving some basic information of simulated data
wdm="./Data/SimulatedData/NoComp" # the folder saving simulated communities
# if test with-compeition datasets, change to the folder FiltComp, like "./Data/SimulatedData/FiltComp".
tool.code.wd="./Code/Tools" # the folder saving tools.r

######
source(paste0(tool.code.wd,"/tools.r"))
library(ape)
library(vegan)
library(bigmemory)
library(NST)
library(iCAMP)


###########################
# 2 # parameter settings.

nworker=4 # parallel thread number 
dirpattern="d0" # set some pattern if you want to test the performance with a certian type of simulated situations 

iwd(wdmm)
Vexp=lazyopen("expect.process.csv",row.names=NULL)

iwd(wdm)
precode="b" # a precode for certain setting of the parameters

pron=c("HoS","HeS","HD","DL","DR")
pron2=c("Sel","HD","DL","DR")

###########################
# 3 # functions.
crect.plotn<-function(v)
{
  crctgrp=matrix(c("LA","LALA",
                   "LB","LBLB",
                   "HA","HAHA",
                   "HB","HBHB",
                   "LBLA","LALB",
                   "HALA","LAHA",
                   "HBLA","LAHB",
                   "HALB","LBHA",
                   "HBLB","LBHB",
                   "HBHA","HAHB"),nrow = 10,ncol=2,byrow = TRUE)
  id.er=which(v %in% crctgrp[,1])
  v[id.er]=crctgrp[match(v[id.er],crctgrp[,1]),2]
  v
}
crct.pron<-function(v)
{
  pronct=matrix(c("Vari_selection","Homo_selection","Dispersal_limit","Dispersal_homogen","Undominated",
                  "Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Homogenizing.Dispersal","Drift.and.Others",
                  "HeS","HoS","DL","HD","DR",
                  "HeS","HoS","DL","HD","DR"),ncol=2)
  id.er=which(v %in% pronct[,1])
  v[id.er]=pronct[match(v[id.er],pronct[,1]),2]
  v
}


sim.test<-function(sim.file,wd,...)
{
  #############
  ses.cut=1.96
  rc.cut=0.95
  #######################
  wd=iwd(wd)
  prefix=paste0(precode,sub(".comm.tree.pd.*.rda$","",substring(sim.file,regexpr("\\/[^\\/]*$",sim.file)+1)))
  print(prefix)
  sena=substr(prefix,2,3)
  abinsz=substr(prefix,regexpr("b[0-9]",prefix)+1,regexpr("d[0-9]",prefix)-1)
  ads=substr(prefix,regexpr("d[0-9]",prefix),regexpr("d[0-9]",prefix)+2)
  situ=substring(prefix,regexpr("\\.S[0-9]",prefix)+1)
  if(grepl("d[0-9].*C[0-9]",prefix)){compr=substr(prefix,regexpr("C[0-9]",prefix)+1,regexpr("\\.S[0-9]",prefix)-1)}else{compr=0}

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
  
  iwd(wd)
  treat=data.frame(plot=substr(rownames(comm),1,2),stringsAsFactors = FALSE)
  rownames(treat)=rownames(comm)
  
  # NP
  snm.comp=iCAMP::snm.comm(comm = comm, treat = treat, meta.group = NULL, rand.time = 1000, taxon = NULL)
  save(snm.comp,file=paste0(prefix,".NeutralModel.rda"))
  idu=which(snm.comp$ratio.summary$index=="obs")
  snmres=data.frame(plot=as.vector(snm.comp$ratio.summary$treatment.id[idu]),
                    STest=snm.comp$ratio.summary$Neutral.wt[idu],stringsAsFactors = FALSE)
  snmres[,1]=crect.plotn(snmres[,1])
  
  # tNST
  tnst=NST::tNST(comm=comm, group=treat, meta.group = NULL, meta.com = NULL, dist.method = "bray", 
                 abundance.weighted = TRUE, rand = 1000, output.rand = FALSE, 
                 nworker = nworker, LB = FALSE, null.model = "PF", between.group = TRUE, 
                 SES = TRUE, RC = TRUE)
  save(tnst,file = paste0(prefix,".tNST.rda"))
  tnstres=rbind(tnst$index.grp[,c("group","NST.i.bray"),drop=FALSE],
                tnst$index.pair.grp[,c("group","NST.i.bray"),drop=FALSE])
  tnstres[,1]=crect.plotn(gsub(".vs.","",tnstres[,1]))
  colnames(tnstres)=c("plot","STest")
  rownames(tnstres)=c()
  
  RC3c=tnst$index[,c("name1","name2","RC.bray"),drop=FALSE]
  RC=col3.dist(m.3col = RC3c,to.dist = FALSE)
  
  # pNST
  pnst=NST::pNST(comm=comm, pd.desc = pd.big$pd.file, pd.wd = pd.big$pd.wd,
                 pd.spname = pd.big$tip.label, group=treat, meta.group = NULL,
                 abundance.weighted = TRUE, rand = 1000, output.rand = FALSE,
                 taxo.null.model = NULL, phylo.shuffle = TRUE, 
                 exclude.conspecifics = FALSE, nworker = nworker, LB = FALSE,
                 between.group = TRUE, SES = TRUE, RC = FALSE)
  save(pnst,file = paste0(prefix,".pNST.rda"))
  
  pnstres=rbind(pnst$index.grp[,c("group","NST.i.dis"),drop=FALSE],
                pnst$index.pair.grp[,c("group","NST.i.dis"),drop=FALSE])
  pnstres[,1]=crect.plotn(gsub(".vs.","",pnstres[,1]))
  colnames(pnstres)=c("plot","STest")
  rownames(pnstres)=c()
  bNTI3c=pnst$index[,c("name1","name2","bNTI.wt"),drop=FALSE]
  bNTI=col3.dist(m.3col = bNTI3c,to.dist = FALSE)
  
  # QPEN
  spc=iCAMP::match.name(both.list = list(RC=RC,bNTI=bNTI))
  RC=spc$RC
  bNTI=spc$bNTI
  
  qpeno=iCAMP::qpen(comm = comm, bNTI = bNTI, RC = RC, ab.weight = TRUE,  
                   rand.time = 1000, sig.bNTI = ses.cut, sig.rc = rc.cut,
                   nworker = nworker)
  qpen=qpeno$result
  grp2=paste0(substr(qpen$sample1,1,2),substr(qpen$sample2,1,2))
  grp2=crect.plotn(grp2)
  grp2.lev=sort(unique(grp2))
  
  qpen$process=crct.pron(qpen$process)
  p.ri<-function(v)
  {
    pron=c("HoS","HeS","HD","DL","DR")
    sapply(pron,function(pn){sum(v==pn)})/length(v)
  }
  qpres=t(sapply(1:length(grp2.lev),
                 function(i)
                 {
                   idi=which(grp2==grp2.lev[i])
                   p.ri(qpen$process[idi])
                 }))
  colnames(qpres)=paste0(colnames(qpres),"obs")
  qpst=data.frame(plot=grp2.lev,STest=rowSums(qpres[,3:5]),stringsAsFactors = FALSE)
  qpres2=data.frame(Selobs=rowSums(qpres[,1:2]),qpres[,3:5,drop=FALSE],stringsAsFactors = FALSE)
  
  Vexpi=Vexp[which(Vexp$Situ==substr(situ,1,regexpr("\\.",situ)-1)),,drop=FALSE]
  
  STrescb=rbind(snmres,tnstres,pnstres,qpst)
  STres=data.frame(Approach=c(rep("NP",nrow(snmres)),rep("tNST",nrow(tnstres)),rep("pNST",nrow(pnstres)),rep("QPEN",nrow(qpst))),
                   Scenario=sena,ABinSize=abinsz,Act.ds=ads,CompetitionRatio=compr,Situation=situ,
                   STrescb[,1,drop=FALSE],STexp=Vexpi$STexp[match(STrescb$plot,Vexpi$Group)],
                   STrescb[,2,drop=FALSE],stringsAsFactors = FALSE)
  QPres=data.frame(Approach="QPEN",Scenario=sena,ABinSize=abinsz,Act.ds=ads,CompetitionRatio=compr,Situation=situ,
                   plot=grp2.lev,Vexpi[match(grp2.lev,Vexpi$Group),paste0(pron,"exp"),drop=FALSE],
                   qpres[,paste0(pron,"obs"),drop=FALSE],stringsAsFactors = FALSE)
  QPres2=data.frame(Approach="QPEN",Scenario=sena,ABinSize=abinsz,Act.ds=ads,CompetitionRatio=compr,Situation=situ,
                   plot=grp2.lev,Vexpi[match(grp2.lev,Vexpi$Group),paste0(pron2,"exp"),drop=FALSE],
                   qpres2[,paste0(pron2,"obs"),drop=FALSE],stringsAsFactors = FALSE)
  
  out=list(STres=STres,QPres=QPres,QPres2=QPres2)
  save(out,file = paste0(wd,"/",prefix,".OtherApp.rda"))
  out
}

###############
simtestx<-function(prefix,...)
{
  wd=paste0(wdm,"/",prefix)
  wd=iwd(wd)
  sim.files=sort(list.files(pattern = ".comm.tree.pd.*.rda"))
  tests=list()
  for(i in 1:length(sim.files))
  {
    tests[[i]]=sim.test(sim.file=sim.files[i],wd=wd)
  }
  names(tests)=sim.files
  save(tests,file=paste0(wd,"/",prefix,".OtherAppAll.rda"))
  stevlu<-function(app,...)
  {
    asl=lapply(1:length(tests),function(i){tests[[i]]$STres[which(tests[[i]]$STres$Approach==app),,drop=FALSE]})
    asm=Reduce(rbind,asl)
    cccx=ccc(x=asm$STexp,y=asm$STest)
    c(qACC=cccx$accuracy,qPRC=cccx$precision)
  }
  app.lev=unique(tests[[1]]$STres$Approach)
  STevalu=data.frame(tests[[1]]$STres[1:length(app.lev),2:5,drop=FALSE],Approach=app.lev,t(sapply(app.lev,stevlu)),stringsAsFactors = FALSE)
  
  
  qpm=Reduce(rbind,lapply(1:length(tests),function(i){tests[[i]]$QPres}))
  qpm2=Reduce(rbind,lapply(1:length(tests),function(i){tests[[i]]$QPres2}))
  qpevlu<-function(Pexpm,Pobsm,pronu,...)
  {
    source(paste0(code.wd,"/quant.perf.r"))
    source(paste0(code.wd,"/quali.perf.r"))
    qtp=quant.perf(expm = Pexpm,obsm = Pobsm,pron = pronu)
    PNmatrix<-function(Pm)
    {
      prmax=apply(Pm,1,max)
      pnm=Pm
      pnm[]="N"
      pnm[Pm==prmax]="P"
      pnm
    }
    qlp=quali.perf(PNmexp = PNmatrix(Pexpm),PNm = PNmatrix(Pobsm),
                   ab.v = rep(1,nrow(Pexpm)),pron = pronu)
    rbind(qtp,qlp)
  }
  qpe=qpevlu(Pexpm = qpm[,paste0(pron,"exp"),drop=FALSE],
             Pobsm = qpm[,paste0(pron,"obs"),drop=FALSE],
             pronu = pron)
  qpe2=qpevlu(Pexpm = qpm2[,paste0(pron2,"exp"),drop=FALSE],
             Pobsm = qpm2[,paste0(pron2,"obs"),drop=FALSE],
             pronu = pron2)
  colnames(qpe2)=paste0(colnames(qpe2),"2")
  qpem=data.frame(tests[[1]]$STres[1:6,2:5,drop=FALSE],Approach="QPEN",qpe,qpe2,stringsAsFactors = FALSE)
  list(STevalu=STevalu,QPevalu=qpem)
}

###########################
# 4 # calculation
iwd(wdm)
prefixss=list.dirs(path=wdm,full.names = FALSE,recursive = FALSE)
prefixs=prefixss[grep(pattern = dirpattern,prefixss)]
prefixs
length(prefixs)
#############################
simsum=list()
for(ps in 1:length(prefixs))
{
  simsum[[ps]]=simtestx(prefixs[ps])
}
STsum=Reduce(rbind,lapply(1:length(simsum),function(i){simsum[[i]]$STevalu}))
QPsum=Reduce(rbind,lapply(1:length(simsum),function(i){simsum[[i]]$QPevalu}))
rownames(STsum)=c()
rownames(QPsum)=c()
time.code=format(Sys.time(),"%y%m%d%H%M")
save.file(STsum,filename = paste0("ST.OtherAppPerform.",time.code),folder = wdm)
save.file(QPsum,filename = paste0("QPEN.Perform",time.code),folder = wdm)

# End #