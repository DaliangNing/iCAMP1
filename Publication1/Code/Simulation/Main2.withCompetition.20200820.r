# model 2017.11.15
# retest 2018.10.15
# retest 2018.10.17
# retest 2018.11.7  # bin=6
# refine 2020.4.24
# add competition 2020.4.27
# revise to server version 2020.4.30
# revise to GitHub version 2020.8.20
rm(list=ls())

###########################
# 0 # file/folder paths/names and package loading
# you may need to change them to the paths on your computer before testing the code.

wdmm="./Data/SimulatedData/" # a path to the folder you want to save the results.
wdinput=paste0(wdmm,"/SpeciesPool") # define subfolder which has regional pool
code.wd="./Code/Simulation" # a path to the folder you saved the code.
tool.code.wd="./Code/Tools" 
save.wdm=paste0(wdmm,"/FiltComp") 
if(dir.exists(save.wdm)){dir.create(save.wdm)}

###
setwd(tool.code.wd)
source("tools.r")
library(ape)
library(vegan)
library(bigmemory)
library(iCAMP)

###########################
# 1 # load tree and trait(OPEN) data
#######################################
wdinput=iwd(wdinput)
tp.file="JS.tree.pd.rda"
tp.in=lazyopen(tp.file)
tree=tp.in$tree
pd=tp.in$pd
tree;max(tree$edge.length)
maxpd=max(pd)
op.file="JS.opens.rda"
ops=lazyopen(op.file)

simmu<-function(prefix,op,bin.size,compet.sel,...)
{
  save.wd=paste0(save.wdm,"/",prefix)
  if(!dir.exists(save.wd)){dir.create(save.wd)}
  
  # setting gobal parameters
  sampname=c(paste0("LA",1:6),paste0("LB",1:6),paste0("HA",1:6),paste0("HB",1:6))
  Nk=20000 # individual number in each community
  J=Nk*1000
  theta=5000
  nworker=4
  spnumi=100
  d.cut=(maxpd-ds)/2
  
  #correct.special=TRUE
  #unitsum=TRUE
  ######################
  if(file.exists(paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda")))
  {
    MZSM.meta=lazyopen(paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda"))
  }else{
    MZSM.meta=sads::rmzsm(n=length(op),J=Nk*1000,theta=theta)
    save(MZSM.meta,file=paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda"))
  }
  
  source(paste0(code.wd,"/simulation.functions.r"))
  save.sim<-function(prefixi,comm,pd,tree,save.wd,ABP=NULL)
  {
    #library(mCAM)
    if(!is.null(ABP))
    {
      spc=iCAMP::match.name(cn.list=list(comm=comm),rn.list = list(ABP=ABP),
                           both.list = list(pd=pd),tree.list=list(tree=tree),
                           silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      comm.use=spc$comm
      ABP.use=spc$ABP
      sim.data=list(comm=comm.use,pd=pd.use,tree=tree.use,ABP=ABP.use)
      filen=paste0(save.wd,"/",prefixi,".comm.tree.pd.ABP.rda")
    }else{
      spc=mCAM::match.name(cn.list=list(comm=comm),both.list = list(pd=pd),tree.list=list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      comm.use=spc$comm
      sim.data=list(comm=comm.use,pd=pd.use,tree=tree.use)
      filen=paste0(save.wd,"/",prefixi,".comm.tree.pd.rda")
    }
    save(sim.data,file = filen)
    filen
  }
  
  ABPcount<-function(comm,sel.spname=NULL,flt.spname=NULL,cmp.spname=NULL,disp.spname=NULL,drift.spname=NULL,sp.bin=NULL,binid.col=3)
  {
    ABP=data.frame(matrix(NA,nrow = ncol(comm),ncol = 3),stringsAsFactors = FALSE)
    rownames(ABP)=colnames(comm)
    colnames(ABP)=c("ABin","AProcess","AProcessB")
    ABP[match(sel.spname,colnames(comm)),2]="Selection"
    ABP[match(flt.spname,colnames(comm)),3]="Filtering"
    ABP[match(cmp.spname,colnames(comm)),3]="Competition"
    ABP[match(disp.spname,colnames(comm)),2:3]="Dispersal"
    ABP[match(drift.spname,colnames(comm)),2:3]="Drift"
    if(sum(is.na(ABP[,2:3]))>0){warning("Wrong! ABP has unexpected NA.")}
    if(!is.null(sp.bin))
    {
      ABP[match(sel.spname,colnames(comm)),1]=sp.bin[match(sel.spname,rownames(sp.bin)),binid.col]
      ABP[match(disp.spname,colnames(comm)),1]=sp.bin[match(disp.spname,rownames(sp.bin)),binid.col]
      ABP[match(drift.spname,colnames(comm)),1]=sp.bin[match(drift.spname,rownames(sp.bin)),binid.col]
      ABP[,1]=paste0("A",ABP[,1])
      if(sum(is.na(ABP))>0){warning("Wrong! ABP has unexpected NA.")}
    }
    ABP
  }
  
  
  try.fun<-function(fun,...)
  {
    tt=try(fun(...))
    t=1
    while(class(tt)=="try-error" & t<=6)
    {
      tt=try(fun(...))
      t=t+1
    }
    tt
  }
  
  #######################################
  # 2 # Scenario 1: selection
  #######################################
  
  s1.sim<-function(prefixi,...)
  {
    source(paste0(code.wd,"/taxa.bin.phy.r"))
    taxbin=taxa.bin.phy(tree,pd,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
    sp.bin.sel=taxbin$sp.bin[,3,drop=FALSE]
    
    spnamea=rownames(pd)
    commss=list()
    delta=rep(NA,10)
    for(t in 1:10)
    {
      commss[[t]]=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op, Nk.sel=Nk,sampname=sampname,code.wd=code.wd)
      comm=commss[[t]]$comm.sel
      delta[t]=length(spnamea)-ncol(comm)
    }
    commsel=commss[[which.min(delta)]]
    comm=commsel$comm.sel
    
    tx=1
    while(min(delta)>0 & tx<10)
    {
      message("testing tx=",tx," delta=",min(delta),". ",date())
      spnamea=colnames(comm)
      spc=match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      taxbin=taxa.bin.phy(tree = tree.use,pd = pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      sp.bin.sel=taxbin$sp.bin[,3,drop=FALSE]
      
      commss=list()
      delta=rep(NA,200)
      for(t in 1:200)
      {
        commss[[t]]=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op, Nk.sel=Nk,sampname=sampname,code.wd=code.wd)
        comm=commss[[t]]$comm.sel
        delta[t]=length(spnamea)-ncol(comm)
      }
      min(delta)
      commsel=commss[[which.min(delta)]]
      comm=commsel$comm.sel
      tx=tx+1
    }
    if(length(spnamea)!=ncol(comm)){stop("length(spnamea)!=ncol(comm) failed.")}
    dim(comm);dim(commsel$comm.flt);dim(commsel$comm.cmp)
    
    ABP=ABPcount(comm,sel.spname=colnames(comm),flt.spname=colnames(commsel$comm.flt),
                 cmp.spname=colnames(commsel$comm.cmp),disp.spname=NULL,drift.spname=NULL,
                 sp.bin=sp.bin.sel,binid.col=1)
    save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  }
  prefixi=paste0(prefix,".S1.Sel")
  try.fun(s1.sim,prefixi=prefixi)
  
  #######################################
  # 3 # Scenario 2: extreme dispersal
  #######################################
  # Sloan's dispersal model, pool.A and B no overlap, richness fixed
  ################################
  m1=0.01
  m2=0.99
  comm=dispersal.ZSM.sloan(op,Nk,sampname,spnumi=100,meta.ab=MZSM.meta,m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
  dim(comm)
  prefixi=paste0(prefix,".S2.Disp")
  ABP=ABPcount(comm = comm, disp.spname = colnames(comm))
  save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  
  #######################################
  # 4 # Scenario 3: drift
  #######################################
  # Sloan's dispersal model with MZSM metacommunity
  #######################################
  m1=0.5
  m2=0.01
  sig2.W=4
  comm=drift.ZSM.sloan(op,Nk,sampname,spnumi=100,meta.ab=MZSM.meta,
                       m1=m1,m2=m2,FH=0.95,FL=0.05,sig2.W=sig2.W)
  dim(comm)
  # test #
  prefixi=paste0(prefix,".S3.Drift")
  ABP=ABPcount(comm = comm, drift.spname = colnames(comm))
  save.sim(prefixi,comm,pd,tree,save.wd, ABP)
  
  #######################################
  # 5 # Scenario 4-6: some species under extreme dispersal, some under drift
  #######################################
  
  s46.sim<-function(a.disp,a.drift,prefixi,...)
  {
    id.disp=sample(1:length(op),length(op)*a.disp)
    meta.disp=MZSM.meta[id.disp]
    meta.drift=MZSM.meta[-id.disp]
    op.disp=rep(1,length(id.disp));names(op.disp)=paste0("disp",1:length(id.disp))
    op.drift=rep(1,length(op)-length(id.disp));names(op.drift)=paste0("drift",1:length(op.drift))
    
    Nk.disp=Nk*a.disp
    Nk.drift=Nk*a.drift
    spnumi.disp=spnumi*a.disp
    spnumi.drift=spnumi*a.drift
    
    m1=0.01
    m2=0.99
    nc.disp=1
    spn.disp=0
    t=1
    while(nc.disp!=spn.disp & t<10)
    {
      comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                    spnumi=spnumi.disp,meta.ab=meta.disp,
                                    m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
      (nc.disp=ncol(comm.disp))
      
      comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                 spnumi=spnumi.drift,meta.ab=meta.drift,
                                 m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
      (nc.drift=ncol(comm.drift))
      
      comm=cbind(comm.disp,comm.drift)
      spnuma=ncol(comm)
      spnamea=sample(tree$tip.label,spnuma)
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.disp+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      (tab.disp=tabspa[which.min(abs(tabspa-nc.disp))])
      spn.disp=tab.disp[[1]]
      spn.drift=spnuma-spn.disp
      
      k=1
      while(nc.disp!=spn.disp & k<100)
      {
        comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                      spnumi=spnumi.disp,meta.ab=meta.disp,
                                      m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
        nc.disp=ncol(comm.disp)
        message("Now nc.disp.delta=",nc.disp-spn.disp," t=",t," k=",k,". ",date())
        k=k+1
      }
      message("====Now nc.disp.delta=",nc.disp-spn.disp," t=",t,". ",date())
      t=t+1
    }
    if(nc.disp!=spn.disp){stop("nc.disp!=spn.disp failed")}
    
    tablet=as.numeric(substr(names(tab.disp)[1],2,regexpr("_",names(tab.disp)[1])-1))
    bin.sel=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.disp)[1],regexpr("_",names(tab.disp)[1])+1))]]
    
    spid.disp=which(taxbin$sp.bin$bin.id.new %in% bin.sel)
    colnames(comm.disp)=rownames(taxbin$sp.bin)[spid.disp]
    
    k=1
    while(nc.drift!=spn.drift & k<400)
    {
      comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                 spnumi=spnumi.drift,meta.ab=meta.drift,
                                 m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
      nc.drift=ncol(comm.drift)
      message("Now nc.drift.delta=",nc.drift-spn.drift," k=",k,". ",date())
      k=k+1
    }
    if(nc.drift!=spn.drift){stop("nc.drift!=spn.drift failed")}
    
    colnames(comm.drift)=rownames(taxbin$sp.bin)[-spid.disp]
    comm=cbind(comm.disp,comm.drift)
    rowSums(comm.disp);rowSums(comm.drift);rowSums(comm)
    
    ABP=ABPcount(comm,disp.spname = colnames(comm.disp),
                 drift.spname = colnames(comm.drift),
                 sp.bin = taxbin$sp.bin,binid.col = 3)
    save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  }
  
  a.disp=0.25 #0.5 #0.25
  a.drift=0.75 #0.5 #0.75
  prefixi=paste0(prefix,".S4.Disp25Drift75")
  try.fun(s46.sim,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  
  a.disp=0.5
  a.drift=0.5
  prefixi=paste0(prefix,".S5.Disp50Drift50")
  #s46.sim(a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  try.fun(s46.sim,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  
  a.disp=0.75
  a.drift=0.25
  prefixi=paste0(prefix,".S6.Disp75Drift25")
  try.fun(s46.sim,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  
  
  #####################################
  # 6 # Scenario 7-9: some species under selection, some under drift
  #####################################
  
  s79.sim<-function(a.sel,a.drift,prefixi,...)
  {
    id.drift=sample(1:length(op),round(length(op)*a.drift))
    meta.drift=MZSM.meta[id.drift]
    op.drift=rep(1,length(id.drift));names(op.drift)=paste0("drift",1:length(id.drift))
    Nk.drift=Nk*a.drift
    spnumi.drift=spnumi*a.drift
    
    nc.drift=1
    spn.drift=0
    t=1
    while(nc.drift!=spn.drift & t<10)
    {
      comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                 spnumi=spnumi.drift,meta.ab=meta.drift,
                                 m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
      (nc.drift=ncol(comm.drift))
      
      
      spnuma=nc.drift+(length(op)-length(op.drift))
      spnamea=sample(tree$tip.label,spnuma)
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      #save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.drift+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      
      (tab.drift=tabspa[which.min(abs(tabspa-nc.drift))])
      spn.drift=tab.drift[[1]]
      k=1
      while(nc.drift!=spn.drift & k<400)
      {
        message("Now nc.drift.delta=",nc.drift-spn.drift," k=",k,". ",date())
        comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                   spnumi=spnumi.drift,meta.ab=meta.drift,
                                   m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
        nc.drift=ncol(comm.drift)
        k=k+1
      }
      message("====Now nc.drift.delta=",nc.drift-spn.drift," t=",t,". ",date())
      t=t+1
    }
    if(nc.drift!=spn.drift){stop("nc.drift!=spn.drift failed.")}
    
    tablet=as.numeric(substr(names(tab.drift)[1],2,regexpr("_",names(tab.drift)[1])-1))
    bin.drift=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.drift)[1],regexpr("_",names(tab.drift)[1])+1))]]
    
    spid.drift=which(taxbin$sp.bin$bin.id.new %in% bin.drift)
    colnames(comm.drift)=rownames(taxbin$sp.bin)[spid.drift]
    
    sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
    sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% bin.drift)),,drop=FALSE]
    Nk.sel=Nk*a.sel
    comsel=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                   Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
    comm.sel=comsel$comm.sel
    
    comm=cbind(comm.sel,comm.drift)
    dim(comm)
    spnuma
    
    t=1
    while((ncol(comm)!=spnuma)& t<10)
    {
      spnamea=c(colnames(comm.drift),colnames(comm.sel))
      comm.drift.old=comm.drift
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.drift+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      (tab.drift=tabspa[which.min(abs(tabspa-nc.drift))])
      spn.drift=tab.drift[[1]]
      k=1
      while(nc.drift!=spn.drift & k<400)
      {
        comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                   spnumi=spnumi.drift,meta.ab=meta.drift,
                                   m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
        nc.drift=ncol(comm.drift)
        message("Now nc.drift.delta=",nc.drift-spn.drift," k=",k,". ",date())
        k=k+1
      }
      if(nc.drift!=spn.drift)
      {
        comm.drift=comm.drift.old
      }else{
        tablet=as.numeric(substr(names(tab.drift)[1],2,regexpr("_",names(tab.drift)[1])-1))
        bin.drift=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.drift)[1],regexpr("_",names(tab.drift)[1])+1))]]
        
        spid.drift=which(taxbin$sp.bin$bin.id.new %in% bin.drift)
        colnames(comm.drift)=rownames(taxbin$sp.bin)[spid.drift]
        
        sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
        sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% bin.drift)),,drop=FALSE]
        Nk.sel=Nk*a.sel
      }
      comss=list()
      deltax=rep(NA,200)
      for(txx in 1:200)
      {
        comss[[txx]]=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                             Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
        deltax[txx]=length(spnamea)-(ncol(comss[[txx]]$comm.sel)+ncol(comm.drift))
      }
      min(deltax)
      comsel=comss[[which.min(abs(deltax))]]
      comm.sel=comsel$comm.sel
      comm=cbind(comm.sel,comm.drift)
      dim(comm)
      spnuma=length(spnamea)
      print(dim(comm))
      print(spnuma)
      message("-----t=",t,". dert=",ncol(comm)-spnuma,". ",date())
      t=t+1
    }
    if(ncol(comm)!=spnuma){stop("(ncol(comm)!=spnuma) failed")}
    ABP=ABPcount(comm,sel.spname = colnames(comm.sel),
                 flt.spname = colnames(comsel$comm.flt),
                 cmp.spname = colnames(comsel$comm.cmp),
                 drift.spname = colnames(comm.drift),
                 sp.bin = taxbin$sp.bin,binid.col = 3)
    save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  }
  
  a.sel=0.25
  a.drift=0.75
  prefixi=paste0(prefix,".S7.Sel25Drift75")
  try.fun(s79.sim,a.sel=a.sel,a.drift=a.drift,prefixi=prefixi)
  
  a.sel=0.5
  a.drift=0.5
  prefixi=paste0(prefix,".S8.Sel50Drift50")
  try.fun(s79.sim,a.sel=a.sel,a.drift=a.drift,prefixi=prefixi)
  
  a.sel=0.75
  a.drift=0.25
  prefixi=paste0(prefix,".S9.Sel75Drift25")
  try.fun(s79.sim,a.sel=a.sel,a.drift=a.drift,prefixi=prefixi)
  
  #####################################
  # 7 # Scenario 10-12: some species under selection, some under extreme dispersal
  #####################################
  
  s1012.sim<-function(a.sel,a.disp,prefixi,...)
  {
    sig2.L=0.015
    sig2.H=0.015
    id.disp=sample(1:length(op),round(length(op)*a.disp))
    meta.disp=MZSM.meta[id.disp]
    op.disp=rep(1,length(id.disp));names(op.disp)=paste0("disp",1:length(id.disp))
    Nk.disp=Nk*a.disp
    spnumi.disp=spnumi*a.disp
    
    m1=0.01
    m2=0.99
    
    nc.disp=1
    spn.disp=0
    t=1
    while(nc.disp!=spn.disp & t<10)
    {
      comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                    spnumi=spnumi.disp,meta.ab=meta.disp,
                                    m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
      (nc.disp=ncol(comm.disp))
      
      
      spnuma=nc.disp+(length(op)-length(op.disp))
      spnamea=sample(tree$tip.label,spnuma)
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.disp+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      
      (tab.disp=tabspa[which.min(abs(tabspa-nc.disp))])
      spn.disp=tab.disp[[1]]
      k=1
      while(nc.disp!=spn.disp & k<400)
      {
        comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                      spnumi=spnumi.disp,meta.ab=meta.disp,
                                      m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
        nc.disp=ncol(comm.disp)
        message("Now nc.disp.delta=",nc.disp-spn.disp," k=",k,". ",date())
        k=k+1
      }
      message("====Now nc.disp.delta=",nc.disp-spn.disp," t=",t,". ",date())
      t=t+1
    }
    if(nc.disp!=spn.disp){stop("nc.disp!=spn.disp failed")}
    
    tablet=as.numeric(substr(names(tab.disp)[1],2,regexpr("_",names(tab.disp)[1])-1))
    bin.disp=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.disp)[1],regexpr("_",names(tab.disp)[1])+1))]]
    
    spid.disp=which(taxbin$sp.bin$bin.id.new %in% bin.disp)
    colnames(comm.disp)=rownames(taxbin$sp.bin)[spid.disp]
    
    sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
    sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% bin.disp)),,drop=FALSE]
    Nk.sel=Nk*a.sel
    comsel=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                   Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
    comm.sel=comsel$comm.sel
    
    comm=cbind(comm.sel,comm.disp)
    dim(comm)
    spnuma
    
    t=1
    while((ncol(comm)!=spnuma) & t<10)
    {
      spnamea=c(colnames(comm.disp),colnames(comm.sel))
      comm.disp.old=comm.disp
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.disp+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      
      (tab.disp=tabspa[which.min(abs(tabspa-nc.disp))])
      
      spn.disp=tab.disp[[1]]
      k=1
      while(nc.disp!=spn.disp & k<400)
      {
        comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                      spnumi=spnumi.disp,meta.ab=meta.disp,
                                      m1=m1,m2=m2,distinct=TRUE,fix.rich = TRUE)
        nc.disp=ncol(comm.disp)
        message("Now nc.disp.delta=",nc.disp-spn.disp," k=",k,". ",date())
        k=k+1
      }
      if(nc.disp!=spn.disp)
      {
        comm.disp=comm.disp.old
      }else{
        tablet=as.numeric(substr(names(tab.disp)[1],2,regexpr("_",names(tab.disp)[1])-1))
        bin.disp=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.disp)[1],regexpr("_",names(tab.disp)[1])+1))]]
        
        spid.disp=which(taxbin$sp.bin$bin.id.new %in% bin.disp)
        colnames(comm.disp)=rownames(taxbin$sp.bin)[spid.disp]
        
        sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
        sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% bin.disp)),,drop=FALSE]
        Nk.sel=Nk*a.sel
      }
      comss=list()
      deltax=rep(NA,200)
      for(txx in 1:200)
      {
        comss[[txx]]=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                             Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
        deltax[txx]=length(spnamea)-(ncol(comss[[txx]]$comm.sel)+ncol(comm.disp))
      }
      min(deltax)
      comsel=comss[[which.min(abs(deltax))]]
      comm.sel=comsel$comm.sel
      
      comm=cbind(comm.sel,comm.disp)
      spnuma=length(spnamea)
      print(dim(comm))
      print(spnuma)
      message("-----t=",t,". dert=",ncol(comm)-spnuma,". ",date())
      t=t+1
    }
    if(ncol(comm)!=spnuma){stop("(ncol(comm)!=spnuma) failed.")}
    ABP=ABPcount(comm,sel.spname = colnames(comm.sel),
                 flt.spname = colnames(comsel$comm.flt),
                 cmp.spname = colnames(comsel$comm.cmp),
                 disp.spname = colnames(comm.disp),
                 sp.bin = taxbin$sp.bin,binid.col = 3)
    save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  }
  
  a.sel=0.25
  a.disp=0.75
  prefixi=paste0(prefix,".S10.Sel25Disp75")
  try.fun(s1012.sim,a.sel=a.sel,a.disp=a.disp,prefixi=prefixi)
  
  a.sel=0.5
  a.disp=0.5
  prefixi=paste0(prefix,".S11.Sel50Disp50")
  try.fun(s1012.sim,a.sel=a.sel,a.disp=a.disp,prefixi=prefixi)
  
  a.sel=0.75
  a.disp=0.25
  prefixi=paste0(prefix,".S12.Sel75Disp25")
  try.fun(s1012.sim,a.sel=a.sel,a.disp=a.disp,prefixi=prefixi)
  
  #####################################
  # 8 # Scenario 13-15: some species under selection, some under extreme dispersal, and others under drift
  #####################################
  
  s1315.sim<-function(a.sel,a.disp,a.drift,prefixi,...)
  {
    id.disp=sample(1:length(op),round(length(op)*a.disp))
    meta.disp=MZSM.meta[id.disp]
    op.disp=rep(1,length(id.disp));names(op.disp)=paste0("disp",1:length(id.disp))
    Nk.disp=Nk*a.disp
    spnumi.disp=spnumi*a.disp
    
    id.drift=sample((1:length(op))[-id.disp],round(length(op)*a.drift))
    meta.drift=MZSM.meta[id.drift]
    op.drift=rep(1,length(id.drift));names(op.drift)=paste0("drift",1:length(id.drift))
    Nk.drift=Nk*a.drift
    spnumi.drift=spnumi*a.drift
    
    
    m1.disp=0.01
    m2.disp=0.99
    comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                  spnumi=spnumi.disp,meta.ab=meta.disp,
                                  m1=m1.disp,m2=m2.disp,distinct=TRUE,fix.rich = TRUE)
    (nc.disp=ncol(comm.disp))
    
    m1.drift=0.5
    m2.drift=0.01
    comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                               spnumi=spnumi.drift,meta.ab=meta.drift,
                               m1=m1.drift,m2=m2.drift,FH=0.95,FL=0.05,sig2.W=4)
    (nc.drift=ncol(comm.drift))
    
    spnuma=nc.disp+nc.drift+(length(op)-length(op.disp)-length(op.drift))
    spnamea=sample(tree$tip.label,spnuma)
    spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
    tree.use=spc$tree
    pd.use=spc$pd
    source(paste0(code.wd,"/taxa.bin.phy.r"))
    taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
    save(taxbin,file=paste0(prefixi,".taxabin.rda"))
    tab.sp=table(taxbin$sp.bin$bin.id.new)
    
    bin.num.m=min(ceiling((a.disp+0.1)*length(tab.sp)),length(tab.sp))
    combx<-tab.spx<-list()
    for(x in 1:bin.num.m)
    {
      nnx=choose(length(tab.sp),x)
      if(nnx<1000)
      {
        combx[[x]]=combn(length(tab.sp),x)
        tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
      }else{
        combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
        tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
      }
      names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
    }
    tabspa=unlist(tab.spx)
    
    (tab.disp=tabspa[which.min(abs(tabspa-nc.disp))])
    
    spn.disp=tab.disp[[1]]
    k=1
    while(nc.disp!=spn.disp & k<400)
    {
      comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                    spnumi=spnumi.disp,meta.ab=meta.disp,
                                    m1=m1.disp,m2=m2.disp,distinct=TRUE,fix.rich = TRUE)
      nc.disp=ncol(comm.disp)
      message("Now nc.disp.delta=",nc.disp-spn.disp," k=",k,". ",date())
      k=k+1
    }
    
    tablet=as.numeric(substr(names(tab.disp)[1],2,regexpr("_",names(tab.disp)[1])-1))
    bin.disp=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.disp)[1],regexpr("_",names(tab.disp)[1])+1))]]
    
    spid.disp=which(taxbin$sp.bin$bin.id.new %in% bin.disp)
    colnames(comm.disp)=rownames(taxbin$sp.bin)[spid.disp]
    
    tab.sp.old=tab.sp
    tab.sp=tab.sp.old[!(names(tab.sp.old) %in% bin.disp)]
    
    bin.num.m=min(ceiling((a.drift+0.1)*length(tab.sp)),length(tab.sp))
    combx<-tab.spx<-list()
    for(x in 1:bin.num.m)
    {
      nnx=choose(length(tab.sp),x)
      if(nnx<1000)
      {
        combx[[x]]=combn(length(tab.sp),x)
        tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
      }else{
        combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
        tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
      }
      names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
    }
    tabspa=unlist(tab.spx)
    
    (tab.drift=tabspa[which.min(abs(tabspa-nc.drift))])
    
    spn.drift=tab.drift[[1]]
    k=1
    while(nc.drift!=spn.drift & k<400)
    {
      comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                 spnumi=spnumi.drift,meta.ab=meta.drift,
                                 m1=m1.drift,m2=m2.drift,FH=0.95,FL=0.05,sig2.W=4)
      nc.drift=ncol(comm.drift)
      message("Now nc.drift.delta=",nc.drift-spn.drift," k=",k,". ",date())
      k=k+1
    }
    if(nc.drift!=spn.drift){stop("nc.drift!=spn.drift failed.")}
    
    tablet=as.numeric(substr(names(tab.drift)[1],2,regexpr("_",names(tab.drift)[1])-1))
    bin.drift=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.drift)[1],regexpr("_",names(tab.drift)[1])+1))]]
    
    spid.drift=which(taxbin$sp.bin$bin.id.new %in% bin.drift)
    colnames(comm.drift)=rownames(taxbin$sp.bin)[spid.drift]
    
    sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
    sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% c(bin.drift,bin.disp))),,drop=FALSE]
    Nk.sel=Nk*a.sel
    comsel=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                   Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
    comm.sel=comsel$comm.sel
    
    comm=cbind(comm.sel,comm.disp,comm.drift)
    dim(comm)
    spnuma
    
    t=1
    while((ncol(comm)!=spnuma)& t<10)
    {
      spnamea=c(colnames(comm.disp),colnames(comm.drift),colnames(comm.sel))
      comm.disp.old=comm.disp
      comm.drift.old=comm.drift
      spc=mCAM::match.name(name.check = spnamea,both.list = list(pd=pd),tree.list = list(tree=tree),silent = TRUE)
      tree.use=spc$tree
      pd.use=spc$pd
      source(paste0(code.wd,"/taxa.bin.phy.r"))
      taxbin=taxa.bin.phy(tree.use,pd.use,d.cut=d.cut,bin.size.limit = bin.size,nworker = nworker,code.wd = code.wd)
      save(taxbin,file=paste0(prefixi,".taxabin.rda"))
      tab.sp=table(taxbin$sp.bin$bin.id.new)
      
      bin.num.m=min(ceiling((a.disp+0.1)*length(tab.sp)),length(tab.sp))
      combx<-tab.spx<-list()
      for(x in 1:bin.num.m)
      {
        nnx=choose(length(tab.sp),x)
        if(nnx<1000)
        {
          combx[[x]]=combn(length(tab.sp),x)
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }else{
          combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
          tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
        }
        names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
      }
      tabspa=unlist(tab.spx)
      (tab.disp=tabspa[which.min(abs(tabspa-nc.disp))])
      
      spn.disp=tab.disp[[1]]
      k=1
      while(nc.disp!=spn.disp & k<400)
      {
        comm.disp=dispersal.ZSM.sloan(op=op.disp,Nk=Nk.disp,sampname=sampname,
                                      spnumi=spnumi.disp,meta.ab=meta.disp,
                                      m1=m1.disp,m2=m2.disp,distinct=TRUE,fix.rich = TRUE)
        nc.disp=ncol(comm.disp)
        message("Now nc.disp.delta=",nc.disp-spn.disp," k=",k,". ",date())
        k=k+1
      }
      if(nc.disp!=spn.disp)
      {
        comm.disp=comm.disp.old
      }else{
        tablet=as.numeric(substr(names(tab.disp)[1],2,regexpr("_",names(tab.disp)[1])-1))
        bin.disp=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.disp)[1],regexpr("_",names(tab.disp)[1])+1))]]
        
        spid.disp=which(taxbin$sp.bin$bin.id.new %in% bin.disp)
        colnames(comm.disp)=rownames(taxbin$sp.bin)[spid.disp]
        
        tab.sp.old=tab.sp
        tab.sp=tab.sp.old[!(names(tab.sp.old) %in% bin.disp)]
        
        bin.num.m=min(ceiling((a.drift+0.1)*length(tab.sp)),length(tab.sp))
        combx<-tab.spx<-list()
        for(x in 1:bin.num.m)
        {
          nnx=choose(length(tab.sp),x)
          if(nnx<1000)
          {
            combx[[x]]=combn(length(tab.sp),x)
            tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
          }else{
            combx[[x]]=sapply(1:1000,function(xx){sample(length(tab.sp),x)})
            tab.spx[[x]]=sapply(1:ncol(combx[[x]]),function(i){sum(tab.sp[combx[[x]][,i]])})
          }
          names(tab.spx[[x]])=paste0("x",x,"_",1:length(tab.spx[[x]]))
        }
        tabspa=unlist(tab.spx)
        
        (tab.drift=tabspa[which.min(abs(tabspa-nc.drift))])
        
        spn.drift=tab.drift[[1]]
        k=1
        while(nc.drift!=spn.drift & k<400)
        {
          message("Now nc.drift.delta=",nc.drift-spn.drift," k=",k,". ",date())
          comm.drift=drift.ZSM.sloan(op=op.drift,Nk=Nk.drift,sampname=sampname,
                                     spnumi=spnumi.drift,meta.ab=meta.drift,
                                     m1=m1.drift,m2=m2.drift,FH=0.95,FL=0.05,sig2.W=4)
          nc.drift=ncol(comm.drift)
          k=k+1
        }
        if(nc.drift!=spn.drift)
        {
          comm.drift=comm.drift.old
          comm.disp=comm.disp.old
        }else{
          tablet=as.numeric(substr(names(tab.drift)[1],2,regexpr("_",names(tab.drift)[1])-1))
          bin.drift=names(tab.sp)[combx[[tablet]][,as.numeric(substring(names(tab.drift)[1],regexpr("_",names(tab.drift)[1])+1))]]
          
          spid.drift=which(taxbin$sp.bin$bin.id.new %in% bin.drift)
          colnames(comm.drift)=rownames(taxbin$sp.bin)[spid.drift]
          
          sp.bin.all=taxbin$sp.bin[,3,drop=FALSE]
          sp.bin.sel=sp.bin.all[which(!(sp.bin.all[,1] %in% c(bin.drift,bin.disp))),,drop=FALSE]
          Nk.sel=Nk*a.sel
        }
      }
      comss=list()
      deltax=rep(NA,200)
      for(txx in 1:200)
      {
        comss[[txx]]=com.sel(sp.bin.sel=sp.bin.sel, compet.sel=compet.sel, op=op,
                             Nk.sel=Nk.sel,sampname=sampname,code.wd=code.wd)
        deltax[txx]=length(spnamea)-(ncol(comss[[txx]]$comm.sel)+ncol(comm.disp)+ncol(comm.drift))
      }
      min(deltax)
      comsel=comss[[which.min(abs(deltax))]]
      comm.sel=comsel$comm.sel
      
      comm=cbind(comm.sel,comm.disp,comm.drift)
      spnuma=length(spnamea)
      print(dim(comm))
      print(spnuma)
      message("-----t=",t,". dert=",ncol(comm)-spnuma,". ",date())
      t=t+1
    }
    if(ncol(comm)!=spnuma){stop("ncol(comm)!=spnuma failed.")}
    ABP=ABPcount(comm,sel.spname = colnames(comm.sel),
                 flt.spname = colnames(comsel$comm.flt),
                 cmp.spname = colnames(comsel$comm.cmp),
                 disp.spname = colnames(comm.disp),
                 drift.spname = colnames(comm.drift),
                 sp.bin = taxbin$sp.bin,binid.col = 3)
    save.sim(prefixi,comm,pd,tree,save.wd,ABP)
  }
  
  a.sel=0.25
  a.disp=0.25
  a.drift=0.5
  prefixi=paste0(prefix,".S13.Sel25Disp25Drift50")
  try.fun(s1315.sim,a.sel=a.sel,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  
  a.sel=0.25
  a.disp=0.5
  a.drift=0.25
  prefixi=paste0(prefix,".S14.Sel25Disp50Drift25")
  try.fun(s1315.sim,a.sel=a.sel,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  
  a.sel=0.5
  a.disp=0.25
  a.drift=0.25
  prefixi=paste0(prefix,".S15.Sel50Disp25Drift25")
  try.fun(s1315.sim,a.sel=a.sel,a.disp=a.disp,a.drift=a.drift,prefixi=prefixi)
  save.wd
}

###########################################
# To set which type of trait to use
# b6 meas actual bin size limit is 6 
# d02 means ds=0.2
bin.size=12
ds=0.2

# Low phylogenetic signal scenario # JS
op=ops$JS$openv
compet.sels=c(0.25,0.5,0.75,1)
for(cs in 2:length(compet.sels))
{
  compet.sel=compet.sels[cs]
  prefix=paste0("JSb",bin.size,"d02C",compet.sel*100)
  simmu(prefix=prefix,op=op,bin.size=bin.size,compet.sel=compet.sel)
}

# Medium phylogenetic signal scenario # BM
op=ops$BM$openv
compet.sels=c(0.25,0.5,0.75,1)
for(cs in 1:length(compet.sels))
{
  compet.sel=compet.sels[cs]
  prefix=paste0("BMb",bin.size,"d02C",compet.sel*100)
  simmu(prefix=prefix,op=op,bin.size=bin.size,compet.sel=compet.sel)
}

# High phylogenetic signal scenario # AD
op=ops$ACDC$openv
compet.sels=c(0.25,0.5,0.75,1)
for(cs in 1:length(compet.sels))
{
  compet.sel=compet.sels[cs]
  prefix=paste0("ADb",bin.size,"d02C",compet.sel*100)
  simmu(prefix=prefix,op=op,bin.size=bin.size,compet.sel=compet.sel)
}

# End #