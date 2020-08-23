evaluate.sim<-function(icamp.res.detail,ABP,comm,sel=FALSE)
{
  icbin=iCAMP::icamp.bins(icamp.detail = icamp.res.detail)
  Wtuvk=icbin$Wtuvk
  head(Wtuvk)
  pairsamps=Wtuvk[,c("name1","name2"),drop=FALSE]
  head(pairsamps)
  pron=c("HoS","HeS","HD","DL","DR")
  if(sel){pron2=c("Sel","HD","DL","DR")}
  # 1 # species level
  # 1.1 # expectation and abundance
  Puvi.l1=lapply(1:nrow(pairsamps),
                function(uv)
                {
                  if(sel)
                  {
                    cbind(name1=pairsamps[uv,1],name2=pairsamps[uv,2],species=rownames(ABP),ptype=ABP[,2:3])
                  }else{
                    cbind(name1=pairsamps[uv,1],name2=pairsamps[uv,2],species=rownames(ABP),ptype=ABP[,2])
                  }
                })
  Puvi.m1=Reduce(rbind,Puvi.l1)
  env1=substr(Puvi.m1[,1],1,1)
  env2=substr(Puvi.m1[,2],1,1)
  site1=substr(Puvi.m1[,1],2,2)
  site2=substr(Puvi.m1[,2],2,2)
  Pexp.uvi=rep(NA,nrow(Puvi.m1))
  Pexp.uvi[which((Puvi.m1[,4]=="Selection")&(env1==env2))]="HoS"
  Pexp.uvi[which((Puvi.m1[,4]=="Selection")&(env1!=env2))]="HeS"
  Pexp.uvi[which((Puvi.m1[,4]=="Dispersal")&(site1==site2))]="HD"
  Pexp.uvi[which((Puvi.m1[,4]=="Dispersal")&(site1!=site2))]="DL"
  Pexp.uvi[which(Puvi.m1[,4]=="Drift")]="DR"
  if(sel)
  {
    Pexp2.uvi=Pexp.uvi
    Pexp2.uvi[which(Puvi.m1[,4]=="Selection")]="Sel"
  }
  Ab.uvi=comm[cbind(match(Puvi.m1[,1],rownames(comm)),match(Puvi.m1[,3],colnames(comm)))]+comm[cbind(match(Puvi.m1[,2],rownames(comm)),match(Puvi.m1[,3],colnames(comm)))]
  # 1.2 # identified process for each turnover of each species
  sp.bin=icamp.res.detail$taxabin$sp.bin[,3,drop=FALSE]
  bin.uvi=sp.bin[match(Puvi.m1[,3],rownames(sp.bin)),1]
  samp2.uvi=paste0(Puvi.m1[,1],"__",Puvi.m1[,2])
  samp2.Wtuvk=paste0(Wtuvk[,"name1"],"__",Wtuvk[,"name2"])
  Pobs.uvi=Wtuvk[cbind(match(samp2.uvi,samp2.Wtuvk),match(paste0("bin",bin.uvi),colnames(Wtuvk)))]
  if(sel){Pobs2.uvi=Pobs.uvi;Pobs2.uvi[which(Pobs.uvi %in% c("HoS","HeS"))]="Sel"}
  Puvi=data.frame(Puvi.m1,Abundance=Ab.uvi,ExpectedPuvi=Pexp.uvi,ObsPuvi=Pobs.uvi,stringsAsFactors = FALSE)
  if(sel){Puvi=data.frame(Puvi,ExpectedPuvi2=Pexp2.uvi,ObsPuvi2=Pobs2.uvi,stringsAsFactors = FALSE)}
  quali.perf<-function(exp.v,PNmexp=NULL,obs.v,PNm=NULL,ab.v,pron=c("HoS","HeS","HD","DL","DR"))
  {
    if(length(unique(c(length(exp.v),nrow(PNmexp),length(obs.v),nrow(PNm),length(ab.v),0)))!=2){stop("Different length!")}
    TFm<-matrix(NA,nrow=length(ab.v),ncol=length(pron))
    if(is.null(PNmexp))
    {
      PNmexp<-TFm
      PNmexp[]<-"N"
      PNmexp[cbind(1:nrow(PNmexp),match(exp.v,pron))]="P"
    }else{
      PNmexp=PNmexp[,match(pron,colnames(PNmexp))]
    }
    if(is.null(PNm))
    {
      PNm<-TFm
      PNm[]<-"N"
      PNm[cbind(1:nrow(PNm),match(obs.v,pron))]="P"
    }else{
      PNm=PNm[,match(pron,colnames(PNm))]
    }
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
    #ACC=(TPab+TNab)/(TPab+TNab+FPab+FNab)
    #PRC=TPab/(TPab+FPab)
    #SPC=TNab/(TNab+FPab)
    #SST=TPab/(TPab+FNab)
    #out=rbind(TPab,TNab,FPab,FNab,ACC,PRC,SPC,SST)
    out=rbind(TPab,TNab,FPab,FNab)
    colnames(out)=c(pron,"Total")
    out
  }
  
  Puvi.perf=quali.perf(exp.v = Puvi$ExpectedPuvi,obs.v = Puvi$ObsPuvi, ab.v = Puvi$Abundance)
  if(sel)
  {
    Puvi.perf2=quali.perf(exp.v = Puvi$ExpectedPuvi2 ,obs.v = Puvi$ObsPuvi2, ab.v = Puvi$Abundance, pron = pron2)
  }
  # 2 # Each turnover
  # 2.1 # bin level
  bins=colnames(Wtuvk)[4:ncol(Wtuvk)]
  
  Puvk.l1=lapply(1:nrow(pairsamps),
                 function(uv)
                 {
                   cbind(name1=pairsamps[uv,1],name2=pairsamps[uv,2],bin=bins)
                 })
  Puvk.m1=Reduce(rbind,Puvk.l1)
  bin.ra=icamp.res.detail$bin.weight
  samp2.uvk=paste0(Puvk.m1[,1],"__",Puvk.m1[,2])
  comb2m<-function(v1,v2,comb1,sep="")
  {
    comb2=paste(v1,v2,sep = sep)
    id.xx=which(!(comb2 %in% comb1))
    if(length(id.xx)>0)
    {
      comb2b=paste(v2,v1,sep = sep)
      comb2[id.xx]=comb2b[id.xx]
    }
    comb2
  }
  samp2.binra=comb2m(bin.ra[,1],bin.ra[,2],samp2.uvk,sep = "__")
  RA.uvk=as.numeric(bin.ra[cbind(match(samp2.uvk,samp2.binra),match(Puvk.m1[,3],colnames(bin.ra)))])
  
  qp.ab<-function(process.v,ab.v,pron=c("HoS","HeS","HD","DL","DR"))
  {
    if(sum(ab.v)==0)
    {
      outx=rep(0,length(pron))
    }else{
      outx=sapply(1:length(pron),
                  function(i)
                  {
                    sum(ab.v[which(process.v==pron[i])])/sum(ab.v)
                  })
    }
    names(outx)=pron
    outx
  }
  
  trac=seq(from=1,to=nrow(Puvk.m1),by=500)
  Pexp.uvk=t(sapply(1:nrow(Puvk.m1),
                    function(i)
                    {
                      if(i %in% trac){message("i=",i," in ",nrow(Puvk.m1),". ",date())}
                      spk=rownames(sp.bin)[sp.bin[,1]==sub("bin","",Puvk.m1[i,3])]
                      id.uvk=which(samp2.uvi==paste0(Puvk.m1[i,1],"__",Puvk.m1[i,2]) & (Puvi$species %in% spk))
                      qp.ab(process.v = Puvi$ExpectedPuvi[id.uvk],ab.v = Puvi$Abundance[id.uvk])
                    }))
  MixExp.uvk=as.integer(rowSums(Pexp.uvk!=0)>1)
  Pexp.maxk=apply(Pexp.uvk,1,max)
  WPNexp.uvk=Pexp.uvk
  WPNexp.uvk[]="N"
  WPNexp.uvk[(Pexp.uvk==Pexp.maxk)]="P"
  Wobs.uvk= Wtuvk[cbind(match(samp2.uvk,samp2.Wtuvk),match(Puvk.m1[,3],colnames(Wtuvk)))]
  Wtuvk.perf=quali.perf(exp.v=NULL,PNmexp=WPNexp.uvk,obs.v=Wobs.uvk,ab.v=RA.uvk)
  colnames(Pexp.uvk)=paste0(colnames(Pexp.uvk),".PtExpect")
  colnames(WPNexp.uvk)=paste0(colnames(WPNexp.uvk),".WtExpect")
  
  if(sel)
  {
    Pexp2.uvk=cbind(Sel=rowSums(Pexp.uvk[,1:2,drop=FALSE]),Pexp.uvk[,3:5,drop=FALSE])
    Pexp2.maxk=apply(Pexp2.uvk,1,max)
    WPNexp2.uvk=Pexp2.uvk
    WPNexp2.uvk[]="N"
    WPNexp2.uvk[(Pexp2.uvk==Pexp2.maxk)]="P"
    Wobs2.uvk=Wobs.uvk
    Wobs2.uvk[which(Wobs.uvk %in% c("HoS","HeS"))]="Sel"
    Wtuvk.perf2=quali.perf(exp.v=NULL,PNmexp=WPNexp2.uvk,obs.v=Wobs2.uvk,ab.v=RA.uvk,pron = pron2)
    colnames(Pexp2.uvk)=paste0(colnames(Pexp2.uvk),".PtExpect2")
    colnames(WPNexp2.uvk)=paste0(colnames(WPNexp2.uvk),".WtExpect2")
  }
  Puvk=data.frame(Puvk.m1,RelativeAbundance=RA.uvk,Pexp.uvk,MixedProcesses=MixExp.uvk,WPNexp.uvk,ObsPuvk=Wobs.uvk,stringsAsFactors = FALSE)
  if(sel)
  {
    Puvk=data.frame(Puvk,Pexp2.uvk,WPNexp2.uvk,ObsPuvk2=Wobs2.uvk,stringsAsFactors = FALSE)
  }
  # 2.2 # community level
  Puv.m1=pairsamps
  Pexp.uv=t(sapply(1:nrow(pairsamps),
                   function(uv)
                   {
                     id.uv=which(samp2.uvk==samp2.Wtuvk[uv])
                     (colSums(Pexp.uvk[id.uv,]*RA.uvk[id.uv]))/sum(RA.uvk[id.uv])
                   }))
  Pexp.max=apply(Pexp.uv,1,max)
  WPNexp.uv=Pexp.uv
  WPNexp.uv[]="N"
  WPNexp.uv[(Pexp.uv==Pexp.max)]="P"
  colnames(WPNexp.uv)=gsub(".PtExpect$","",colnames(WPNexp.uv))
  
  Ptuv=icbin$Ptuv
  samp2.ptuv=comb2m(Ptuv[,2],Ptuv[,3],samp2.Wtuvk,sep="__")
  Pobs.uv=as.matrix(Ptuv[match(samp2.Wtuvk,samp2.ptuv),match(pron,colnames(Ptuv))])
  Pobs.max=apply(Pobs.uv,1,max)
  WPNobs.uv=Pobs.uv
  WPNobs.uv[]="N"
  WPNobs.uv[(Pobs.uv==Pobs.max)]="P" 
  head(WPNobs.uv)  
  Wtuv.perf=quali.perf(exp.v=NULL,PNmexp=WPNexp.uv,obs.v=NULL,PNm=WPNobs.uv,ab.v=rep(1,nrow(WPNobs.uv)))
  
  if(sel)
  {
    Pexp2.uv=cbind(Sel.PtExpect=rowSums(Pexp.uv[,1:2,drop=FALSE]),Pexp.uv[,3:5,drop=FALSE])
    colnames(Pexp2.uv)=paste0(colnames(Pexp2.uv),"2")
    Pexp2.max=apply(Pexp2.uv,1,max)
    WPNexp2.uv=Pexp2.uv
    WPNexp2.uv[]="N"
    WPNexp2.uv[(Pexp2.uv==Pexp2.max)]="P"
    colnames(WPNexp2.uv)=gsub(".PtExpect2$","",colnames(WPNexp2.uv))
    
    Pobs2.uv=cbind(Sel=rowSums(Pobs.uv[,1:2,drop=FALSE]),Pobs.uv[,3:5,drop=FALSE])
    Pobs2.max=apply(Pobs2.uv,1,max)
    WPNobs2.uv=Pobs2.uv
    WPNobs2.uv[]="N"
    WPNobs2.uv[(Pobs2.uv==Pobs2.max)]="P" 
    head(WPNobs2.uv)  
    Wtuv.perf2=quali.perf(exp.v=NULL,PNmexp=WPNexp2.uv,obs.v=NULL,PNm=WPNobs2.uv,ab.v=rep(1,nrow(WPNobs2.uv)),pron = pron2)
  }
  
  colnames(WPNexp.uv)=paste0(colnames(WPNexp.uv),".WtExpect")
  colnames(WPNobs.uv)=paste0(colnames(WPNobs.uv),".WtObs")
  colnames(Pobs.uv)=paste0(colnames(Pobs.uv),".PtObs")
  Puv=data.frame(pairsamps,Pexp.uv,WPNexp.uv,Pobs.uv,WPNobs.uv,stringsAsFactors = FALSE)
  
  if(sel)
  {
    colnames(WPNexp2.uv)=paste0(colnames(WPNexp2.uv),".WtExpect2")
    colnames(WPNobs2.uv)=paste0(colnames(WPNobs2.uv),".WtObs2")
    colnames(Pobs2.uv)=paste0(colnames(Pobs2.uv),".PtObs2")
    Puv=data.frame(Puv,Pexp2.uv,WPNexp2.uv,Pobs2.uv,WPNobs2.uv,stringsAsFactors = FALSE)
  }
  
  
  # 3 # Each plot and each comparison between plots
  # 3.1 # bin level
  pairplots=c("HAHA","HBHB","LALA","LBLB","HAHB","HALA","HALB","HBLA","HBLB","LALB")
  plot2.uvk=comb2m(substr(Puvk[,1],1,2),substr(Puvk[,2],1,2),pairplots,sep="")
  Pgk.m1=data.frame(plots=rep(pairplots,each=length(bins)),bin=rep(bins,times=length(pairplots)))
  
  Pexpobs.gk=t(sapply(1:nrow(Pgk.m1),
                   function(i)
                   {
                     idik=which(plot2.uvk==Pgk.m1[i,1] & Puvk[,"bin"]==Pgk.m1[i,2])
                     out.exp=colSums(Puvk[idik,match(paste0(pron,".PtExpect"),colnames(Puvk))]*Puvk$RelativeAbundance[idik])/sum(Puvk$RelativeAbundance[idik])
                     out.obs=qp.ab(process.v = Puvk$ObsPuvk[idik],ab.v = Puvk$RelativeAbundance[idik])
                     c(RA.bin=sum(Puvk$RelativeAbundance[idik]),out.exp,out.obs)
                   }))
  colnames(Pexpobs.gk)[7:11]=paste0(colnames(Pexpobs.gk)[7:11],".PtObs")
  Pgk=data.frame(Pgk.m1,Pexpobs.gk,stringsAsFactors = FALSE)
  if(sel)
  {
    Pgk=data.frame(Pgk,Sel.PtExpect2=Pgk$HoS.PtExpect+Pgk$HeS.PtExpect,Sel.PtObs2=Pgk$HoS.PtObs+Pgk$HeS.PtObs,stringsAsFactors = FALSE)
  }
  
  # 3.2 # community level
  Pexpobs.g=t(sapply(1:length(pairplots),
                     function(i)
                     {
                       idi=which(Pgk$plots==pairplots[i])
                       colSums(Pgk[idi,4:ncol(Pgk)]*Pgk$RA.bin[idi])/sum(Pgk$RA.bin[idi])
                     }))
  Pg=data.frame(plots=pairplots,Pexpobs.g,stringsAsFactors = FALSE)
  
  output=list(Pg=Pg,Pgk=Pgk,Puv=Puv,Puvk=Puvk,Puvi=Puvi,
              TNPuv=Wtuv.perf,TNPuvk=Wtuvk.perf,TNPuvi=Puvi.perf)
  if(sel)
  {
    output=list(Pg=Pg,Pgk=Pgk,Puv=Puv,Puvk=Puvk,Puvi=Puvi,
                TNPuv=Wtuv.perf,TNPuvk=Wtuvk.perf,TNPuvi=Puvi.perf,
                TNPuv2=Wtuv.perf2,TNPuvk2=Wtuvk.perf2,TNPuvi2=Puvi.perf2)
  }
  output
}