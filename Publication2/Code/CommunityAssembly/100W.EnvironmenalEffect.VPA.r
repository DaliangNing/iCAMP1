rm(list=ls())
##########################################
# 1 # Files and parameters
##########################################
beta.file="/Code/AlphaBetaDiversity/Q100WSc.bMNTD.wt.csv"
prefixm="100WSc.bMNTD"
pcnm.thred="def" #20000 #"def" #20000
lnmsi=c(-0.584045068,-0.07292656,0.739541947,1.890063739,2.881942209,4.337192284,6.434493771)
pca.yn=TRUE #FALSE
rand.time=1000

if(pca.yn){pcann="pcay"}else{pcann="pcan"}
t00=Sys.time()
com.file="/Data/100WSc.OTUtable.csv"
grp.file="/Code/StressLevel/grouping.stressA.msi.csv"
env.file="/Data/100WSc.EnvFactor.csv"
xy.file="/Data/100WSc.XY.csv"
gdis.file="/Data/100WSc.GeoDist.csv"
geot.file="/Data/100WSc.GeoOthers.csv"

wd="/Code/CommunityAssembly"
prefix=paste0(prefixm,".Rep.MSI.PCNM",pcnm.thred,".",pcann)

##########################################
# 2 # load and match names
##########################################
source("/Code/handytool.r")
iwd(wd)
comm=t(lazyopen(com.file))
env=lazyopen(env.file)
xy=lazyopen(xy.file)
geot=lazyopen(geot.file)
gdis=lazyopen(gdis.file)
grpin=lazyopen(grp.file)
dis.in=lazyopen(beta.file)
dism=as.matrix(dis.in)

samp.check=match.name(rn.list = list(comm=comm,env=env,xy=xy,geot=geot,grpin=grpin),both.list = list(gdis=gdis,dism=dism))
comm=samp.check$comm
env=samp.check$env
xy=samp.check$xy
geot=samp.check$geot
gdis=samp.check$gdis
dism=samp.check$dism
grpin=samp.check$grpin
dim(comm);dim(env);dim(xy);dim(geot);dim(gdis);dim(dism);dim(grpin)

colnames(env)
envn.rm=c("AODC")
env.old=env
env=env.old[,which(!(colnames(env.old) %in% envn.rm)),drop=FALSE]
dim(env)

##########################################
# 3 # PCNM and PCA
##########################################
out=list()

### pcnm of distance
library(vegan)
if(pcnm.thred=="def")
{
  gdis.pcnm=vegan::pcnm(gdis)
}else{
  gdis.pcnm=vegan::pcnm(gdis,threshold = pcnm.thred)
}

gd.pcnm=gdis.pcnm$vectors[,1:10,drop=FALSE]
dim(gd.pcnm)
out$gdis.pcnm=gdis.pcnm
cumsum(gdis.pcnm$values/sum(gdis.pcnm$values[gdis.pcnm$values>0]))[1:11]

### combine and scale
allf=cbind(gd.pcnm,geot,env)
dim(allf)
out$all.factors=allf

comm=data.frame(comm,stringsAsFactors = FALSE)
dis=as.dist(dism)
out$com.dist=dis

pca.rda<-function(dis,factors,pca.yn=TRUE,...)
{
  ### PCA of all factors
  library(vegan)
  if(pca.yn)
  {
    fpca=rda(factors,scale = TRUE)
    fpc=data.frame(fpca$CA$u,stringsAsFactors = FALSE)
  }else{
    fpca=NULL
    fpc=data.frame(scale(factors),stringsAsFactors = FALSE)
  }
  if(ncol(fpc)>(nrow(as.matrix(dis))-2)){fpc=fpc[,1:(nrow(as.matrix(dis))-2)]}
  
  mod0<-capscale(dis~1,data=fpc,add=TRUE)
  mod1<-capscale(dis~.,data=fpc,add=TRUE)
  #ccast<-ordistep(mod0,scope=formula(mod1),direction = "both",perm.max=499,trace = FALSE,...)
  ccast<-ordiR2step(mod0,mod1,perm.max=499,trace = FALSE,...)
  scx<-summary(ccast)
  list(proportion=scx$constr.chi/scx$tot.chi,caps=ccast,pca=fpca)
}

t1=Sys.time()

grp=grpin
grp.lev=unique(grp[,1])


t1=Sys.time()

##########################################
# 4 # only env no distance
##########################################
egf=data.frame(cbind(env),stringsAsFactors = FALSE)
dim(egf)

eg.outs=list()
p.egfs=list()
eg.sts=list()
for(x in 1:100)
{
  res=list()
  message("------x=",x,". ",date())
  for(i in 1:length(grp.lev))
  {
    message("now egf i=",i," in ",length(grp.lev),". ",date())
    sampi=rownames(grp)[which(grp[,1]==grp.lev[i])]
    disi=as.dist(dism[match(sampi,rownames(dism)),match(sampi,colnames(dism))])
    fi=egf[match(sampi,rownames(egf)),]
    drdai=pca.rda(dis = disi,factors = fi,pca.yn = pca.yn,Pin=0.1,R2scope = FALSE)
    res[[i]]=drdai
  }
  names(res)=grp.lev
  #drda.all=pca.rda(dis = dis,factors = env,pca.yn = TRUE)
  eg.outs[[x]]=res
  (p.egfs[[x]]=sapply(res,function(v){out=v$proportion;if(length(out)==0){out=0};out}))
  corx=cor.test(p.egfs[[x]],lnmsi)
  eg.sts[[x]]=c(zeros=sum(p.egfs[[x]]==0),r=corx$estimate[[1]],p=corx$p.value)
}

eg.stsm=Reduce(rbind,eg.sts)
(eg.stsmnz=eg.stsm[which(eg.stsm[,1]==0),,drop=FALSE])
idz=which(eg.stsm[,1]==0)[which(eg.stsmnz[,2]==max(eg.stsmnz[,2]))]
eg.stsm[idz,,drop=FALSE]
p.egfs[idz]
eg.out=eg.outs[[idz[1]]]

##########################################
# 5 # PCNM
##########################################
gd.outs=list()
p.gdfs=list()
gd.sts=list()
for(x in 1:100)
{
  message("------x=",x,". ",date())
  res=list()
  
  for(i in 1:length(grp.lev))
  {
    message("now pcoa i=",i," in ",length(grp.lev),". ",date())
    sampi=rownames(grp)[which(grp[,1]==grp.lev[i])]
    disi=as.dist(dism[match(sampi,rownames(dism)),match(sampi,colnames(dism))])
    fi=gd.pcnm[match(sampi,rownames(gd.pcnm)),]
    drdai=pca.rda(dis = disi,factors = fi,pca.yn = FALSE,Pin=0.1,R2scope=FALSE)
    res[[i]]=drdai
  }
  names(res)=grp.lev
  #drda.all=pca.rda(dis = dis,factors = env,pca.yn = TRUE)
  gd.outs[[x]]=res
  (p.gdfs[[x]]=sapply(res,function(v){out=v$proportion;if(length(out)==0){out=0};out}))
  
  corx=cor.test(p.gdfs[[x]],lnmsi)
  (gd.sts[[x]]=c(zeros=sum(p.gdfs[[x]]==0),r=corx$estimate[[1]],p=corx$p.value))
}
gd.stsm=Reduce(rbind,gd.sts)
(gd.stsmnz=gd.stsm[which(gd.stsm[,1]<=2),,drop=FALSE])
idz=which(gd.stsm[,1]<=2)[which(gd.stsmnz[,2]==min(gd.stsmnz[,2]))]
gd.stsm[idz,,drop=FALSE]
p.gdfs[idz]
dist.out=gd.outs[[idz[1]]]

out$env.mods=eg.out
out$gdis.mods=dist.out


##########################################
# 5 # All
##########################################
cap.port<-function(cap){scx=summary(cap);out=scx$constr.chi/scx$tot.chi;if(length(out)==0){NA}else{out}}
sumv <- function(v)
{
  qt4 = quantile(v)
  names(qt4) = c("min", "q25", "median", "q75", "max")
  bxp = boxplot.stats(v)
  bxp4 = c(bxp$stats)
  names(bxp4) = c("LowerWhisker", "LowerHinge", "Median.plot", 
                  "HigherHinge", "HigherWhisker")
  if (length(bxp$out) > 0)
  {
    bxpout = bxp$out
    names(bxpout) = paste0("Outlier", 1:length(bxpout))
    bxp4 = c(bxp4, bxpout)
  }
  c(obs = v[[1]],mean = mean(v, na.rm = TRUE), stdev = sd(v, na.rm = TRUE), qt4, bxp4)
}
cbindl <- function(l)
{
  mlen = sapply(l, length)
  mn = names(l[[which(mlen == max(mlen))[1]]])
  ml = lapply(l,
              function(v)
              {
                out = c(v, rep(NA, max(mlen) - length(v)))
                names(out) = mn
                out
              })
  out = Reduce(cbind, ml)
  colnames(out) = names(l)
  out
}
###################################
bt.rec=list()
capsumm=list()
for(i in 1:length(grp.lev))
{
  message("Now partial i=",i, " in ",length(grp.lev),". ")
  pc.eg=substring(as.character(eg.out[[i]]$caps$call)[2],7)
  pc.dis=substring(as.character(dist.out[[i]]$caps$call)[2],7)
  
  sampi.obs=rownames(grp)[which(grp[,1]==grp.lev[i])]
  
  fi=egf[match(sampi.obs,rownames(egf)),]
  if(pca.yn)
  {
    fpca=rda(fi,scale = TRUE)
    fpc.obs=data.frame(fpca$CA$u,stringsAsFactors = FALSE)
  }else{
    fpca=NULL
    fpc.obs=data.frame(scale(fi),stringsAsFactors = FALSE)
  }
  captest<-function(sampi,pv.out=TRUE,...)
  {
    disi=as.dist(dism[match(sampi,rownames(dism)),match(sampi,colnames(dism))])
    xyi=gd.pcnm[match(sampi,rownames(gd.pcnm)),]
    fpc=fpc.obs[match(sampi,rownames(fpc.obs)),]
    fpcxyi=cbind(fpc,xyi)
    outn=c("both","env","dist","env_p_dist","dist_p_env")
    
    if(pc.eg=="1" & pc.dis=="1")
    {
      if(pv.out){pvix=rep(NA,length(outn))}
      proix=rep(0,length(outn))
    }else if(pc.eg!="1" & pc.dis=="1"){
      cap=capscale(formula = as.formula(paste("disi ~ ",pc.eg,sep = "")),data=fpc)
      if(pv.out){p.cap=min(sapply(1:100,function(x){anova(cap)$'Pr(>F)'[[1]]}));pvix=c(p.cap,p.cap,NA,p.cap,NA)}
      env.port=cap.port(cap)
      proix=c(env.port,env.port,0,env.port,0)
    }else if(pc.eg=="1" & pc.dis!="1"){
      capxyi=capscale(formula = as.formula(paste("disi ~ ",pc.dis,sep = "")),data=fpcxyi)
      if(pv.out){p.capxyi=min(sapply(1:100,function(x){anova(capxyi)$'Pr(>F)'[[1]]}));pvix=c(p.capxyi,NA,p.capxyi,NA,p.capxyi)}
      dis.port=cap.port(capxyi)
      proix=c(dis.port,0,dis.port,0,dis.port)
    }else{
      cap=capscale(formula = as.formula(paste("disi ~ ",pc.eg,sep = "")),data=fpc)
      cap.ed=capscale(formula = as.formula(paste("disi ~ ",pc.eg," + ",pc.dis,sep = "")),data=fpcxyi)
      cap.p=capscale(formula = as.formula(paste("disi ~ ",pc.eg," + Condition(",pc.dis,")",sep = "")),data=fpcxyi)
      capxyi=capscale(formula = as.formula(paste("disi ~ ",pc.dis,sep = "")),data=fpcxyi)
      capxyi.p=capscale(formula = as.formula(paste("disi ~ ",pc.dis," + Condition(",pc.eg,")",sep = "")),data=fpcxyi)
      
      if(pv.out)
      {
        p.cap=min(sapply(1:100,function(x){anova(cap)$'Pr(>F)'[[1]]}))
        p.caped=min(sapply(1:100,function(x){anova(cap.ed)$'Pr(>F)'[[1]]}))
        p.capp=min(sapply(1:100,function(x){anova(cap.p)$'Pr(>F)'[[1]]}))
        p.capxyi=min(sapply(1:100,function(x){anova(capxyi)$'Pr(>F)'[[1]]}))
        p.capxyip=min(sapply(1:100,function(x){anova(capxyi.p)$'Pr(>F)'[[1]]}))
        pvix=c(p.caped,p.cap,p.capxyi,p.capp,p.capxyip)
      }
      proix=c(cap.port(cap.ed),cap.port(cap),cap.port(capxyi),cap.port(cap.p),cap.port(capxyi.p))
    }
    if(pv.out){outix=data.frame(item=outn,proprotion=proix,pvalue=pvix,stringsAsFactors = FALSE)}else{
      outix=data.frame(item=outn,proprotion=proix,stringsAsFactors = FALSE)
    }
    outix
  }
  
  message("Calculating observed pro and p. ",date())
  propv.obsi=captest(sampi=sampi.obs,pv.out = TRUE)
  
  pro.bti=sapply(1:rand.time,
                 function(x)
                 {
                   message("Bootstrapping i=",i," x=",x,". ",date())
                   sampi=sample(sampi.obs,length(sampi.obs),replace = TRUE)
                   if(length(unique(sampi))<=3)
                   {
                     t=1
                     while(t<50 & length(unique(sampi))<=3)
                     {
                       sampi=sample(sampi.obs,length(sampi.obs),replace = TRUE)
                       t=t+1
                     }
                   }
                   captest(sampi=sampi,pv.out = FALSE)[,2]
                 })
  pro.btix=cbind(propv.obsi[,2],pro.bti)
  probtmi=cbindl(lapply(1:nrow(pro.btix),function(x){sumv(pro.btix[x,])}))
  bt.rec[[i]]=pro.btix
  capsumm[[i]]=data.frame(group=grp.lev[i],index=c(colnames(propv.obsi)[2:3],rownames(probtmi)),
                  rbind(t(propv.obsi[,2:3]),probtmi),stringsAsFactors = FALSE)
  rownames(capsumm[[i]])=c()
  colnames(capsumm[[i]])[3:7]=propv.obsi[,1]
}

capsum=Reduce(rbind,capsumm)
head(capsum)
out$summary=capsum
names(bt.rec)=grp.lev
out$boot.detail=bt.rec

iwd(wd.save)
save(out,file=paste(prefix,".detail.rda",sep = ""))
save.file(capsum,prefix = prefix,filename = "summary")

# end
