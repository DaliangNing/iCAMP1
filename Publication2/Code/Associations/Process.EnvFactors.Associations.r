rm(list=ls())
env.file="/Data/100WSc.EnvFactor.csv"
res.file="/Code/CommunityAssembly/iCAMP.Results/Q100W.ProcessImp.Pairwise.bNRIi.RCa.csv"
dis.file="/Data/100WSc.GeoDist.csv"
geo2.file="/Data/100WSc.GeoOthers.csv"
r.cut=0.1
p.cut=0.05

envLN=TRUE
prefix="EnvlogPIcrl"
source("/Code/handytool.r")
wd="/Code/Associations"
####################################
# 1 # load data
####################################
# load process importance
library(coda.base)
wd=iwd(wd)
resm=lazyopen(res.file)
dim(resm)
samps=unique(as.vector(as.matrix(resm[,1:2])))
if(prefix=="EnvlogPIcrl")
{
  library(coda.base)
  resm2=resm[,3:ncol(resm),drop=FALSE]
  resm2[resm2==0]=min(resm2[resm2>0])/20
  min(resm2)
  resm3=coda.base::coordinates(resm2,basis = "clr")
  colnames(resm3)=colnames(resm2)
  head(resm3)
  
  resm.old=resm
  resm=cbind(resm[,1:2,drop=FALSE],resm3)
}
# load env and transform data
env.in=lazyopen(env.file)
sum(!(rownames(env.in) %in% samps))
sum(is.na(match(samps,rownames(env.in))))
env.in=env.in[match(samps,rownames(env.in)),,drop=FALSE]
dim(env.in)
colnames(env.in)
env=env.in[,2:ncol(env.in),drop=FALSE]
colnames(env)
minenv=sapply(1:ncol(env),function(i){min(env[,i])})
id.neg=which(minenv<0)
envtrs=env
envtrs[,id.neg]=envtrs[,id.neg]-matrix(minenv[id.neg]*1.05,nrow = nrow(envtrs),ncol=length(id.neg),byrow = TRUE) #(minenv[id.neg]*1.05)
minpenv=sapply(1:ncol(env),function(i){min(env[,i][env[,i]>0])})
id.pos=which(minenv>=0 & (colnames(env)!="pH"))
colnames(env)[id.pos]
envtrs[,id.pos]=envtrs[,id.pos]+matrix(minpenv[id.pos]/20,nrow=nrow(envtrs),ncol=length(id.pos),byrow = TRUE)
min(sapply(1:ncol(envtrs),function(i){min(envtrs[,i])}))
envlog=log(envtrs)
head(envlog)
envlog[,"pH"]=env[,"pH"]
sum(is.na(envlog));min(envlog);max(envlog)
envsc=scale(env)
envlns=scale(envlog)

# load geographic 
gdis=lazyopen(dis.file)
sum(is.na(match(samps,rownames(gdis))))
gdis=gdis[match(samps,rownames(gdis)),match(samps,colnames(gdis))]
dim(gdis)
gdism=dist.3col(gdis)
dim(gdism)
colnames(gdism)[3]="Geodist"

geom=lazyopen(geo2.file)
sum(is.na(match(samps,rownames(geom))))
geom=geom[match(samps,rownames(geom)),,drop=FALSE]
dim(geom)
min(geom)
geosc=scale(geom)
geolog=log(geom)
geologsc=scale(geolog)
##############################
if(envLN)
{
  ckid=match.name(rn.list=list(envlns,geologsc))
  gdism[,3]=log(gdism[,3])
}else{
  ckid=match.name(rn.list=list(envsc,geosc))
}
###############################
fv=cbind(ckid$rn.list.1,ckid$rn.list.2)
colnames(fv)
ncol(fv);length(unique(colnames(fv)))

fvdm1=sapply(1:ncol(fv),
             function(i)
             {
               dist.3col(dist(fv[,i]))[,3]
             })
dim(fvdm1)
colnames(fvdm1)=colnames(fv)
fvdm1[1:5,1:5]
fvdm=cbind(dist.3col(dist(fv[,1]))[,1:2,drop=FALSE],fvdm1)
fvdm[1:5,1:5]

ckid2=match.2col(check.list = list(a=resm,b=gdism,c=fvdm))
fm=cbind(ckid2$a[,3:ncol(ckid2$a),drop=FALSE],ckid2$b[,3:ncol(ckid2$b),drop=FALSE],ckid2$c[,3:ncol(ckid2$c),drop=FALSE])
dim(fm)
fma=cbind(ckid2$a[,1:2,drop=FALSE],fm)
ncol(fm);length(unique(colnames(fm)))
sum(!(colnames(fv) %in% colnames(fm)))

fv[1:5,1:5]
save.file(fv,prefix=prefix,filename = "EnvGeoScaled")
fm[1:5,1:10]
save.file(data.frame(ckid2$a[,1:2,drop=FALSE],fm,stringsAsFactors = FALSE),
          prefix=prefix,filename = "ProcEnvGeoScaled.Dist")
######################################
# 2 # association matrix
######################################
# 2.1 # cross validated linear model for vectors
######################################
source("/Code/Associations/cv.lm.r")
cor.meth="pearson"
fvr=matrix(1,nrow=ncol(fv),ncol=ncol(fv))
rownames(fvr)<-colnames(fvr)<-colnames(fv)
fvp=fvr;fvp[]=0
fvrpout=list()
k=1
for(i in 1:(ncol(fv)-1))
{
  fvi=fv[,i]
  for(j in (i+1):ncol(fv))
  {
    message("--- i=",i," j=",j," k=",k," ",date())
    fvj=fv[,j]
    dataij=data.frame(x=fvi,y=fvj,stringsAsFactors = FALSE)
    cvij=cv.lm(y~x,data = dataij, method = "monte", mt.nv = 0.1, mt.rand = 1000,perm.test.rand = 1000)
    fvrpout[[k]]=data.frame(factor1=colnames(fv)[i],factor2=colnames(fv)[j],
                   R2.obs=cvij$R2[1],R2.adj=cvij$R2[2],R2cv=cvij$R2cv,P.Ftest=cvij$p.values['lm.p.Ftest'],
                   P.perm=cvij$p.values['lm.p.perm'],stringsAsFactors = FALSE)
    k=k+1
    fvp[i,j]<-fvp[j,i]<-cvij$p.values['lm.p.perm']
    fvr[i,j]<-fvr[j,i]<-cvij$R2cv
  }
}
fvrpoutm=data.frame(data.table::rbindlist(fvrpout),stringsAsFactors = FALSE)
fvrpoutm$P.adj=p.adjust(fvrpoutm$P.perm,method = "fdr")
head(fvrpoutm)
dim(fvrpoutm)
save.file(fvrpoutm,prefix = prefix, filename = "EnvGeoCor")

fvpa=fvp
for(i in 1:(ncol(fv)-1))
{
  for(j in (i+1):ncol(fv))
  {
    fvpa[i,j]<-fvpa[j,i]<-fvrpoutm$P.adj[which(fvrpoutm$factor1==colnames(fv)[i] & fvrpoutm$factor2==colnames(fv)[j])]
  }
}
fvcor=fvr
fvcor[fvr<(r.cut^2)]=NA
fvcor[fvpa>p.cut]=NA
sum(is.na(fvcor))/sum(!is.na(fvcor))

######################################
# 2.2 # cross validated Mantel test for matrixes
######################################
source("/Code/Associations/mantel.n.r")
source("/Code/Associations/mantel.cv.r")

fmr=matrix(NA,nrow=ncol(fm),ncol=ncol(fm))
diag(fmr)=1
rownames(fmr)<-colnames(fmr)<-colnames(fm)
fmp=fmr
diag(fmp)=0
fmrpout=list()
k=1
for(i in 1:(ncol(fm)-1))
{
  message("---i=",i,". ",date())
  fmi=cbind(fma[,1:2,drop=FALSE],fm[,i,drop=FALSE])
  fmdi=col3.dist(fmi,to.dist = FALSE)
  for(j in (i+1):ncol(fm))
  {
    if(sum(colnames(fm)[c(i,j)] %in% colnames(fv))<2)
    {
      message("---i=",i," j=",j," k=",k,". ",date())
      fmj=cbind(fma[,1:2,drop=FALSE],fm[,j,drop=FALSE])
      fmdj=col3.dist(fmj,to.dist = FALSE)
      corij=mantel.n(xmatrix = fmdi,ymatrix = fmdj,method = cor.meth)
      mtcvij=mantel.cv(xmatrix = fmdi,ymatrix = fmdj,method = 'monte',mt.nv = 0.1, mt.rand = 200)

      fmrpout[[k]]=data.frame(factor1=colnames(fm)[i],factor2=colnames(fm)[j],
                     r=corij$statistic, R2.obs=mtcvij$R2[[1]],
                     R2.adj=mtcvij$R2[[2]],R2cv=mtcvij$R2cv,
                     P.Mantel=corij$signif,stringsAsFactors = FALSE)
      k=k+1
      fmp[i,j]<-fmp[j,i]<-corij$signif
      fmr[i,j]<-fmr[j,i]<-mtcvij$R2cv
    }
  }
}

fmrpoutm=data.table::rbindlist(fmrpout)
fmrpoutm$P.adj=p.adjust(fmrpoutm$P.Mantel,method = "fdr")
dim(fmrpoutm)
head(fmrpoutm)
save.file(fmrpoutm, prefix = prefix, filename = "ProcEnvCor")

fmpa=fmp
for(i in 1:(ncol(fm)-1))
{
  for(j in (i+1):ncol(fm))
  {
    if(sum(colnames(fm)[c(i,j)] %in% colnames(fv))<2)
    {
      fmpa[i,j]<-fmpa[j,i]<-fmrpoutm$P.adj[which(fmrpoutm$factor1==colnames(fm)[i] & fmrpoutm$factor2==colnames(fm)[j])]
    }
  }
}
#####################

fmr[match(rownames(fvr),rownames(fmr)),match(colnames(fvr),colnames(fmr))]=fvr
fmp[match(rownames(fvp),rownames(fmp)),match(colnames(fvp),colnames(fmp))]=fvp
fmpa[match(rownames(fvpa),rownames(fmpa)),match(colnames(fvpa),colnames(fmpa))]=fvpa
fmcor=fmr
fmcor[1:5,1:5]
save.file(fmcor,prefix=prefix,filename = "CV.R2")
save.file(fmp,prefix=prefix,filename = "Pvalue")
save.file(fmpa,prefix=prefix,filename = "P.adj.value")

fmcor[fmr<(r.cut^2)]=NA
fmcor[fmpa>p.cut]=NA
sum(is.na(fmcor))/sum(!is.na(fmcor))
save.file(fmcor,prefix=prefix,filename = "CV.R2.sig")
save(fmcor,file = paste0(prefix,".CVSigCor.rda"))

# end
