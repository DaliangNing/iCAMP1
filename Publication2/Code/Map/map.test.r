# contact # Daliang Ning (ningdaliang@ou.edu)

rm(list=ls())
wd="/Code/Map"
env.file="/Data/100WSc.Env.csv"
ref.file="/Data/100WSc.Env.Reference.csv"
xy.file="/Data/100WSc.XY.csv"

#########################################
# 1 # load package and data, and basic setting for cross validation
#########################################
source("/Code/handytool.r")
wd=iwd(wd)
env=lazyopen(env.file)
ref=lazyopen(ref.file)
xy=lazyopen(xy.file)
envr=env[,match(rownames(ref),colnames(env))]
dim(env);dim(ref);dim(envr)

sampc=match.name(rn.list = list(env,xy))

SI=t(t(envr)/ref[,1])
MSI=apply(SI, 1, max)
MSI.env=colnames(envr)[apply(SI,1,which.max)]

nv=0.1 # ratio of test data in cross validation
rand=100
n.test=round(nrow(env)*nv)
n.train=nrow(env)-n.test
perm.n = permute::shuffleSet(nrow(env), 1000)
library(akima)
x.num=100
y.num=300

seta=(-52/180)*pi
x.old=xy$splane_eas
y.old=xy$splane_nor
x.old.m=x.old-min(x.old)
y.old.m=y.old-min(y.old)
x=((y.old.m*tan(seta))+x.old.m)*cos(seta)
y=((y.old.m/tan(seta))-x.old.m)*sin(seta)

(max(x)-min(x))/(max(y)-min(y))

xo=seq(min(x), max(x), length = x.num)
yo=seq(min(y), max(y), length = y.num)
xo1=sort(unique(c(xo,x)))
yo1=sort(unique(c(yo,y)))

#############################
# handy functions
#############################
lnx<-function(v){minv=min(v,na.rm = TRUE);if(minv<0){v=v-minv};out=log(v);if(sum(v==0,na.rm=TRUE)>0){out[v==0]=min(density(log(v[v>0]))$x)};out}

#############################
# 2 # MSI: preliminary test
#############################

SIZ=sapply(1:ncol(SI),
           function(i)
           {
             int=interp(x,y,SI[,i],xo=xo1,yo=yo1,duplicate = "mean")
             int$z
           },simplify = "array")
MSIZ=apply(SIZ, c(1,2), max)
MSIZ.obs=sapply(1:length(x),function(i){MSIZ[match(x[i],xo1),match(y[i],yo1)]})
perm.n = permute::shuffleSet(length(MSIZ.obs), 1000)
MSIZ.perm=sapply(1:nrow(perm.n),function(i){MSIZ.obs[perm.n[i,]]})
colnames(MSIZ.perm)=paste0("Rand",1:ncol(MSIZ.perm))

n.test=round(length(x)*nv)
n.train=length(x)-n.test
trsq=seq(from=1,to=rand,by=2)

datt=lapply(1:rand,
            function(j)
            {
              message("j=",j,". ",date())
              id.train=sort(sample(length(x),size = n.train))
              id.test=which(!((1:length(x)) %in% id.train))
              SIZt=sapply(1:ncol(SI),
                          function(i)
                          {
                            int=interp(x[id.train],y[id.train],SI[id.train,i],xo=xo1,yo=yo1,duplicate = "mean")
                            int$z
                          },simplify = "array")
              MSIZt=apply(SIZt, c(1,2), max, na.rm=TRUE)
              MSIZt[which(MSIZt==(-Inf),arr.ind = TRUE)]=NA
              z.test=sapply(id.test,function(i){MSIZt[match(x[i],xo1),match(y[i],yo1)]})
              z.obs=MSIZ.obs[id.test]
              z.perm=MSIZ.perm[id.test,,drop=FALSE]
              idi=which(!is.na(z.test))
              if(length(idi)>0)
              {
                out=data.frame(obs=z.obs[idi],test=z.test[idi],z.perm[idi,,drop=FALSE],stringsAsFactors = FALSE)
              }else{
                out=NULL
              }
              out
            })
datm=Reduce(rbind,datt)
rownames(datm)=c()
R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
R2cv=R2fun(yv=datm$obs,fv=datm$test)
R2cvln=R2fun(yv=log(datm$obs),fv=log(datm$test))
c(R2cv,R2cvln)
dim(datm)
EPS <- sqrt(.Machine$double.eps)
R2.perm=sapply(3:ncol(datm),
               function(i)
               {
                 c(R2fun(yv=datm[,i],fv=datm$test),
                   R2fun(yv=log(datm[,i]),fv=log(datm$test)))
               })
pv=(rowSums(R2.perm>=(matrix(c(R2cv,R2cvln),nrow = nrow(R2.perm),ncol = ncol(R2.perm))-EPS))+1)/(nrow(perm.n)+1)
out=data.frame(Index="MSI",Transform=c("No","Log"),R2=c(R2cv,R2cvln),P=pv)
save.file(out,filename = "MSI.map.R2.p")

#####################################
# 3 # Processes scores
#####################################
isc.res=lazyopen("/Code/CommunityAssembly/iCAMP.Results/Q100W.iCAMP.detail.rda")
isc=isc.res$processes$bNRIiRCa
head(isc)

samps=rownames(env)

x.num=100
y.num=300

psc=t(sapply(1:length(samps),
             function(i)
             {
               idi=which(isc[,1]==samps[i] | isc[,2]==samps[i])
               psci=isc[idi,3:7]
               sampi=isc[idi,1]
               sampi[sampi==samps[i]]=isc[idi,2][sampi==samps[i]]
               sapply(1:ncol(psci),
                      function(j)
                      {
                        intij=interp(x[rownames(xy)!=samps[i]],y[rownames(xy)!=samps[i]],
                                     psci[,j],xo=seq(min(x), max(x), length = x.num),
                                     yo=seq(min(y), max(y), length = x.num),duplicate = "mean")
                        mean(intij$z,na.rm = TRUE)
                      })
             }))
rownames(psc)=samps
colnames(psc)=colnames(isc)[3:7]
head(psc)
sampc=match.name(name.check = rownames(env),rn.list = list(psc=psc))


#####################################
# 4 # Other data preparation
#####################################
x.num=100
y.num=300
xo=seq(min(x), max(x), length = x.num)
yo=seq(min(y), max(y), length = y.num)
xo1=sort(unique(c(xo,x)))
yo1=sort(unique(c(yo,y)))

envn.rm=c("Date_Collected","H","OH")
env2=env[,which(!(colnames(env) %in% envn.rm)),drop=FALSE]
lnx<-function(v){minv=min(v,na.rm = TRUE);if(minv<0){v=v-minv};out=log(v);if(sum(v==0,na.rm = TRUE)>0){out[v==0]=min(density(log(v[v>0]))$x)};out}
env2ln=data.frame(sapply(1:ncol(env2),function(i){lnx(env2[,i])}))
rownames(env2ln)=rownames(env2)
colnames(env2ln)=colnames(env2)
env2ln$pH=env2$pH

comm=t(lazyopen("/Data/100WSc.OTUtable.csv"))
dim(comm)

#################################
# 5 # functions
#################################

R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
lnx<-function(v,detailout=FALSE)
{
  minv=min(v,na.rm = TRUE)
  if(minv<0){v=v-minv}else{minv=0}
  out=log(v)
  if(sum(v==0,na.rm = TRUE)>0){zerov=min(density(log(v[which(v>0)]))$x);out[v==0]=zerov}else{zerov=NULL}
  if(detailout){out=list(result=out,zerov=zerov,minv=minv)}
  out
}
lnxrev<-function(v,minv=0,zerov=NULL){EPS <- sqrt(.Machine$double.eps); out=exp(v)+minv; if(!is.null(zerov)){out[which(abs(v-zerov)<=EPS)]=minv}; out}
zintf<-function(XYZm,id.train,...)
{
  id.test=which(!((1:nrow(XYZm)) %in% id.train))
  idxo=match(XYZm$x,xo1)
  idyo=match(XYZm$y,yo1)
  # method 1 # classic bilinear interpolation (CBI)
  intxyzf<-function(xf,yf,zf,xof,yof,...)
  {
    intkf=interp(x=xf[id.train],y=yf[id.train],z=zf[id.train],xo=xof,yo=yof,duplicate = "mean",linear = TRUE)
    sapply(id.test,function(i){intkf$z[idxo[i],idyo[i]]})
  }
  zt.CBI=intxyzf(xf=XYZm$x,yf=XYZm$y,zf=XYZm$z,xof=xo1,yof=yo1)
  zln=lnx(XYZm$z,T)
  zt.CBIln=lnxrev(intxyzf(xf=XYZm$x,yf=XYZm$y,zf=zln$result,xof=xo1,yof=yo1),minv = zln$minv,zerov = zln$zerov)
  out=data.frame(ID=id.test, obs=XYZm$z[id.test], pred.CBI=zt.CBI, pred.CBILn=zt.CBIln,
             stringsAsFactors = FALSE)
  out
}
R2cvf<-function(z.obs,znamei,rand=100,...)
{
  n.test=round(length(x)*nv)
  n.train=length(x)-n.test
  trsq=seq(from=1,to=rand,by=10)
  datt=lapply(1:rand,
              function(j)
              {
                if(j %in% trsq) {message("zname=",znamei,"  j=",j,". ",date())}
                id.train=sort(sample(length(x),size = n.train))
                zintf(XYZm=data.frame(x=x,y=y,z=z.obs),id.train=id.train)
              })
  datm=Reduce(rbind,datt)
  rownames(datm)=c()
  datm=datm[which(rowSums(is.na(datm))==0),,drop=FALSE]
  R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
  R2cv=sapply(3:ncol(datm),function(k){R2fun(datm$obs,datm[,k])})
  names(R2cv)=colnames(datm)[3:ncol(datm)]
  R2cv
}

#################################
# 5 # Find reliably interpolated variables: env.sig
#################################
R2env=list()
asv.cut=45
comm2=comm[,which(colSums(comm>0)>asv.cut),drop=FALSE]
dim(comm2)
sampc=match.name(rn.list = list(comm2,psc,env2))
envm=cbind(MSI=MSI,env2,comm2)
envm.ln=cbind(MSI=log(MSI),env2ln,log(comm2+1))

for(i in 1:ncol(envm))
{
  message("i=",i," in ",ncol(envm),". ",date())
  R2env[[i]]=R2cvf(z.obs=envm[,i],znamei=colnames(envm)[i])
  print(R2cv[which.max(R2cv)])
}
names(R2env)=colnames(envm)
save(R2env,file = "Factors.EachMap.R2cv.rda")
R2envm=Reduce(rbind,R2env)
rownames(R2envm)=colnames(envm)  
R2envm=R2envm[,1:2,drop=FALSE] 
R2mx=apply(R2envm,1,max)
R2mx.method=apply(R2envm,1,function(v){colnames(R2envm)[which.max(v)]})
id.sig=which(R2mx>0.25)
envname.sig=rownames(R2envm)[id.sig]
method.sig=c(R2mx.method[id.sig],trans_x="pred.CBI",trans_y="pred.CBI")
env.sig=data.frame(envm[,match(envname.sig,colnames(envm)),drop=FALSE],
                   trans_x=(x-min(x)),trans_y=(y-min(y)),check.names = FALSE)
env.sig.ln=data.frame(envm.ln[,match(envname.sig,colnames(envm.ln)),drop=FALSE],
                      trans_x=lnx(x-min(x)),trans_y=lnx(y-min(y)),check.names = FALSE)
env.sig.mm=data.frame(env.sig,env.sig.ln,check.names = FALSE)
colnames(env.sig.mm)=c(colnames(env.sig),paste0(colnames(env.sig.ln),".Ln"))
dim(env.sig);dim(env.sig.ln);dim(env.sig.mm);length(method.sig)
sig.sum=data.frame(Factors=names(R2mx[id.sig]),Method.Sig=R2mx.method[id.sig],R2.sig=R2mx[id.sig])
save.file(sig.sum,filename = "Factors.EachMap.R2cvMax.EnvSig.summary2")

#################################
# 6 # use env.sig to predict a certain variable zi
#################################
# 6.1 # functions
#################################
mapintf<-function(XYZ.train,XY.test,method,...)
{
  idxo=match(XY.test$x,xo1)
  idyo=match(XY.test$y,yo1)
  idxot=match(XYZ.train$x,xo1)
  idyot=match(XYZ.train$y,yo1)
  if(substring(method,nchar(method)-1)=="Ln"){zln=lnx(XYZ.train$z,T)}
  intxyzf<-function(xtr,ytr,ztr,xof,yof,...)
  {
    intkf=interp(x=xtr,y=ytr,z=ztr,xo=xof,yo=yof,duplicate = "mean",linear = TRUE)
    sapply(1:length(idxo),function(i){intkf$z[idxo[i],idyo[i]]})
  }
  dwmif<-function(xtr,ytr,ztr,xts,yts,method=c("Rev","LnRev","SqRev"),LNZ=FALSE,add.lm=FALSE,...)
  {
    dwmf<-function(xtarget,ytarget,XYZpm,method=c("Rev","LnRev","SqRev"))
    {
      EPS <- sqrt(.Machine$double.eps)
      distv=(((XYZpm$x-xtarget)^2)+((XYZpm$y-ytarget)^2))^0.5
      id.zero=which(distv<=EPS)
      if(length(id.zero)>0){out=mean(XYZpm$z[id.zero],na.rm=TRUE)}else{
        if(method[1]=="Rev"){distwt=1/distv}
        if(method[1]=="LnRev"){distwt=1/(log(distv))}
        if(method[1]=="SqRev"){distwt=1/(distv^2)}
        out=sum(XYZpm$z*(distwt/sum(distwt)))
      }
      out
    }
    z.dwma=sapply(1:length(xts),function(k){dwmf(xtarget = xts[k],ytarget = yts[k],XYZpm = data.frame(x=xtr,y=ytr,z=ztr),method = method)})
    if(add.lm)
    {
      z.dwm=sapply(1:length(ztr),function(k){dwmf(xtarget = xtr[k],ytarget = ytr[k],XYZpm = data.frame(x=xtr[-k],y=ytr[-k],z=ztr[-k]),method=method)})
      zdwlm=lm(ztr~z.dwm)
      out1=as.vector(predict(zdwlm,newdata=data.frame(z.dwm=z.dwma)))
    }else{out1=z.dwma}     
    if(LNZ)
    {
      outi=lnxrev(out1,minv = zln$minv,zerov = zln$zerov)
    }else{outi=out1}
    outi
  }
  xoln=lnx(xo1);yoln=lnx(yo1)
  xosq=(xo1-min(xo1))^2;yosq=(yo1-min(yo1))^2
  if(method=="pred.CBI"){out=intxyzf(xtr = XYZ.train$x,ytr = XYZ.train$y,ztr = XYZ.train$z,xof = xo1, yof = yo1)}
  if(method=="pred.CBILn"){out=lnxrev(intxyzf(xtr = XYZ.train$x,ytr = XYZ.train$y,ztr = zln$result,xof = xo1, yof = yo1),minv = zln$minv,zerov = zln$zerov)}
  if(method=="pred.LCBI"){out=intxyzf(xtr = xoln[idxot],ytr = yoln[idyot],ztr = XYZ.train$z,xof = xoln, yof = yoln)}
  if(method=="pred.LCBILn"){out=lnxrev(intxyzf(xtr = xoln[idxot],ytr = yoln[idyot],ztr = zln$result,xof = xoln, yof = yoln),minv = zln$minv,zerov = zln$zerov)}
  if(method=="pred.SCBI"){out=intxyzf(xtr = xosq[idxot],ytr = yosq[idyot],ztr = XYZ.train$z,xof = xosq, yof = yosq)}
  if(method=="pred.SCBILn"){out=lnxrev(intxyzf(xtr = xosq[idxot],ytr = yosq[idyot],ztr = zln$result,xof = xosq, yof = yosq),minv = zln$minv,zerov = zln$zerov)}
  if(grepl(pattern = "pred.*DWMI",method))
  {#print("find DWMI")
    if(grepl(pattern = "pred\\.[R,L,S]DWMI",method))
    {
      xtr=XYZ.train$x;ytr=XYZ.train$y;xts=XY.test$x;yts=XY.test$y
    }else if(grepl(pattern = "pred\\.L[R,L,S]DWMI",method)){
      xtr=xoln[idxot];ytr=yoln[idyot];xts=xoln[idxo];yts=yoln[idyo]
    }else if(grepl(pattern = "pred\\.S[R,L,S]DWMI",method)){
      xtr=xosq[idxot];ytr=yosq[idyot];xts=xosq[idxo];yts=yosq[idyo]
    }
    if(grepl(pattern = "Ln$",method)){ztr=zln$result;LNZ=TRUE}else{ztr=XYZ.train$z;LNZ=FALSE}
    if(grepl(pattern = "DWMI2",method)){add.lm=TRUE}else{add.lm=FALSE}
    if(grepl(pattern = "RDWMI",method)){method="Rev"}else if(grepl(pattern = "LDWMI",method)){method="LnRev"}else if(grepl(pattern = "SDWMI",method)){method="SqRev"}
    out=dwmif(xtr=xtr,ytr=ytr,ztr=ztr,xts=xts,yts=yts,method=method,LNZ=LNZ,add.lm=add.lm)
  }
  out
}
modr2f<-function(yyi,xxmi,mod,cvmethod="monte")
{
  source("C:/Users/Daliang/Dropbox/Data_analysis_service/ieg.R.package/DevelopingFunctions/Daliang.CrossValidateLm/cv.pls.rf.r")
  if(mod=="lm"){modxi=cv.LM(Y=yyi,X=xxmi,cvmethod=cvmethod)}
  if(mod=="pls"){modxi=cv.PLS(Y=yyi,X=xxmi,cvmethod=cvmethod)}
  if(mod=="RF"){modxi=cv.RF(Y=yyi,X=xxmi,cvmethod=cvmethod)}
  modxi
}
lgst<-function(v)
{
  denv=density(v[!is.na(v)])
  minv=min(denv$x)
  maxv=max(denv$x)
  vt=(v-minv)/(maxv-minv)
  list(result=log(vt/(1-vt)),minv=minv,maxv=maxv)
}
lgstrev<-function(lgstv,minv,maxv)
{
  vt=1/(exp(-lgstv)+1)
  (vt*(maxv-minv))+minv
}
detrn.lgst<-function(detailm,minv,maxv)
{
  obsn=lgstrev(detailm$obs,minv = minv,maxv = maxv)
  testn=lgstrev(detailm$test,minv = minv, maxv = maxv)
  list(R2=R2fun(obsn,testn),detailn=data.frame(ID=detailm$ID,obs=obsn,test=testn))
}
modstep.cv<-function(yy,xx,n.limit,mod="lm")
{
  id.sel=NULL
  id.cad=1:ncol(xx)
  R2stp=list()
  R2stp[[1]]=-Inf
  derR=1
  k=1
  yyt=lgst(yy)
  
  while(derR>0.001 & length(id.sel)<n.limit)
  {
    R2xi=sapply(id.cad,
                function(xxi)
                {
                  modr2f(yyi=yyt$result,xxmi=xx[,c(id.sel,xxi),drop=FALSE],mod=mod,cvmethod="k-fold")$R2cv
                })
    if(max(R2xi)>R2stp[[k]])
    {
      id.sel=c(id.sel,id.cad[which.max(R2xi)])
      id.cad=which(!((1:ncol(xx)) %in% id.sel))
      k=k+1
      R2stp[[k]]=max(R2xi)
      derR=R2stp[[k]]-R2stp[[k-1]]
      message("Now R2=",max(R2xi)," id.sel number=",length(id.sel),". ",date())
    }else{
      message("Now R2 decrease, break, id.sel number=",length(id.sel),". ",date())
      derR=max(R2xi)-R2stp[[k]]
    }
  }
  
  modop=modr2f(yyi=yyt$result,xxmi=xx[,id.sel,drop=FALSE],mod=mod,cvmethod="k-fold")
  modop2=modr2f(yyi=yyt$result,xxmi=xx[,id.sel,drop=FALSE],mod=mod,cvmethod="monte")
  det.kfold=detrn.lgst(detailm = modop$detail,minv = yyt$minv,maxv = yyt$maxv)
  det.mod=detrn.lgst(detailm = modop$detail.mod,minv = yyt$minv,maxv = yyt$maxv)
  det.monte=detrn.lgst(detailm = modop2$detail,minv = yyt$minv,maxv = yyt$maxv)
  list(R2cv=det.kfold$R2,R2cv.monte=det.monte$R2,R2=det.mod$R2,
       detail=det.kfold$detailn,detail.monte=det.monte$detailn,detail.mod=det.mod$detailn,
       Y=yy,X=xx[,id.sel,drop=FALSE],Y.lgst=yyt)
}
modstep2.cv<-function(yy,xx,xxln,n.limit,mod="lm")
{
  if(sum(colnames(xx)!=colnames(xxln))>0){stop("colnames not match")}
  id.sel1=NULL
  id.sel2=NULL
  id.cad=1:ncol(xx)
  R2stp=list()
  R2stp[[1]]=-Inf
  derR=1
  k=1
  yyt=lgst(yy)
  
  if(mod=="lm"){n.limit.use=n.limit}else{n.limit.use=Inf}
  while(derR>0.001 & length(c(id.sel1,id.sel2))<n.limit.use)
  {
    R2xi=sapply(id.cad,
                function(xxi)
                {
                  c(modr2f(yyi=yyt$result,xxmi=data.frame(xx[,c(id.sel1,xxi),drop=FALSE],xxln[,id.sel2,drop=FALSE],check.names = FALSE),mod=mod,cvmethod="k-fold")$R2cv,
                    modr2f(yyi=yyt$result,xxmi=data.frame(xx[,id.sel1,drop=FALSE],xxln[,c(id.sel2,xxi),drop=FALSE],check.names = FALSE),mod=mod,cvmethod="k-fold")$R2cv)
                })
    R2xi.max=max(R2xi,na.rm = TRUE)
    if(R2xi.max>R2stp[[k]])
    {
      idr2m=which(R2xi==R2xi.max,arr.ind = TRUE)
      id.sela=id.cad[idr2m[1,2]]
      if(idr2m[1,1]==1){id.sel1=c(id.sel1,id.sela)}
      if(idr2m[1,1]==2){id.sel2=c(id.sel2,id.sela)}
      id.cad=which(!((1:ncol(xx)) %in% c(id.sel1,id.sel2)))
      k=k+1
      R2stp[[k]]=max(R2xi,na.rm=TRUE)
      derR=R2stp[[k]]-R2stp[[k-1]]
      message("Now R2=",R2stp[[k]]," id.sel number=",length(c(id.sel1,id.sel2)),". ",date())
    }else{
      message("Now R2 decrease, break, id.sel number=",length(c(id.sel1,id.sel2)),". ",date())
      derR=R2xi.max-R2stp[[k]]
    }
  }
  modop=modr2f(yyi=yy,xxmi=data.frame(xx[,id.sel1,drop=FALSE],xxln[,id.sel2,drop=FALSE],check.names = FALSE),mod=mod,cvmethod="k-fold")
  modop2=modr2f(yyi=yy,xxmi=data.frame(xx[,id.sel1,drop=FALSE],xxln[,id.sel2,drop=FALSE],check.names = FALSE),mod=mod,cvmethod="monte")
  det.kfold=detrn.lgst(detailm = modop$detail,minv = yyt$minv,maxv = yyt$maxv)
  det.mod=detrn.lgst(detailm = modop$detail.mod,minv = yyt$minv,maxv = yyt$maxv)
  det.monte=detrn.lgst(detailm = modop2$detail,minv = yyt$minv,maxv = yyt$maxv)
  xx.sel=data.frame(xx[,id.sel1,drop=FALSE],xxln[,id.sel2,drop=FALSE],check.names = FALSE)
  colnames(xx.sel)=c(colnames(xx)[id.sel1],paste0(colnames(xxln)[id.sel2],rep(".Ln",length(id.sel2))))
  list(R2cv=det.kfold$R2,R2cv.monte=det.monte$R2,R2=det.mod$R2,
       detail=det.kfold$detailn,detail.monte=det.monte$detailn,detail.mod=det.mod$detailn,
       Y=yy,X=xx.sel,Y.lgst=yyt)
}
zintfn<-function(id.train,...)
{
  id.test=which(!((1:length(x)) %in% id.train))
  env.sig.test=data.frame(sapply(1:ncol(env.sig),
                                 function(j)
                                 {
                                   mapintf(XYZ.train=data.frame(x=x[id.train],y=y[id.train],z=env.sig[id.train,j]),
                                           XY.test=data.frame(x=x[id.test],y=y[id.test]),method=method.sig[j])
                                 }))
  colnames(env.sig.test)=colnames(env.sig)
  id.use=which(rowSums(is.na(env.sig.test))==0)
  envim=rbind(env.sig,env.sig.test)
  envimln=data.frame(sapply(envim,lnx),check.names = FALSE)
  rownames(envimln)=rownames(envim)
  id.asvs=which(colnames(envimln) %in% colnames(comm))
  if(length(id.asvs)>0){envimln[,id.asvs]=log(envim[,id.asvs]+1)}
  if("pH" %in% colnames(envimln)){envimln[,"pH"]=envim[,"pH"]}
  envsigln=data.frame(envimln[1:nrow(env.sig),,drop=FALSE],check.names = FALSE)
  envsigtestln=data.frame(envimln[nrow(env.sig)+(1:nrow(env.sig.test)),,drop=FALSE],check.names = FALSE)
  envsigmm=data.frame(env.sig,envsigln,check.names = FALSE)
  envsigtestmm=data.frame(env.sig.test,envsigtestln,check.names = FALSE)
  colnames(envsigmm)<-colnames(envsigtestmm)<-c(colnames(env.sig),paste0(colnames(envsigln),".Ln"))
  if(length(id.use)>0)
  {
    regrf<-function(zii,envii,envii.test,suffixi,envname.PF,LNZ=FALSE,...)
    {
      lmx=lm(zii[id.train]~.,data = envii[id.train,match(envname.PF,colnames(envii)),drop=FALSE])
      z.lmt=predict(lmx,newdata=envii.test[id.use,match(envname.PF,colnames(envii.test)),drop=FALSE])
      plsx=try(ropls::opls(x=envii[id.train,,drop=FALSE],y=zii[id.train],fig.pdfC='none',info.txtC='none',predI=NA,orthoI=0,permI = 20),silent = TRUE)
      if(class(plsx)=="try-error"){plsx=try(ropls::opls(x=envii[id.train,,drop=FALSE],y=zii[id.train],fig.pdfC='none',info.txtC='none',predI=1,orthoI=0,permI = 20))}
      z.plst=predict(plsx,newdata=envii.test[id.use,,drop=FALSE])
      rfx=randomForest::randomForest(x=envii[id.train,,drop=FALSE],y=zii[id.train])
      z.rft=predict(rfx,newdata=envii.test[id.use,,drop=FALSE])
      if(LNZ)
      {
        outii=data.frame(test.lm=lnxrev(z.lmt,minv = zlni$minv,zerov = zlni$zerov),
                         test.pls=lnxrev(z.plst,minv = zlni$minv,zerov = zlni$zerov),
                         test.rf=lnxrev(z.rft,minv = zlni$minv,zerov = zlni$zerov))
      }else{
        outii=data.frame(test.lm=z.lmt,test.pls=z.plst,test.rf=z.rft)
      }
      colnames(outii)=paste0(colnames(outii),".",suffixi)
      outii
    }
    testZE=regrf(zii=zi,envii=env.sig,envii.test=env.sig.test,suffixi="ZE",envname.PF=envname.PF.ZE,LNZ = FALSE)
    testLnZE=regrf(zii=zlni$result,envii=env.sig,envii.test=env.sig.test,suffixi="LnZE",envname.PF=envname.PF.LnZE,LNZ = TRUE)
    testZLnE=regrf(zii=zi,envii=envsigln,envii.test=envsigtestln,suffixi="ZLnE",envname.PF=envname.PF.ZLnE,LNZ = FALSE)
    testLnZLnE=regrf(zii=zlni$result,envii=envsigln,envii.test=envsigtestln,suffixi="LnZLnE",envname.PF=envname.PF.LnZLnE,LNZ = TRUE)
    testZME=regrf(zii=zi,envii=envsigmm,envii.test=envsigtestmm,suffixi="ZME",envname.PF=envname.PF.ZME,LNZ = FALSE)
    testLnZME=regrf(zii=zlni$result,envii=envsigmm,envii.test=envsigtestmm,suffixi="LnZME",envname.PF=envname.PF.LnZME,LNZ = TRUE)
    outi=data.frame(ID=id.test[id.use],obs=zi[id.test[id.use]],testZE,testLnZE,testZLnE,testLnZLnE,testZME,testLnZME)
  }else{outi=NULL}
  outi
}
#################################
# 6.2 # calculate each zi
#################################
# 1
zname="V.select"
zi=psc[,zname]

# 2
zname="Disp.limit"
zi=psc[,zname]

# 3
zname="pH"
zi=env2[,zname]

# 4
zname="NO3"
zi=env2[,zname]

# 5
zname="U"
zi=env2[,zname]

# 6
zname="DO"
zi=env2[,zname]

# 7
zname="MSI"
zi=envm[,zname]


#################
zname %in% colnames(env.sig)

###################
zlni=lnx(zi,T)
rand=100

########
R2cvmap=R2cvf(z.obs = zi,znamei=zname)
max(R2cvmap)
save.file(data.frame(Z.name=zname,t(R2cvmap)),filename = paste0(zname,".MapIntpr.R2cv"))

########
modsl=list()
# lm: linear models
modsl$lm.ZE=modstep.cv(yy=zi,xx=env.sig,n.limit = n.train-2)
modsl$lm.ZLnE=modstep.cv(yy=zi,xx=env.sig.ln,n.limit = n.train-2)
modsl$lm.ZME=modstep2.cv(yy=zi,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2)
modsl$lm.ZME2=modstep.cv(yy=zi,xx=env.sig.mm,n.limit = n.train-2)
outlnztf<-function(outi,...)
{
  newdetail=data.frame(ID=outi$detail$ID,
                       obs=lnxrev(outi$detail$obs,minv = zlni$minv,zerov = zlni$zerov),
                       test=lnxrev(outi$detail$test,minv = zlni$minv,zerov = zlni$zerov))
  newR2cv=R2fun(newdetail$obs,newdetail$test)
  newdetail.monte=data.frame(ID=outi$detail.monte$ID,
                             obs=lnxrev(outi$detail.monte$obs,minv = zlni$minv,zerov = zlni$zerov),
                             test=lnxrev(outi$detail.monte$test,minv = zlni$minv,zerov = zlni$zerov))
  newR2cv.monte=R2fun(newdetail.monte$obs,newdetail.monte$test)
  outi$detail=newdetail
  outi$R2cv=newR2cv
  outi$detail.monte=newdetail.monte
  outi$R2cv.monte=newR2cv.monte
  outi
}
modsl$lm.LnZE=modstep.cv(yy=zlni$result,xx=env.sig,n.limit = n.train-2)
modsl$lm.LnZE=outlnztf(modsl$lm.LnZE)
modsl$lm.LnZLnE=modstep.cv(yy=zlni$result,xx=env.sig.ln,n.limit = n.train-2)
modsl$lm.LnZLnE=outlnztf(modsl$lm.LnZLnE)
modsl$lm.LnZME=modstep2.cv(yy=zlni$result,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2)
modsl$lm.LnZME=outlnztf(modsl$lm.LnZME)
modsl$lm.LnZME2=modstep.cv(yy=zlni$result,xx=env.sig.mm,n.limit = n.train-2)
modsl$lm.LnZME2=outlnztf(modsl$lm.LnZME2)

# PLS
#source("C:/Users/Daliang/Dropbox/Data_analysis_service/ieg.R.package/DevelopingFunctions/Daliang.CrossValidateLm/cv.pls.rf.r")
modsl$pls.ZE=modstep.cv(yy=zi,xx=env.sig,n.limit = n.train-2,mod = "pls") #cv.PLS(Y=zi,X=env.sig,nv=0.1,rand = 100)
modsl$pls.ZLnE=modstep.cv(yy=zi,xx=env.sig.ln,n.limit = n.train-2,mod = "pls") #cv.PLS(Y=zi,X=env.sig.ln,nv=0.1,rand = 100)
modsl$pls.ZME=modstep2.cv(yy=zi,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2,mod = "pls") #cv.PLS(Y=zi,X=env.sig.mm,nv=0.1,rand = 100)
modsl$pls.ZME2=modstep.cv(yy=zi,xx=env.sig.mm,n.limit = n.train-2,mod = "pls")

plsmdsi=modstep.cv(yy=zlni$result,xx=env.sig,n.limit = n.train-2,mod = "pls") # cv.PLS(Y=zlni$result,X=env.sig,nv=0.1,rand = 100)
modsl$pls.LnZE=outlnztf(plsmdsi)
plsmdsi2=modstep.cv(yy=zlni$result,xx=env.sig.ln,n.limit = n.train-2,mod = "pls") #cv.PLS(Y=zlni$result,X=env.sig.ln,nv=0.1,rand = 100)
modsl$pls.LnZLnE=outlnztf(plsmdsi2)
plsmdsi3=modstep2.cv(yy=zlni$result,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2,mod = "pls") # cv.PLS(Y=zlni$result,X=env.sig.mm,nv=0.1,rand = 100)
modsl$pls.LnZME=outlnztf(plsmdsi3)
plsmdsi4=modstep.cv(yy=zlni$result,xx=env.sig.mm,n.limit = n.train-2,mod = "pls")
modsl$pls.LnZME2=outlnztf(plsmdsi4)

# RF: random forest
modsl$RF.ZE=modstep.cv(yy=zi,xx=env.sig,n.limit = n.train-2,mod = "RF") #cv.RF(Y=zi,X=env.sig,nv=0.1,rand = 100)
modsl$RF.ZLnE=modstep.cv(yy=zi,xx=env.sig.ln,n.limit = n.train-2,mod = "RF") #cv.RF(Y=zi,X=env.sig.ln,nv=0.1,rand = 100)
modsl$RF.ZME=modstep2.cv(yy=zi,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2,mod = "RF") #cv.RF(Y=zi,X=env.sig.mm,nv=0.1,rand = 100)
modsl$RF.ZME2=modstep.cv(yy=zi,xx=env.sig.mm,n.limit = n.train-2,mod = "RF")

RFmdsi=modstep.cv(yy=zlni$result,xx=env.sig,n.limit = n.train-2,mod = "RF") # cv.RF(Y=zlni$result,X=env.sig,nv=0.1,rand = 100)
modsl$RF.LnZE=outlnztf(RFmdsi)
RFmdsi2=modstep.cv(yy=zlni$result,xx=env.sig.ln,n.limit = n.train-2,mod = "RF") #cv.RF(Y=zlni$result,X=env.sig.ln,nv=0.1,rand = 100)
modsl$RF.LnZLnE=outlnztf(RFmdsi2)
RFmdsi3=modstep2.cv(yy=zlni$result,xx=env.sig,xxln = env.sig.ln,n.limit = n.train-2,mod = "RF") # cv.RF(Y=zlni$result,X=env.sig.mm,nv=0.1,rand = 100)
modsl$RF.LnZME=outlnztf(RFmdsi3)
RFmdsi4=modstep.cv(yy=zlni$result,xx=env.sig.mm,n.limit = n.train-2,mod = "RF")
modsl$RF.LnZME2=outlnztf(RFmdsi4)
save(modsl,file = paste0(zname,".PredModSelection.rda"))
#l
R2cvs=sapply(modsl,function(v){v$R2cv})
max(R2cvs,na.rm=TRUE)
id.opts=order(R2cvs,decreasing = TRUE)[1:5]
#id.opt=which(R2cvs==max(R2cvs,na.rm=TRUE))
#id.opt=which(names(R2cvs)=="lm.ZME2")
#id.opt=which(names(R2cvs)=="lm.LnZLnE")
(mod.opt.names=names(modsl)[id.opts])
mod.opts=modsl[id.opts]
sapply(mod.opts,function(v){v$R2cv})
#

perm.n = permute::shuffleSet(length(x), 1000)
z.perm=sapply(1:nrow(perm.n),function(i){zi[perm.n[i,]]})
colnames(z.perm)=paste0("Rand",1:ncol(z.perm))
EPS <- sqrt(.Machine$double.eps)
R2.perm=sapply(1:ncol(z.perm),
               function(i)
               {
                 sapply(1:length(mod.opts),
                        function(j)
                        {
                          c(R2fun(yv=z.perm[mod.opts[[j]]$detail$ID,i],fv=mod.opts[[j]]$detail$test),
                            R2fun(yv=z.perm[mod.opts[[j]]$detail.monte$ID,i],fv=mod.opts[[j]]$detail.monte$test),
                            R2fun(yv=z.perm[mod.opts[[j]]$detail.mod$ID,i],fv=mod.opts[[j]]$detail.mod$test))
                        })
               },simplify = "array")
pvm=sapply(1:length(mod.opts),
           function(j)
           {
             R2obsm=matrix(c(mod.opts[[j]]$R2cv,mod.opts[[j]]$R2cv.monte,mod.opts[[j]]$R2),nrow=3,ncol = ncol(z.perm))
             (rowSums((R2.perm[,j,]-R2obsm)>=(-EPS),na.rm = TRUE)+1)/((rowSums(!is.na(R2.perm[,j,,drop=FALSE])))+1)
           })
pvm
#pv=(sums(R2.perm>=(matrix(c(mod.opt$R2cv-EPS),na.rm = TRUE)+1)/(sum(!is.na(R2.perm))+1)
#moda=lm(Yt~.,data=data.frame(Yt=mod.opt$Y,mod.opt$X,check.names = FALSE))
#Y.moda=predict(moda)
#mod.opt$R2=R2fun(mod.opt$Y,Y.moda)
opt.sum=data.frame(Z.name=zname,Mod=mod.opt.names,
                   R2cv.kfold=sapply(mod.opts,function(v){v$R2cv}),
                   R2cv.monte=sapply(mod.opts,function(v){v$R2cv.monte}),
                   R2cv.mod=sapply(mod.opts,function(v){v$R2}),
                   P.kfold=pvm[1,],P.monte=pvm[2,],P.mod=pvm[3,],
                   Predictors=sapply(mod.opts,function(v){paste(colnames(v$X),collapse = "+")}))
opt.sum
save.file(opt.sum,filename = paste0(zname,".OptMod.summary"))

###################
# 7 # draw map
###################
col2 <- colorRampPalette(c(cbcol('dblue'),cbcol('red',light.dark.amount=-0.95),cbcol('red')))
####################
# 7.1 # traditional map of zi
###################

intzi=interp(x,y,zi,xo=xo,yo=yo,duplicate = "mean",linear = TRUE)
range(zi)
range(intzi$z,na.rm = TRUE)

tiff(file=paste0(zname,".CBI.Map.tiff"),width = 18,height = 18,units = "cm",
     compression = "lzw",res=800)
par(bty="n")
filled.contour(intzi$x,intzi$y,intzi$z,color.palette = col2,xlim=c(-7000,1000),zlim=c(min(intzi$z,na.rm=TRUE),max(intzi$z,na.rm=TRUE)),asp = 1,plot.axes=NULL,axes = TRUE)
dev.off()

###################
# 7.2 # map each factors in mod.opt
###################
mto3c<-function(m,rowx,coly,na.rm=TRUE)
{
  m=as.matrix(m)
  outi=data.frame(x=rowx[as.vector(row(m))],y=coly[as.vector(col(m))],z=as.vector(m))
  if(na.rm){outi=outi[which(!is.na(outi$z)),,drop=FALSE]}
  outi
}
c3tom<-function(c3,rowx,coly)
{
  outf=matrix(NA,nrow=length(rowx),ncol=length(coly))
  idxm=match(c3[,1],rowx)
  idym=match(c3[,2],coly)
  if(sum(is.na(idxm))+sum(is.na(idym))>0){stop("c3tom not match")}
  outf[cbind(row=idxm,col=idym)]=c3[,3]
  outf
}
zm=mto3c(m=intzi$z,rowx=xo,coly=yo,na.rm = TRUE)
dim(zm)
head(zm)

if(zname %in% colnames(env.sig))
{
  method.sig.opt=method.sig[zname]
  #method.sig.opt="pred.CBI"
  pred.zi=mapintf(XYZ.train=data.frame(x=x,y=y,z=zi),XY.test=zm[,1:2,drop=FALSE],method=method.sig.opt)
  z.predim=c3tom(c3=data.frame(zm[,1:2,drop=FALSE],z=pred.zi),rowx = xo,coly = yo)
  dim(z.predim)
  fig.setting=list(x=xo,y=yo,z=z.predim,colp=col2,xlim=c(-7000,1000),zlim=c(min(z.predim,na.rm=TRUE),max(z.predim,na.rm=TRUE)),
                   maxz=max(zi),minz=min(zi),obs=data.frame(x=x,y=y,z=zi),zname=zname)
  save(fig.setting,file = paste0(zname,".0.",method.sig.opt,".Map.Setting.rda"))
  tiff(file=paste0(zname,".0.",method.sig.opt,".Map.tiff"),width = 18,height = 18,units = "cm",
       compression = "lzw",res=800)
  par(bty="n")
  filled.contour(fig.setting$x,fig.setting$y,fig.setting$z,color.palette = fig.setting$colp,
                 xlim=fig.setting$xlim,zlim=fig.setting$zlim,asp = 1,plot.axes=NULL,axes = TRUE)
  dev.off()
}


for(j in 1:length(mod.opts))
{
  mod.opt=mod.opts[[j]]
  mod.opt.name=names(mod.opts)[j]
  message("j=",j," mod.name=",mod.opt.name,". ",date())
  xnames=gsub("\\.Ln$","",colnames(mod.opt$X))
  env.predm=data.frame(sapply(1:length(xnames),function(i){mapintf(XYZ.train=data.frame(x=x,y=y,z=env.sig[,xnames[i]]),XY.test=zm[,1:2,drop=FALSE],method=method.sig[xnames[i]])}))
  dim(env.predm)
  colnames(env.predm)=colnames(mod.opt$X)
  # log transform of env.predm as needed
  lnenvf<-function(v,vname,...)
  {
    if(vname=="pH")
    {
      outf=v
    }else if(vname %in% colnames(comm)){
      outf=log(v+1)
    }else if(vname %in% colnames(env2)){
      outf=lnx(c(v,env2[,vname]))[1:length(v)]
    }else{outf=log(v)}
    outf
  }
  if(grepl("ME$",mod.opt.name) | grepl("ME2$",mod.opt.name)){id.ln=grep("\\.Ln$",colnames(mod.opt$X))}else if(grepl("LnE$",mod.opt.name)){id.ln=1:ncol(env.predm)}
  if(length(id.ln)>0){for(idi in id.ln){env.predm[,idi]=lnenvf(v=env.predm[,idi],vname=xnames[idi])}}
  # training model
  if(grepl("^lm",mod.opt.name)){modtr=lm(Y~.,data = data.frame(Y=mod.opt$Y.lgst$result,mod.opt$X,check.names = FALSE))}
  if(grepl("^pls",mod.opt.name))
  {
    modtr=try(ropls::opls(x=mod.opt$X,y=mod.opt$Y.lgst$result,fig.pdfC='none',info.txtC='none',predI=NA,orthoI=0,permI = 20),silent = TRUE)
    if(class(modtr)=="try-error"){modtr=try(ropls::opls(x=mod.opt$X,y=mod.opt$Y.lgst$result,fig.pdfC='none',info.txtC='none',predI=1,orthoI=0,permI = 20))}
  }
  if(grepl("^RF",mod.opt.name)){modtr=randomForest::randomForest(x=mod.opt$X,y=mod.opt$Y.lgst$result)}
  # prediction
  env.predm1=env.predm[which(rowSums(is.na(env.predm))==0),,drop=FALSE]
  z.predi1=lgstrev(predict(modtr,newdata=env.predm1),minv = mod.opt$Y.lgst$minv,maxv = mod.opt$Y.lgst$maxv)
  if(grepl("LnZ",mod.opt.name)){z.predi1=lnxrev(z.predi1,minv = zlni$minv,zerov = zlni$zerov)}
  z.predi=rep(NA,nrow(env.predm))
  z.predi[which(rowSums(is.na(env.predm))==0)]=z.predi1
  length(z.predi)
  range(z.predi,na.rm=TRUE)
  range(zi)
  z.predim=c3tom(c3=data.frame(zm[,1:2,drop=FALSE],z=z.predi),rowx = xo,coly = yo)
  dim(z.predim)
  fig.setting=list(x=xo,y=yo,z=z.predim,colp=col2,xlim=c(-7000,1000),zlim=c(min(z.predim,na.rm=TRUE),max(z.predim,na.rm=TRUE)),
                   maxz=max(zi),minz=min(zi),obs=data.frame(x=x,y=y,z=zi),zname=zname)
  save(fig.setting,file = paste0(zname,".",j,".",mod.opt.name,".Map.Setting.rda"))
  #filled.contour(fig.setting$x,fig.setting$y,fig.z,color.palette = fig.setting$colp,
  #               xlim=fig.setting$xlim,zlim=c(min(fig.z,na.rm=TRUE),max(fig.z,na.rm = TRUE)),asp = 1,plot.axes=NULL,axes = TRUE)
  
  #filled.contour(xo,yo,z.predim,color.palette = col2,xlim=c(-7000,1000),asp = 1,plot.axes=NULL,axes = TRUE)
  tiff(file=paste0(zname,".",j,".",mod.opt.name,".ModPred.Map.tiff"),width = 18,height = 18,units = "cm",
       compression = "lzw",res=800)
  par(bty="n")
  filled.contour(fig.setting$x,fig.setting$y,fig.setting$z,color.palette = fig.setting$colp,
                 xlim=fig.setting$xlim,zlim=fig.setting$zlim,asp = 1,plot.axes=NULL,axes = TRUE)
  dev.off()
}

#end



