# version 2020.8.20 for GitHub backup

list.of.packages <- c("ape", "seqinr","vegan","bigmemory","MASS","iCAMP","NST")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

iwd<-function(wd=getwd())
{
  wd=sub("file:///","",wd)
  wd=gsub("\\\\","/",wd)
  if(length(list.dirs(wd))==0)
  {
    wd=substr(wd,1,regexpr("\\/[^\\/]*$",wd)-1)
  }
  setwd(wd)
  message("set wd as ",getwd(),"\n",date())
  getwd()
}

open.file<-function(filename,folder=NA,header=TRUE,row.names=1,load.obj=1,stringsAsFactors=FALSE,check.names=FALSE,...)
{
  # version 20160730 by Daliang Ning (ningdaliang@gmail.com)
  # version 20171104 by Daliang Ning (ningdaliang@gmail.com)
  # version 20200612 by Daliang Ning (ningdaliang@gmail.com), add tsv
  if(is.na(folder)){folder=getwd()}
  suffix=regexpr("\\.[^\\.]*$", filename)
  filetype=substr(filename,suffix,nchar(filename))
  if(filetype %in% c(".nwk",".tre",".tree"))
  {
    #library(ape)
    file=ape::read.tree(file=paste(folder,"/",filename,sep = ""),...)
  }else if(filetype==".fna"|filetype==".fasta"){
    #library(seqinr)
    file=seqinr::read.fasta(file = paste(folder,"/",filename,sep = ""),seqtype = "DNA",as.string = FALSE,forceDNAtolower = FALSE,set.attributes = TRUE,legacy.mode = TRUE,seqonly = FALSE,strip.desc = FALSE,bfa = FALSE)
  }else if(filetype==".RData"){
    if(load.obj=="all")
    {
      old.ls=ls(envir = .GlobalEnv)
      sp=load(paste(folder,"/",filename,sep = ""),envir = .GlobalEnv)
      new.ls=sp
      message("All objects are loaded, including ",sp)
      message("Replaced objects include ",intersect(old.ls,new.ls))
      file=sp
    }else{
      sp=load(paste(folder,"/",filename,sep = ""))
      message("Only one object is loaded: ",sp[load.obj])
      file=get(sp[load.obj])
    }
  }else if(filetype==".rda")
  {
    sp=load(paste(folder,"/",filename,sep = ""))
    file=get(sp)
  }else{
    if(filetype %in% c(".txt",".tsv")){delim="\t"}else if(filetype==".csv"){delim=","}else if(filetype==".tabular"){delim=""}else{message("The file type of ",filename," is not included in this function.")}
    file=read.table(file = paste(folder,"/",filename,sep = ""),header=header,sep = delim,
                    row.names = row.names,stringsAsFactors=stringsAsFactors,check.names = check.names,
                    comment.char = "",...)
  }
  file
}

lazyopen<-function(filepath,load.obj = "all",...)
{
  if(is.null(filepath)){file=NULL}else{
    filepath=gsub("\\\\","/",filepath)
    if(regexpr("/",filepath)>0)
    {
      filepath=sub("file:///","",filepath)
      locate=as.vector(gregexpr("/",filepath)[[1]]);locate=locate[length(locate)]
      file=open.file(filename=substring(filepath,locate+1),folder = substr(filepath,1,locate-1),load.obj = load.obj,...)
    }else{
      file=open.file(filename=filepath,folder = getwd(),load.obj = load.obj,...)
    }
  }
  file
}

save.file<-function(obj,prefix=numeric(0),filename=NA,folder=NA,
                    type=c("csv","txt","RData","nwk","fasta"),
                    rowname.title="ID", replace=TRUE,row.names=TRUE,...)
{
  # version 20171109
  if(is.null(obj)){message(deparse(substitute(obj))," is NULL, nothing to save.");return(invisible())}
  if(is.na(folder)){folder=getwd()}
  if(is.na(filename)){filename=as.character(substitute(obj))}
  if(length(prefix)!=0){prefix=paste(prefix,".",sep = "")}
  tcode=format(Sys.time(),format = "%y%m%H%M%S")
  if(type[1]=="csv")
  {
    if(row.names)
    {
      obj=data.frame(obj,stringsAsFactors = FALSE)
      if(is.null(rownames(obj))) rownames(obj)=1:nrow(obj)
      obj=data.frame(rownames(obj),obj)
      colnames(obj)[1]=rowname.title
    }
    if(file.exists(paste(folder,"/",prefix,filename,".csv",sep = ""))&(!replace))
    {
      write.csv(obj,file = paste(folder,"/",prefix,filename,".",tcode,".csv",sep = ""),
                row.names=FALSE,...)
    }else{
      write.csv(obj,file = paste(folder,"/",prefix,filename,".csv",sep = ""),
                row.names=FALSE,...)
    }
  }else if(type[1]=="RData"){
    if(file.exists(paste(folder,"/",prefix,filename,".RData",sep = ""))&(!replace))
    {
      save(obj,file = paste(folder,"/",prefix,filename,".",tcode,".RData",sep = ""))
    }else{
      save(obj,file = paste(folder,"/",prefix,filename,".RData",sep = ""))
    }
  }else if(type[1]=="txt"){
    if(row.names)
    {
      obj=data.frame(obj,stringsAsFactors = FALSE)
      if(is.null(rownames(obj))) rownames(obj)=1:nrow(obj)
      obj=data.frame(rownames(obj),obj)
      colnames(obj)[1]=rowname.title
    }
    if(file.exists(paste(folder,"/",prefix,filename,".txt",sep = ""))&(!replace))
    {
      write.table(obj,file=paste(folder,"/",prefix,filename,".",tcode,".txt",sep = ""),
                  quote=FALSE,sep="\t",row.names=FALSE,...)
    }else{
      write.table(obj,file=paste(folder,"/",prefix,filename,".txt",sep = ""),
                  quote=FALSE,sep="\t",row.names=FALSE,...)
    }
  }else if(type[1]=="nwk"){
    #library(ape)
    if(file.exists(paste(folder,"/",prefix,filename,".nwk",sep = ""))&(!replace))
    {
      ape::write.tree(obj,file= paste(folder,"/",prefix,filename,".",tcode,".nwk",sep = ""))
    }else{
      ape::write.tree(obj,file= paste(folder,"/",prefix,filename,".nwk",sep = ""))
    }
  }else if(type[1]=="fasta"){
    #library(seqinr)
    if(file.exists(paste(folder,"/",prefix,filename,".fasta",sep = ""))&(!replace))
    {
      seqinr::write.fasta(obj,names = names(obj),file.out = paste(folder,"/",prefix,filename,".",tcode,".fasta",sep = ""))
    }else{
      seqinr::write.fasta(obj,names = names(obj),file.out = paste(folder,"/",prefix,filename,".fasta",sep = ""))
    }
  }else{
    warning("The file ",filename," was not saved, because type is not included in this function.")
  }
  invisible()
}

ccc<-function(x,y,out.vector=FALSE)
{
  d <- y - x
  m1 <- mean(y)
  m2 <- mean(x)
  v1 <- var(y)
  v2 <- var(x)
  n <- length(d)
  e2 <- sum(d^2)/n
  if (length(y) != length(x)) 
    stop("x and y should have the same length!")
  mu_d <- m1 - m2
  d2 <- mu_d^2
  s12 <- v1 * (n - 1)/n
  s22 <- v2 * (n - 1)/n
  U <- mu_d/sqrt(sqrt(s12 * s22))
  V <- sqrt(s12/s22)
  Ca <- 2/(V + 1/V + U^2)
  rc <- 1 - e2/(d2 + s12 + s22)
  r <- (rc/Ca)
  out=list(accuracy=Ca,precision=r,ccc=rc)
  if(out.vector)
  {
    out=unlist(out)
  }
  out
}

col3.dist<-function(m.3col,to.dist=TRUE)
{
  # convert 3-col matrix to distance object
  m.3col=as.matrix(m.3col)
  
  name.factor=as.factor(c(m.3col[,1],m.3col[,2]))
  name.lev=levels(name.factor)
  num=length(name.lev)
  
  if(sum(m.3col[,1]==m.3col[,2])>0)
  {
    if(to.dist){warning("In some samples, distance within sample is not zero. it is better to return Matrix instead of dist.")}
    idd=which(m.3col[,1]==m.3col[,2])
    m.diag=m.3col[idd,,drop=FALSE]
    m.3col=m.3col[-idd,,drop=FALSE]
  }else{idd=NA}
  
  if(nrow(m.3col)>(0.5*num*(num-1)))
  {
    warning("The row number is more than all possible pairwise comparisons, some duplicate(s) were ignored.")
  }else if(nrow(m.3col)<(0.5*num*(num-1)))
  {
    warning("The row number is less than all possible pairwise comparisons, NA was returned.")
  }
  res=matrix(NA,nrow=length(name.lev),ncol=length(name.lev))
  rownames(res)<-colnames(res)<-name.lev
  res.name=paste(name.lev[row(res)],name.lev[col(res)],sep = ".")
  m.name=cbind(c(paste(m.3col[,1],m.3col[,2],sep = "."),paste(m.3col[,2],m.3col[,1],sep = ".")),c(m.3col[,3],m.3col[,3]))
  res[]=as.numeric(m.name[match(res.name,m.name[,1]),2])
  if(to.dist)
  {
    res=as.dist(res)
  }else{
    diag(res)=0
    if(!is.na(idd)){diag(res)[match(m.diag[,1],name.lev)]=as.numeric(m.diag[,3])}
  }
  res
}

#########################################
# beta version of mcMRM
# author: Daliang Ning (ningdaliang@ou.edu)
#####
# I # Calculate R2 R2adj
#####
r2ma<-function(mod,y=NA,p=NA)
{
  if(class(mod)[1]=="try-error"){out=c(R2=NA,R2adj=NA)}else{
    if("lm" %in% class(mod)){y=mod$model[,1]}
    if(grepl("lmerMod",class(mod))){y=mod@frame[,1]}
    yh=fitted(mod)
    n=length(y)
    if("lm" %in% class(mod)){p=ncol(mod$model)-1}
    if(grepl("lmerMod",class(mod))){p=ncol(mod@frame)-1}
    r2=1-sum((y-yh)^2)/sum((y-mean(y))^2)
    if(is.na(p)){out=c(R2=NA,R2adj=NA)}else if((n-p-1)<=0){out=c(R2=r2,R2adj=NA)}else{
      out=c(R2=r2,R2adj=1-((1-r2)*(n-1)/(n-p-1)))
    }
  }
  out
}

######
# II # Forward selection function
######
fw.mod<-function(y,datin,method=c("lm","glm","lmm"),glm.family=quasibinomial,
                 lmm.formula=NULL,lmm.datr=NULL,silent=FALSE,dat.perm=NULL,mod1.perm=TRUE)
{
  EPS=sqrt(.Machine$double.eps)
  datin=data.frame(datin,stringsAsFactors = FALSE)
  pt=ncol(datin)
  df1=lapply(1:ncol(datin),function(i){datin[,i,drop=FALSE]})
  
  comp1=t(sapply(1:length(df1),
                 function(i)
                 {
                   if(!silent) message("--mod1 i=",i," in ",length(df1),".--",date())
                   if(method=="lm"){lmi=try(lm(y~.,data = df1[[i]]))}else if(method=="glm"){
                     lmi=try(glm(y~.,data = df1[[i]],family = glm.family))
                   }else if(method=="lmm"){
                     formu=formula(sub("\\.",paste0(colnames(df1[[i]]),collapse = "+"),deparse(lmm.formula)))
                     lmi=try(lmer(formu,data = data.frame(df1[[i]],y=y,lmm.datr,stringsAsFactors = FALSE)))
                   }
                   if(class(lmi)[1]=="try-error")
                   {
                     out=c(coef=NA,t.value=NA,p.ttest=NA,p.coef.perm=NA,R2=NA,R2adj=NA,p.mod.perm=NA)
                   }else{
                     lmsi=summary(lmi)
                     if(method=="lmm"){coefi=fixef(lmi)[[2]]}else{coefi=lmi$coefficients[[2]]}
                     ti.obs=lmsi$coefficients[,"t value"][[2]]
                     pi.t=lmsi$coefficients[,"Pr(>|t|)"][[2]]
                     r2i=r2ma(lmi)
                     if((!is.null(dat.perm))&(mod1.perm))
                     {
                       r2ti.obs=c(r2i[[1]],ti.obs)
                       r2ti.rd=sapply(1:ncol(dat.perm),
                                      function(j)
                                      {
                                        if((method=="lmm")&(!silent)&(j %in% seq(from=1,to=ncol(dat.perm),by=100))) message("--mod1 i=",i," j=",j," in ",ncol(dat.perm),".--",date())
                                        yrj=dat.perm[,j] #y[permid[[j]]]
                                        if(method=="lm"){lmrj=try(lm(yrj~.,data=df1[[i]]))}else if(method=="glm"){
                                          lmrj=try(glm(yrj~.,data=df1[[i]],family = glm.family))
                                        }else if(method=="lmm"){
                                          lmrj=try(lmer(formu,data = data.frame(df1[[i]],y=yrj,lmm.datr,stringsAsFactors = FALSE)))
                                        }
                                        if(class(lmrj)=="try-error"){outt=rep(NA,2)}else{
                                          outt=c(r2ma(lmrj)[[1]],summary(lmrj)$coefficients[,"t value"][[2]])
                                        }
                                        outt
                                      })
                       pi.pm=(rowSums(r2ti.rd>=(r2ti.obs-EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2ti.rd))+1)
                       pi.pm[which(r2ti.obs<0)]=((rowSums(r2ti.rd<=(r2ti.obs+EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2ti.rd))+1))[which(r2ti.obs<0)]
                     }else{pi.pm=rep(NA,2)}
                     out=c(coef=coefi,t.value=ti.obs,p.ttest=pi.t,p.coef.perm=pi.pm[2],r2i,p.mod.perm=pi.pm[1])
                   }
                   out
                 }))
  rownames(comp1)=colnames(datin)
  comp=comp1[,5:6,drop=FALSE]
  id1=which(comp[,"R2adj"]==max(comp[,"R2adj"]))
  comps=list()
  
  comps[[1]]=cbind(t(sapply(1:length(id1),function(i){out=rep(0,pt);out[id1[i]]=1;out})),
                   comp[id1,,drop=FALSE])
  ids=lapply(1:length(id1),function(i){id1[[i]]})
  
  
  m=2
  while(length(ids[[1]])<ncol(datin) & (length(ids[[1]])<(nrow(datin)-2)) & m<1000)
  {
    if(!silent) message("--Now m=",m," R2adj=",comps[[m-1]][,ncol(comps[[m-1]])][[1]],".--",date())
    idcm=Reduce(rbind,lapply(1:length(ids),function(i){id.rem=(1:pt)[-ids[[i]]];t(sapply(id.rem,function(j){sort(c(ids[[i]],j))}))}))
    idcm=idcm[!duplicated(idcm),,drop=FALSE]
    dfm=lapply(1:nrow(idcm),function(i){datin[,idcm[i,],drop=FALSE]})
    compm=t(sapply(1:length(dfm),
                   function(i)
                   {
                     if(method=="lm"){lmi=try(lm(y~.,data = dfm[[i]]))}else if(method=="glm"){
                       lmi=try(glm(y~.,data = dfm[[i]],family = glm.family))
                     }else if(method=="lmm"){
                       formui=formula(sub("\\.",paste0(colnames(dfm[[i]]),collapse = "+"),deparse(lmm.formula)))
                       lmi=try(lmer(formui,data = data.frame(dfm[[i]],y=y,lmm.datr,stringsAsFactors = FALSE)))
                     }
                     r2ma(lmi)
                   }))
    idm=which(compm[,"R2adj"]==max(compm[,"R2adj"],na.rm=TRUE))
    if(length(idm)>0)
    {
      comps[[m]]=cbind(t(sapply(1:length(idm),function(i){out=rep(0,pt);out[idcm[idm[i],]]=1;out})),
                       compm[idm,,drop=FALSE])
      ids=lapply(1:length(idm),function(i){idcm[idm[[i]],]})
      m=m+1
    }else{m=1000}
  }
  compsm=Reduce(rbind,comps)
  idm=min(which(compsm[,"R2adj"]==max(compsm[,"R2adj"])))
  compso=compsm[1:idm,,drop=FALSE]
  fns=rowSums(compso[,1:(ncol(compso)-2),drop=FALSE])
  
  # remove some rows not leading to the final optimum combination
  if(max(fns)>1)
  {
    idi=max(which(fns==(max(fns)-1)))
    id.rm=integer(0)
    for(i in idi:1)
    {
      if(sum(fns==fns[i])>1)
      {
        idh=which(fns==(fns[i]+1))
        idh=idh[which(!(idh %in% id.rm))]
        jud=prod(sapply(1:length(idh),function(j){sum((compso[i,1:(ncol(compso)-2)]-compso[idh[j],1:(ncol(compso)-2)])>0)}))
        if(jud>0) id.rm=c(id.rm,i)
      }
    }
    if(length(id.rm)>0){compso=compso[-id.rm,,drop=FALSE];fns=fns[-id.rm]}
  }
  
  pre.id<-new.factor<-list()
  pre.id[[1]]=rep(0,sum(fns==1))
  new.factor[[1]]=sapply(which(fns==1),function(i){colnames(datin)[which(compso[i,1:(ncol(compso)-2)]>0)]})
  k=2
  if(max(fns)>1)
  {
    idi=min(which(fns==2))
    for(i in idi:nrow(compso))
    {
      idl=which(fns==(fns[i]-1))
      jud=sapply(idl,function(j){sum((compso[j,1:(ncol(compso)-2)]-compso[i,1:(ncol(compso)-2)])>0)})
      pre.id[[k]]=min(idl[jud==0])
      new.factor[[k]]=colnames(datin)[which((compso[pre.id[[k]],1:(ncol(compso)-2)]-compso[i,1:(ncol(compso)-2)])<0)]
      k=k+1
    }
  }
  compsim=data.frame(Previous.Step=unlist(pre.id),
                     Add.Factor=unlist(new.factor),
                     R2.cum=compso[,ncol(compso)-1],
                     R2adj.cum=compso[,ncol(compso)],
                     stringsAsFactors = FALSE)
  
  id.sel=which(compso[nrow(compso),1:(ncol(compso)-2)]>0)
  dat.sel=datin[,id.sel[order(match(colnames(datin)[id.sel],compsim$Add.Factor))],drop=FALSE]
  if(method=="lm"){mod=lm(y~.,data = dat.sel)}else if(method=="glm"){
    mod=glm(y~.,data = dat.sel,family = glm.family)
  }else if(method=="lmm"){
    formu.sel=formula(sub("\\.",paste0(colnames(dat.sel),collapse = "+"),deparse(lmm.formula)))
    mod=lmer(formu.sel,data = data.frame(dat.sel,y=y,lmm.datr,stringsAsFactors = FALSE))
  }
  mods=summary(mod)
  if(method=="lmm"){coef.obs=fixef(mod)}else{coef.obs=mod$coefficients}
  t.obs=mods$coefficients[,"t value"]
  anv=anova(mod)
  if(method=="lm")
  {
    SumSq=anv$`Sum Sq`[match(compsim$Add.Factor,rownames(anv))]
    TSS=sum(anv$`Sum Sq`)
  }else if(method=="glm"){
    SumSq=anv$Deviance[match(compsim$Add.Factor,rownames(anv))]
    TSS=anv$`Resid. Dev`[1]
  }else if(method=="lmm"){
    SumSq=anv$`Sum Sq`[match(compsim$Add.Factor,rownames(anv))]
    TSS=sum((mod@frame[,1]-mean(mod@frame[,1]))^2)
  }
  coef1=coef.obs[match(compsim$Add.Factor,names(coef.obs))];names(coef1)=c()
  t.value1=t.obs[match(compsim$Add.Factor,names(t.obs))];names(t.value1)=c()
  p.para1=mods$coefficients[,"Pr(>|t|)"][match(compsim$Add.Factor,rownames(mods$coefficients))];names(p.para1)=c()
  out1=data.frame(compsim,coef=coef1,
                  VariancePortion=SumSq/TSS,t.value=t.value1,
                  P.parametric=p.para1,
                  stringsAsFactors = FALSE)
  R2.out=r2ma(mod,y)
  
  if(!is.null(dat.perm))
  {
    r2c.obs=c(model.R2=R2.out[[1]],t.obs)
    r2c.rand=sapply(1:ncol(dat.perm),
                    function(i)
                    {
                      if((!silent)&(i %in% seq(from=1,to=ncol(dat.perm),by=100))){message("---perm i=",i,".---",date())}
                      yri=dat.perm[,i] #y[permid[[i]]]
                      if(method=="lm"){modri=try(lm(yri~.,data = dat.sel))}else if(method=="glm"){
                        modri=try(glm(yri~.,data = dat.sel,family = glm.family))
                      }else if(method=="lmm"){
                        datri=mod@frame
                        datri[,1]=dat.perm[,i]
                        modri=try(lmer(formu.sel,data = datri))
                      }
                      if(class(modri)[1]=="try-error"){out=rep(NA,length(r2c.obs))}else{
                        modsi=summary(modri)
                        tri=modsi$coefficients[,"t value"][match(names(t.obs),rownames(modsi$coefficients))];names(tri)=c()
                        out=c(r2ma(modri,yri)[[1]],tri)
                      }
                      out
                    })
    p.perm=(rowSums(r2c.rand>=(r2c.obs-EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2c.rand))+1)
    p.neg=(rowSums(r2c.rand<=(r2c.obs+EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2c.rand))+1)
    p.perm[which(r2c.obs<0)]=p.neg[which(r2c.obs<0)]
    p.perm1=p.perm[match(compsim$Add.Factor,names(r2c.obs))];names(p.perm1)=c()
    out1=data.frame(out1,P.perm=p.perm1,stringsAsFactors = FALSE)
    P.perm.out=p.perm[[1]]
  }else{P.perm.out=NA}
  
  rownames(out1)=1:nrow(out1)
  list(summary=out1,mod=mod,R2=R2.out[[1]],R2adj=R2.out[[2]],AIC=AIC(mod),P.perm=P.perm.out,mod1=comp1)
}


#######
# III # main function of mcMRM
#######
mcMRM<-function(y.3col,data.x,y.ids=NULL,grp.rand=NULL,grp.const=NULL,scale.yn,
                method=c("lm","glm","lmm"),forward=FALSE,mod1.perm=FALSE,
                glm.family=quasibinomial,lmm.formula=NULL,lmm.dat.rand=NULL,
                rand=1000,silent=FALSE)
{
  # multi-dimension constrained MRM based on lm, glm, and/or lmm
  # grp.rand defines the groups to permutate
  # grp.const defines the groups which should be maintained the same as observed.
  #####
  # 1 # match IDs
  #####
  if(is.null(y.ids)){y.use=y.3col}else{y.use=y.3col[y.ids,,drop=FALSE]}
  n1a=nrow(y.use);n2a=nrow(data.x)
  for(i in 1:nrow(y.3col))
  {
    y.3col[i,1:2]=sort(c(y.3col[i,1],y.3col[i,2]))
  }
  for(i in 1:nrow(y.use))
  {
    y.use[i,1:2]=sort(c(y.use[i,1],y.use[i,2]))
  }
  for(i in 1:nrow(data.x))
  {
    data.x[i,1:2]=sort(c(data.x[i,1],data.x[i,2]))
  }
  samp2y=paste0(y.use[,1],"__",y.use[,2])
  samp2x=paste0(data.x[,1],"__",data.x[,2])
  samp2=intersect(samp2x,samp2y)
  
  if(!is.null(lmm.dat.rand))
  {
    for(i in 1:nrow(lmm.dat.rand))
    {
      lmm.dat.rand[i,1:2]=sort(c(lmm.dat.rand[i,1],lmm.dat.rand[i,2]))
    }
    samp2r=paste0(lmm.dat.rand[,1],"__",lmm.dat.rand[,2])
    samp2=intersect(samp2,samp2r)
    n5a=nrow(lmm.dat.rand)
  }
  y.use=y.use[match(samp2,samp2y),,drop=FALSE]
  data.x=data.x[match(samp2,samp2x),,drop=FALSE]
  if(!is.null(lmm.dat.rand)) lmm.dat.rand=lmm.dat.rand[match(samp2,samp2r),,drop=FALSE]
  samps=unique(as.vector(as.matrix(y.use[,1:2])))
  if(!is.null(grp.rand)){samps=intersect(samps,rownames(grp.rand));n3a=nrow(grp.rand)}
  if(!is.null(grp.const)){samps=intersect(samps,rownames(grp.const));n4a=nrow(grp.const)}
  y.in=y.3col[which((y.3col[,1] %in% samps)&(y.3col[,2] %in% samps)),,drop=FALSE]
  y.use=y.use[which((y.use[,1] %in% samps)&(y.use[,2] %in% samps)),,drop=FALSE]
  if(n1a>nrow(y.use)) message("y.use lost ",n1a-nrow(y.use)," rows")
  data.x=data.x[which((data.x[,1] %in% samps)&(data.x[,2] %in% samps)),,drop=FALSE]
  if(n2a>nrow(data.x)) message("data.x lost ",n2a-nrow(data.x)," rows")
  if(!is.null(grp.rand))
  {
    grp.rand=grp.rand[match(samps,rownames(grp.rand)),,drop=FALSE]
    if(n3a>nrow(grp.rand)) message("grp.rand lost ",n3a-nrow(grp.rand)," rows")
  }
  if(!is.null(grp.const))
  {
    grp.const=grp.const[match(samps,rownames(grp.const)),,drop=FALSE]
    if(n4a>nrow(grp.const)) message("grp.const lost ",n4a-nrow(grp.rand))
  }
  if(!is.null(lmm.dat.rand))
  {
    lmm.dat.rand=lmm.dat.rand[which((lmm.dat.rand[,1] %in% samps)&(lmm.dat.rand[,2] %in% samps)),,drop=FALSE]
    if(n5a>nrow(lmm.dat.rand)) message("lmm.dat.rand lost ",n5a-nrow(lmm.dat.rand)," rows")
  }
  
  samp2ya=paste0(y.in[,1],"__",y.in[,2])
  
  #####
  # 2 # data.obs standardize
  #####
  if(scale.yn)
  {
    dats=data.frame(scale(data.x[,3:ncol(data.x)]),stringsAsFactors = FALSE)
    if(sum(is.na(dats))>0)
    {
      dats=dats[,which(colSums(is.na(dats))==0),drop=FALSE]
      warning("some input data has no variation. removed ",
              paste(colnames(dats)[which(colSums(is.na(dats))>0)],collapse = ", "))
    }
    data.x=cbind(data.x[,1:2,drop=FALSE],as.matrix(dats))
  }
  #####
  # 3 # data permutation
  #####
  if(is.null(grp.rand)){grp.rand=data.frame(sample=samps,stringsAsFactors = FALSE);rownames(grp.rand)=samps}
  if(is.null(rand) | is.na(rand)){pdat=NULL}else{
    permf<-function(randn,...)
    {
      t.lev<-list()
      if(is.null(grp.const))
      {
        for(i in 1:ncol(grp.rand))
        {
          message("permutating group type i=",i," name=",colnames(grp.rand)[i],". ",date())
          t.lev[[i]]=unique(grp.rand[,i])
          permi=permute::shuffleSet(length(t.lev[[i]]),randn)
          if(nrow(permi)<randn) permi=rbind(permi,1:ncol(permi))
          pti=sapply(1:nrow(permi),function(j){t.lev[[i]][permi[j,match(grp.rand[,i],t.lev[[i]])]]})
          if(i==1){pts=pti;ctsi=grp.rand[,i]}else{
            ctsi=sapply(1:nrow(grp.rand),function(j){paste(grp.rand[j,1:i],collapse = "__")})
            ptsi=lapply(1:ncol(pti),
                        function(j)
                        {
                          outj=sapply(1:ncol(pts),function(k){paste0(pts[,k],"__",pti[,j])})
                          cj=matrix((!(outj %in% ctsi)),nrow = nrow(outj),ncol=ncol(outj))
                          outj[,colSums(cj)==0,drop=FALSE]
                        })
            pts=Reduce(cbind,ptsi)
          }
        }
      }else{
        
        ctc=sapply(1:nrow(grp.const),function(i){paste(grp.const[i,],collapse = "__")})
        c.lev=unique(ctc)
        
        for(i in 1:ncol(grp.rand))
        {
          message("permutating group type i=",i," name=",colnames(grp.rand)[i],". ",date())
          ti.lev=unique(grp.rand[,i])
          
          ti.ct=sapply(1:length(ti.lev),
                       function(j)
                       {
                         paste(sort(unique(ctc[grp.rand[,i]==ti.lev[j]])),collapse = "__")
                       })
          tict.lev=unique(ti.ct)
          idj=lapply(1:length(tict.lev),function(j){which(ti.ct==tict.lev[j])})
          idjr=lapply(1:length(tict.lev),
                      function(j)
                      {
                        if(length(idj[[j]])==1){matrix(idj[[j]],1,1)}else{
                          permij=permute::shuffleSet(length(idj[[j]]),randn)
                          if(nrow(permij)<randn) permij=rbind(permij,1:ncol(permij))
                          matrix(idj[[j]][permij],nr=nrow(permij),nc=ncol(permij))
                        }
                      })
          prdn=prod(sapply(idjr,nrow))
          if(prdn>randn)
          {
            tpid=t(sapply(1:randn,
                          function(j)
                          {
                            outj=rep(0,length(ti.lev))
                            for(k in 1:length(idjr))
                            {
                              outj[idj[[k]]]=idjr[[k]][sample(nrow(idjr[[k]]),1),]
                            }
                            outj
                          }))
          }else{
            for(j in 1:length(tict.lev))
            {
              if(j==1){tpid=t(sapply(1:nrow(idjr[[1]]),function(k){out=rep(0,length(ti.lev));out[idj[[1]]]=idjr[[1]][k,];out}))}else{
                tpidl=lapply(1:nrow(idjr[[j]]),
                             function(k)
                             {
                               outk=tpid
                               outk[,idj[[j]]]=matrix(idjr[[j]][k,],nr=nrow(tpid),nc=ncol(idjr[[j]]),byrow = TRUE)
                               outk
                             })
                tpid=Reduce(rbind,tpidl)
              }
            }
          }
          permi=tpid
          pti=sapply(1:nrow(permi),function(j){ti.lev[permi[j,match(grp.rand[,i],ti.lev)]]})
          if(i==1){pts=sapply(1:ncol(pti),function(j){paste0(ctc,"__",pti[,j])});ctsi=paste0(ctc,"__",grp.rand[,i])}else{
            ctsi=sapply(1:nrow(grp.rand),function(j){paste(c(ctc[j],as.vector(as.matrix(grp.rand[j,1:i]))),collapse = "__")})
            ptsi=lapply(1:ncol(pti),
                        function(j)
                        {
                          outj=sapply(1:ncol(pts),function(k){paste0(pts[,k],"__",pti[,j])})
                          cj=matrix((!(outj %in% ctsi)),nrow = nrow(outj),ncol=ncol(outj))
                          outj[,colSums(cj)==0,drop=FALSE]
                        })
            pts=Reduce(cbind,ptsi)
          }
        }
      }
      
      ncol(pts)
      
      ct.lev=unique(ctsi)
      if(length(ct.lev)==length(ctsi))
      {
        pid=matrix(match(pts,ctsi),nr=nrow(pts),nc=ncol(pts))
      }else{
        message("the combination of grp info is not unique ID. may use random match when permutation. ",date())
        pid=sapply(1:ncol(pts),
                   function(j)
                   {
                     ptij=pts[,j]
                     outj=rep(NA,nrow(pts))
                     for(i in 1:length(ct.lev))
                     {
                       ctii=ct.lev[i]
                       idai=which(ctsi==ctii)
                       idij=which(ptij==ctii)
                       idaij=rep(idai,ceiling(length(idij)/length(idai)))
                       outj[idij]=sample(idaij,length(idij))
                     }
                     outj
                   })
      }
      
      id.rm=which(colSums(pid!=(1:nrow(pid)))==0)
      if(length(id.rm)>0){pid=pid[,-id.rm,drop=FALSE];endc=TRUE}else{endc=FALSE}
      
      pdsn=0
      pds=matrix(nr=nrow(y.use),nc=0)
      pid.res<-pid
      t=1
      while(pdsn<randn & ncol(pid.res)>0)
      {
        message("matching perm data t=",t,". ",date())
        if(ncol(pid.res)>randn)
        {
          idp.use=sample(ncol(pid.res),randn)
          pid=pid.res[,idp.use,drop=FALSE]
          pid.res=pid.res[,-idp.use,drop=FALSE]
        }else{
          pid=pid.res
          pid.res=matrix(nr=nrow(pid.res),nc=0)
        }
        
        pdst=sapply(1:ncol(pid),
                    function(i)
                    {
                      pidi=pid[,i]
                      pn2i=matrix(rownames(grp.rand)[pidi[match(as.matrix(y.use[,1:2]),rownames(grp.rand))]],nr=nrow(y.use),nc=2)
                      pn2c=sapply(1:nrow(pn2i),function(j){paste(sort(pn2i[j,]),collapse = "__")})
                      match(pn2c,samp2ya)
                    })
        pdst=pdst[,colSums(is.na(pdst))==0,drop=FALSE]
        pds=cbind(pds,pdst)
        pdsn=ncol(pds)
        t=t+1
      }
      pdat=matrix(y.in[pds,3],nr=nrow(pds),nc=ncol(pds))
      list(pdat=pdat,endc=endc)
    }
    
    
    ptf=permf(randn=rand)
    pdat=ptf$pdat
    endc=ptf$endc
    
    if(ncol(pdat)>=rand)
    {
      pdat=pdat[,sample(ncol(pdat),rand),drop=FALSE]
    }else{
      t=2
      while(ncol(pdat)<rand & (!endc) & t<5)
      {
        message("Repeating perm t=",t,". ",date())
        ptf=permf(randn=rand*t)
        pdat=ptf$pdat
        endc=ptf$endc
        t=t+1
      }
      if(endc){message("All possible permuations are used.")}
      if(ncol(pdat)<rand){warning("Actual permutation time is ",ncol(pdat))}
      if(ncol(pdat)<=10){stop("possible permutation time is no more than 10.")}
    }
  }
  #####
  # 4 # if NOT forward selection
  #####
  res<-list()
  y=y.use[,3]
  datin=data.frame(data.x[,3:ncol(data.x),drop=FALSE],stringsAsFactors = FALSE)
  EPS = sqrt(.Machine$double.eps)
  if("lmm" %in% method)
  {
    library(lme4)
    library(lmerTest)
    if(is.null(lmm.dat.rand))
    {
      if(!is.null(grp.rand))
      {
        idr1=match(y.use[,1],rownames(grp.rand))
        idr2=match(y.use[,2],rownames(grp.rand))
        dtr1=sapply(1:ncol(grp.rand),
                    function(i)
                    {
                      sapply(1:nrow(y.use),function(j){paste(sort(c(grp.rand[idr1[j],i],grp.rand[idr2[j],i])),collapse = "__")})
                    })
        colnames(dtr1)=colnames(grp.rand)
      }else{dtr1=matrix(nr=nrow(y.use),nc=0)}
      if(!is.null(grp.const))
      {
        idr1=match(y.use[,1],rownames(grp.const))
        idr2=match(y.use[,2],rownames(grp.const))
        dtr2=sapply(1:ncol(grp.const),
                    function(i)
                    {
                      sapply(1:nrow(y.use),function(j){paste(sort(c(grp.const[idr1[j],i],grp.const[idr2[j],i])),collapse = "__")})
                    })
        colnames(dtr2)=colnames(grp.const)
      }else{dtr2=matrix(nr=nrow(y.use),nc=0)}
      dtr=data.frame(dtr1,dtr2,stringsAsFactors = FALSE)
    }else{dtr=data.frame(lmm.dat.rand[,3:ncol(lmm.dat.rand),drop=FALSE],stringsAsFactors = FALSE)}
    
    formu=formula(sub("\\.",paste0(colnames(datin),collapse = "+"),deparse(lmm.formula)))
  }
  
  if(!forward)
  {
    sum.mod<-function(mod,dat.obs,dat.perm,silent=FALSE)
    {
      mods=summary(mod)
      if(grepl("lmerMod",class(mod))){coef.obs=fixef(mod)}else{coef.obs=mod$coefficients}
      t.obs=mods$coefficients[,"t value"]
      anv=anova(mod)
      
      if(class(mod)[1]=="lm")
      {
        SumSq=anv$`Sum Sq`[match(colnames(dat.obs),rownames(anv))]
        TSS=sum(anv$`Sum Sq`)
      }else if(class(mod)[1]=="glm"){
        SumSq=anv$Deviance[match(colnames(dat.obs),rownames(anv))]
        TSS=anv$`Resid. Dev`[1]
      }else if(grepl("lmerMod",class(mod))){
        SumSq=anv$`Sum Sq`[match(colnames(dat.obs),rownames(anv))]
        TSS=sum((mod@frame[,1]-mean(mod@frame[,1]))^2)
      }
      coef1=coef.obs[match(colnames(dat.obs),names(coef.obs))];names(coef1)=c()
      t.value1=t.obs[match(colnames(dat.obs),names(t.obs))];names(t.value1)=c()
      p.para1=mods$coefficients[,"Pr(>|t|)"][match(colnames(dat.obs),rownames(mods$coefficients))];names(p.para1)=c()
      out1=data.frame(coef=coef1,VariancePortion=SumSq/TSS,t.value=t.value1,
                      P.parametric=p.para1,stringsAsFactors = FALSE)
      R2.out=r2ma(mod)
      rownames(out1)=colnames(dat.obs)
      if(!is.null(dat.perm))
      {
        r2c.obs=c(model.R2=R2.out[[1]],t.obs)
        r2c.rand=sapply(1:ncol(dat.perm),
                        function(i)
                        {
                          if((!silent)&(i %in% seq(from=1,to=ncol(dat.perm),by=100))){message("---perm i=",i,".---",date())}
                          yri=dat.perm[,i] #y[permid[[i]]]
                          if(class(mod)[1]=="lm"){modri=lm(yri~.,data = dat.obs)}else if(class(mod)[1]=="glm"){
                            modri=glm(yri~.,data=dat.obs,family = mod$family)
                          }else if(grepl("lmerMod",class(mod))){
                            datri=mod@frame
                            datri[,1]=yri
                            modri=lmer(formula(mod),data=datri)
                          }
                          modsi=summary(modri)
                          tri=modsi$coefficients[,"t value"][match(names(t.obs),rownames(modsi$coefficients))];names(tri)=c()
                          c(r2ma(modri)[[1]],tri)
                        })
        p.perm=(rowSums(r2c.rand>=(r2c.obs-EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2c.rand))+1)
        p.neg=(rowSums(r2c.rand<=(r2c.obs+EPS),na.rm = TRUE)+1)/(rowSums(!is.na(r2c.rand))+1)
        p.perm[which(r2c.obs<0)]=p.neg[which(r2c.obs<0)]
        p.perm1=p.perm[match(colnames(dat.obs),names(r2c.obs))];names(p.perm1)=c()
        out1=data.frame(out1,P.perm=p.perm1,stringsAsFactors = FALSE)
        P.perm.out=p.perm[[1]]
      }else{P.perm.out=NA}
      list(summary=out1,mod=mod,R2=R2.out[[1]],R2adj=R2.out[[2]],AIC=AIC(mod),P.perm=P.perm.out)
    }
    
    if("lm" %in% method)
    {
      mod=lm(y~.,data=datin)
      res$lm=sum.mod(mod,dat.obs = datin,dat.perm = pdat,silent = silent)
    }
    if("glm" %in% method)
    {
      mod=glm(y~.,data=datin,family = glm.family)
      res$glm=sum.mod(mod,dat.obs = datin,dat.perm = pdat,silent = silent)
    }
    if("lmm" %in% method)
    {
      mod=lmer(formu,data = data.frame(datin,y=y,dtr,stringsAsFactors = FALSE))
      res$lmm=sum.mod(mod,dat.obs = datin, dat.perm = pdat, silent = silent)
    }
  }
  
  #####
  # 5 # if forward selection
  #####
  if(forward)
  {
    if("lm" %in% method)
    {
      res$lm=fw.mod(y=y,datin=datin,method = "lm",silent=silent,dat.perm = pdat,mod1.perm=mod1.perm)
    }
    if("glm" %in% method)
    {
      res$glm=fw.mod(y=y,datin=datin,method = "glm",glm.family=glm.family,silent=silent,dat.perm = pdat,mod1.perm=mod1.perm)
    }
    if("lmm" %in% method)
    {
      res$lmm=fw.mod(y=y,datin = datin,method = "lmm",lmm.formula = lmm.formula,lmm.datr=dtr,silent=silent,dat.perm=pdat,mod1.perm=mod1.perm)
    }
  }
  res
}

# multi-dimension constrained permutation used for modified MRM. extracted from above function for mcMRM
# author: 
mc.perm<-function(grp.rand,grp.const=NULL,id.2col=NULL,rand=1000,try.time=5)
{
  # multiple-dimension constrained permutation
  if(!is.null(grp.const))
  {
    if(sum(rownames(grp.rand)!=rownames(grp.const))>0){stop("grp.rand and grp.const must have the same ids in the same order.")}
  }
  if(!is.null(id.2col))
  {
    id2c=unique(as.vector(as.matrix(id.2col)))
    if(sum(!(id2c %in% rownames(grp.rand)))>0){stop("id.2col has some ids not defined in grp.rand.")}
  }
  permf<-function(randn,...)
  {
    t.lev<-list()
    if(is.null(grp.const))
    {
      for(i in 1:ncol(grp.rand))
      {
        message("permutating group type i=",i," name=",colnames(grp.rand)[i],". ",date())
        t.lev[[i]]=unique(grp.rand[,i])
        permi=permute::shuffleSet(length(t.lev[[i]]),randn)
        if(nrow(permi)<randn) permi=rbind(permi,1:ncol(permi))
        pti=sapply(1:nrow(permi),function(j){t.lev[[i]][permi[j,match(grp.rand[,i],t.lev[[i]])]]})
        if(i==1){pts=pti;ctsi=grp.rand[,i]}else{
          ctsi=sapply(1:nrow(grp.rand),function(j){paste(grp.rand[j,1:i],collapse = "__")})
          ptsi=lapply(1:ncol(pti),
                      function(j)
                      {
                        outj=sapply(1:ncol(pts),function(k){paste0(pts[,k],"__",pti[,j])})
                        cj=matrix((!(outj %in% ctsi)),nrow = nrow(outj),ncol=ncol(outj))
                        outj[,colSums(cj)==0,drop=FALSE]
                      })
          pts=Reduce(cbind,ptsi)
        }
      }
    }else{
      
      ctc=sapply(1:nrow(grp.const),function(i){paste(grp.const[i,],collapse = "__")})
      c.lev=unique(ctc)
      
      for(i in 1:ncol(grp.rand))
      {
        message("permutating group type i=",i," name=",colnames(grp.rand)[i],". ",date())
        ti.lev=unique(grp.rand[,i])
        
        ti.ct=sapply(1:length(ti.lev),
                     function(j)
                     {
                       paste(sort(unique(ctc[grp.rand[,i]==ti.lev[j]])),collapse = "__")
                     })
        tict.lev=unique(ti.ct)
        idj=lapply(1:length(tict.lev),function(j){which(ti.ct==tict.lev[j])})
        idjr=lapply(1:length(tict.lev),
                    function(j)
                    {
                      if(length(idj[[j]])==1){matrix(idj[[j]],1,1)}else{
                        permij=permute::shuffleSet(length(idj[[j]]),randn)
                        if(nrow(permij)<randn) permij=rbind(permij,1:ncol(permij))
                        matrix(idj[[j]][permij],nr=nrow(permij),nc=ncol(permij))
                      }
                    })
        prdn=prod(sapply(idjr,nrow))
        if(prdn>randn)
        {
          tpid=t(sapply(1:randn,
                        function(j)
                        {
                          outj=rep(0,length(ti.lev))
                          for(k in 1:length(idjr))
                          {
                            outj[idj[[k]]]=idjr[[k]][sample(nrow(idjr[[k]]),1),]
                          }
                          outj
                        }))
        }else{
          for(j in 1:length(tict.lev))
          {
            if(j==1){tpid=t(sapply(1:nrow(idjr[[1]]),function(k){out=rep(0,length(ti.lev));out[idj[[1]]]=idjr[[1]][k,];out}))}else{
              tpidl=lapply(1:nrow(idjr[[j]]),
                           function(k)
                           {
                             outk=tpid
                             outk[,idj[[j]]]=matrix(idjr[[j]][k,],nr=nrow(tpid),nc=ncol(idjr[[j]]),byrow = TRUE)
                             outk
                           })
              tpid=Reduce(rbind,tpidl)
            }
          }
        }
        permi=tpid
        pti=sapply(1:nrow(permi),function(j){ti.lev[permi[j,match(grp.rand[,i],ti.lev)]]})
        if(i==1){pts=sapply(1:ncol(pti),function(j){paste0(ctc,"__",pti[,j])});ctsi=paste0(ctc,"__",grp.rand[,i])}else{
          ctsi=sapply(1:nrow(grp.rand),function(j){paste(c(ctc[j],as.vector(as.matrix(grp.rand[j,1:i]))),collapse = "__")})
          ptsi=lapply(1:ncol(pti),
                      function(j)
                      {
                        outj=sapply(1:ncol(pts),function(k){paste0(pts[,k],"__",pti[,j])})
                        cj=matrix((!(outj %in% ctsi)),nrow = nrow(outj),ncol=ncol(outj))
                        outj[,colSums(cj)==0,drop=FALSE]
                      })
          pts=Reduce(cbind,ptsi)
        }
      }
    }
    
    ncol(pts)
    
    ct.lev=unique(ctsi)
    if(length(ct.lev)==length(ctsi))
    {
      pid=matrix(match(pts,ctsi),nr=nrow(pts),nc=ncol(pts))
    }else{
      message("the combination of grp info is not unique ID. may use random match when permutation. ",date())
      pid=sapply(1:ncol(pts),
                 function(j)
                 {
                   ptij=pts[,j]
                   outj=rep(NA,nrow(pts))
                   for(i in 1:length(ct.lev))
                   {
                     ctii=ct.lev[i]
                     idai=which(ctsi==ctii)
                     idij=which(ptij==ctii)
                     idaij=rep(idai,ceiling(length(idij)/length(idai)))
                     outj[idij]=sample(idaij,length(idij))
                   }
                   outj
                 })
    }
    
    id.rm=which(colSums(pid!=(1:nrow(pid)))==0)
    if(length(id.rm)>0){pid=pid[,-id.rm,drop=FALSE];endc=TRUE}else{endc=FALSE}
    
    if(!is.null(id.2col))
    {
      pdsn=0
      pds=matrix(nr=nrow(id.2col),nc=0)
      pid.res<-pid
      t=1
      while(pdsn<randn & ncol(pid.res)>0)
      {
        message("matching perm data t=",t,". ",date())
        if(ncol(pid.res)>randn)
        {
          idp.use=sample(ncol(pid.res),randn)
          pid=pid.res[,idp.use,drop=FALSE]
          pid.res=pid.res[,-idp.use,drop=FALSE]
        }else{
          pid=pid.res
          pid.res=matrix(nr=nrow(pid.res),nc=0)
        }
        id22=paste0(id.2col[,1],"__",id.2col[,2])
        pdst=sapply(1:ncol(pid),
                    function(i)
                    {
                      pidi=pid[,i]
                      pn2i=matrix(rownames(grp.rand)[pidi[match(as.matrix(id.2col[,1:2]),rownames(grp.rand))]],nr=nrow(id.2col),nc=2)
                      pn2c=paste0(pn2i[,1],"__",pn2i[,2])
                      pn2cb=paste0(pn2i[,2],"__",pn2i[,1])
                      out=match(pn2c,id22)
                      outb=match(pn2cb,id22)
                      out[is.na(out)]=outb[is.na(out)]
                      out
                    })
        pdst=pdst[,colSums(is.na(pdst))==0,drop=FALSE]
        pds=cbind(pds,pdst)
        pdsn=ncol(pds)
        t=t+1
      }
      pdat=pds
    }else{
      pdat=pid
    }
    list(pdat=pdat,endc=endc)
  }
  
  ptf=permf(randn=rand)
  pdat=ptf$pdat
  endc=ptf$endc
  
  if(ncol(pdat)>=rand)
  {
    pdat=pdat[,sample(ncol(pdat),rand),drop=FALSE]
  }else{
    t=2
    while(ncol(pdat)<rand & (!endc) & t<try.time)
    {
      message("Repeating perm t=",t,". ",date())
      ptf=permf(randn=rand*t)
      pdat=ptf$pdat
      endc=ptf$endc
      t=t+1
    }
    if(endc){message("All possible permuations are used.")}
    if(ncol(pdat)<rand){warning("Actual permutation time is ",ncol(pdat))}
    if(ncol(pdat)<=10){warning("possible permutation time is no more than 10.")}
  }
  colnames(pdat)=paste0("Perm",1:ncol(pdat))
  t(pdat)
}

mcMantel<-function(y.3col,y.ids=NULL,x.3col, z.3col=NULL, 
                   grp.rand,grp.const=NULL,try.time=5,
                   method = "pearson", permutations = 999, 
                   na.rm = FALSE, parallel = getOption("mc.cores"), 
                   tail=1)
{
  # multi-dimension constrained Mantel test and partial Mantel test
  # grp.rand defines the groups to permutate
  # grp.const defines the groups which should be maintained the same as observed.
  # version 2020.5.17
  # previous dependence: mc.perm, match.2col in {ieggr}
  #####
  # 1 # match IDs
  #####
  if(is.null(z.3col))
  {
    namec=iCAMP::match.2col(check.list = list(y=y.3col,x=x.3col))
    y3=namec$y
    x3=namec$x
  }else{
    namec=iCAMP::match.2col(check.list = list(y=y.3col,x=x.3col,z=z.3col))
    y3=namec$y
    x3=namec$x
    z3=namec$z
  }
  samps=unique(as.vector(as.matrix(y3[,1:2])))
  if(is.null(grp.rand)){grp.rand=data.frame(sample=samps,stringsAsFactors = FALSE);rownames(grp.rand)=samps}
  if(!is.null(grp.const))
  {
    spc=iCAMP::match.name(rn.list = list(grp.rand=grp.rand,grp.const=grp.const))
    grp.rand=spc$grp.rand
    grp.const=spc$grp.const
  }
  if(sum(!(samps %in% rownames(grp.rand)))){stop("all samples should be defined in grp.rand or grp.const.")}
  grp.rand=grp.rand[match(samps,rownames(grp.rand)),,drop=FALSE]
  if(!is.null(grp.const)){grp.const=grp.const[match(samps,rownames(grp.const)),,drop=FALSE]}
  
  id2y3=paste0(y3[,1],"__",y3[,2])
  if(is.null(y.ids))
  {
    y.ids=1:nrow(y3)
  }else{
    y.use=y.3col[y.ids,,drop=FALSE]
    id2use=paste0(y.use[,1],"__",y.use[,2])
    id2useb=paste0(y.use[,2],"__",y.use[,1])
    y.ids=match(id2use,id2y3)
    y.idsb=match(id2useb,id2y3)
    y.ids[is.na(y.ids)]=y.idsb[is.na(y.ids)]
  }
  
  
  #####
  # 2 # observed correlation
  #####
  
  if (na.rm) {use <- "complete.obs"} else{use <- "all.obs"}
  ryx <- cor(y3[y.ids,3], x3[y.ids,3], method = method, use = use)
  if(!is.null(z.3col))
  {
    ryz <- cor(y3[y.ids,3], z3[y.ids,3], method = method, use = use)
    rxz <- cor(x3[y.ids,3], z3[y.ids,3], method = method, use = use)
    part.cor <- function(ryx, ryz, rxz) {
      (ryx - ryz * rxz)/sqrt(1 - ryz * ryz)/sqrt(1 - rxz * rxz)
    }
    statistic <- part.cor(ryx, ryz, rxz)
  }else{
    statistic <- ryx
  }
  variant <- match.arg(method, eval(formals(cor)$method))
  variant <- switch(variant, pearson = "Pearson's product-moment correlation", 
                    kendall = "Kendall's rank correlation tau", spearman = "Spearman's rank correlation rho", 
                    variant)
  #####
  # 3 # data permutation
  #####
  if(is.na(statistic)){permutations=0}else{
    permat=mc.perm(grp.rand = grp.rand,grp.const = grp.const,id.2col = y3[,1:2,drop=FALSE],
                   rand=permutations,try.time = try.time)
    permutations=nrow(permat)
  }
  #####
  # 4 # calculate r after permutation
  #####
  
  if (permutations>1) {
    EPS <- sqrt(.Machine$double.eps)
    perm <- rep(0, permutations)
    ptest <- function(take, ...)
    {
      y3r=y3[take, ]
      ryrx <- cor(y3r[y.ids,3], x3[y.ids,3], method = method, use = use)
      if(!is.null(z.3col))
      {
        ryrz <- cor(y3r[y.ids,3], z3[y.ids,3], method = method, use = use)
        out <- part.cor(ryrx, ryrz, rxz)
      }else{
        out <- ryrx
      }
      out
    }
    
    if (is.null(parallel)){parallel <- 1}
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
      if (.Platform$OS.type == "unix" && !hasClus) {
        perm <- do.call(rbind, parallel::mclapply(1:permutations, 
                                                  function(i, ...) ptest(permat[i, ], ...), mc.cores = parallel))
      } else {
        if (!hasClus) {
          parallel <- parallel::makeCluster(parallel)
        }
        perm <- parallel::parRapply(parallel, permat, ptest)
        if (!hasClus) 
          parallel::stopCluster(parallel)
      }
    } else {
      perm <- sapply(1:permutations, function(i, ...) ptest(permat[i,], ...))
    }
    
    if(tail==2)
    {
      signif <- (sum((perm^2) >= ((statistic - EPS)^2)) + 1)/(permutations +1)
    }else if (tail==1){
      if(statistic>=0)
      {
        signif <- (sum(perm >= (statistic - EPS)) + 1)/(permutations +1)
      }else{
        signif <- (sum(perm <= (statistic - EPS)) + 1)/(permutations +1)
      }
    }else{signif=NA}
  }else {
    signif <- NA
    perm <- NULL
  }
  list(call = match.call(), method = variant, statistic = statistic, 
       signif = signif, perm = perm, permutations = permutations)
}


# End #