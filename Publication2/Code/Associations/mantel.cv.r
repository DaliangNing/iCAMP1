mantel.cv<-function (xmatrix, ymatrix, method=c("k-fold","monte","repeat"),
                     k.m=5,mt.nv=0.1,mt.rand=100,rp.nv=mt.nv,rp.rand=100,
                     nv.perm=NULL) 
{
  EPS <- sqrt(.Machine$double.eps)
  xmatrix<-as.matrix(xmatrix)
  ymatrix<-as.matrix(ymatrix)
  
  idck=match.name(both.list = list(xmatrix=xmatrix,ymatrix=ymatrix))
  xmatrix=idck$xmatrix
  ymatrix=idck$ymatrix
  
  obn=nrow(xmatrix)
  if(obn<3){stop("Too few samples to do cross validate.")}
  if("k-fold" %in% method)
  {
    if(k.m>obn)
    {
      k.m=obn
      warning("k-fold number was more than observation number, thus set equal to observation number.")
    }
  }
  if("monte" %in% method)
  {
    if(round(mt.nv*obn)<1)
    {
      mt.nv=1/obn
      warning("The Monte Carlo test data sample size setting too low, thus changed to 1.")
    }
    if(round((1-mt.nv)*obn)<2)
    {
      mt.nv=1-(2/obn)
      warning("The Monte Carlo training data sample size setting too low, thus changed to 2.")
    }
  }
  if("repeat" %in% method)
  {
    if(round(rp.nv*obn)<1)
    {
      rp.nv=1/obn
      warning("The Repeated learning-testing test data sample size setting too low, thus changed to 1.")
    }
    if(round((1-rp.nv)*obn)<2)
    {
      rp.nv=1-(2/obn)
      warning("The Repeated learning-testing training data sample size setting too low, thus changed to 2.")
    }
  }
  # 1 # regular lm
  x3col=dist.3col(xmatrix)
  y3col=dist.3col(ymatrix)
  data=data.frame(x=x3col[,3],y=y3col[,3],stringsAsFactors = FALSE)
  gm=lm(y~x,data = data, y=TRUE)
  gms=summary(gm)
  
  r2.obs=1-(sum((gm$y-gm$fitted.values)^2)/sum((gm$y-mean(gm$y))^2));names(r2.obs)="Obs"
  fn=length(gm$coefficients)-1
  r2.adj=r2.obs-((1-r2.obs)*fn/(obn-fn-1));names(r2.adj)="Adjusted"
  r2.out=c(r2.obs,r2.adj)
  
  input.out<-cv.out<-r2cv.out<-vector()
  y.cvpred=list()
  
  if("k-fold" %in% method)
  {
    rand=sample(obn)%%k.m+1
    k.ms=sort(unique(rand))
    y.k=list();yk=1
    for(i in k.ms)
    {
      nv=which(rand==i)
      nt=which(rand!=i)
      idv=which((x3col[,1] %in% rownames(xmatrix)[nv])|(x3col[,2] %in% rownames(xmatrix)[nv]))
      data.t=data.frame(x=dist.3col(xmatrix[nt,nt,drop=FALSE])[,3],
                        y=dist.3col(ymatrix[nt,nt,drop=FALSE])[,3],
                        stringsAsFactors = FALSE)
      data.v=data.frame(x=x3col[idv,3],
                        y=y3col[idv,3],
                        stringsAsFactors = FALSE)
      gm.i=lm(y~x,data=data.t)
      y.k[[yk]]=data.frame(x=data.v$x,y.obs=data.v$y,y.predict=predict(gm.i,newdata = data.v),stringsAsFactors = FALSE)
      yk=yk+1
    }
    y.cvpred$kfold=data.frame(data.table::rbindlist(y.k),stringsAsFactors = FALSE)
    cv.k=(sum((y.cvpred$kfold$y.obs-y.cvpred$kfold$y.predict)^2))/nrow(y.cvpred$kfold)
    r2.cv.k=1-(sum((y.cvpred$kfold$y.obs-y.cvpred$kfold$y.predict)^2)/sum((y.cvpred$kfold$y.obs-mean(y.cvpred$kfold$y.obs))^2))
    names(cv.k)<-names(r2.cv.k)<-"kfold"
    cv.out=c(cv.out,cv.k);r2cv.out=c(r2cv.out,r2.cv.k)
    names(k.m)="kfolder.number";input.out=c(input.out,k.m)
  }
  
  if("monte" %in% method)
  {
    nvn=round(mt.nv*obn)
    ntn=obn-nvn
    meany=mean(y3col[,3])
    m.cv=sapply(1:mt.rand,
                function(i)
                {
                  nt=sort(sample(1:obn,ntn))
                  nv=sort((1:obn)[-nt])
                  idv=which((x3col[,1] %in% rownames(xmatrix)[nv])|(x3col[,2] %in% rownames(xmatrix)[nv]))
                  data.t=data.frame(x=dist.3col(xmatrix[nt,nt,drop=FALSE])[,3],
                                    y=dist.3col(ymatrix[nt,nt,drop=FALSE])[,3],
                                    stringsAsFactors = FALSE)
                  data.v=data.frame(x=x3col[idv,3],
                                    y=y3col[idv,3],
                                    stringsAsFactors = FALSE)
                  
                  gm.i=lm(y~x,data=data.t)
                  y.i=data.frame(x=data.v$x,y.obs=data.v$y,y.predict=predict(gm.i,newdata = data.v),stringsAsFactors = FALSE)
                  cvi=sum((y.i$y.obs-y.i$y.predict)^2)/nrow(y.i)
                  toti=sum((y.i$y.obs-meany)^2)/nrow(y.i)
                  c(cvi,toti,y.i$x,y.i$y.obs,y.i$y.predict)
                })
    cv.m=mean(m.cv[1,])
    r2.cv.m=1-(cv.m/(mean(m.cv[2,])))
    names(cv.m)<-names(r2.cv.m)<-"ManteCarlo"
    nvnm=(nrow(m.cv)-2)/3
    y.cvpred$MonteCarlo=data.frame(x=as.vector(m.cv[3:(nvnm+2),]),
                                   y.obs=as.vector(m.cv[(nvnm+3):(2*nvnm+2),]),
                                   y.predict=as.vector(m.cv[(2*nvnm+3):(3*nvnm+2),]),
                                   stringsAsFactors = FALSE)
    cv.out=c(cv.out,cv.m);r2cv.out=c(r2cv.out,r2.cv.m)
    names(mt.nv)="ManteCarlo.nv";names(mt.rand)="ManteCarlo.rand.time";input.out=c(input.out,mt.nv,mt.rand)
  }
  
  if("repeat" %in% method)
  {
    nvn=round(rp.nv*obn)
    ntn=obn-nvn
    if(!is.null(nv.perm))
    {
      if(nrow(nv.perm)<rp.rand)
      {
        nv.perm=NULL
        warning("Assinged nv.perm is not enough for randomizaiton.")
      }
    }
    
    if(is.null(nv.perm))
    {
      if(choose(obn,nvn)<rp.rand)
      {
        nv.perm=t(combn(obn,nvn))
        rp.rand=nrow(nv.perm)
      }else if(choose(obn,nvn)<(10^5)){
        nv.perm=t(combn(obn,nvn)[,sample(choose(obn,nvn),rp.rand),drop=FALSE])
      }else{
        nv.rand=matrix(replicate(rp.rand*2,sort(sample(obn,nvn))),nr=rp.rand*2,byrow=TRUE)
        nv.perm=nv.rand[!duplicated(nv.rand),,drop=FALSE]
        while(nrow(nv.perm)<rp.rand)
        {
          nv.rand=rbind(nv.rand,matrix(replicate(rp.rand*2,sort(sample(obn,nvn))),nr=rp.rand*2,byrow=TRUE))
          nv.perm=nv.rand[!duplicated(nv.rand),,drop=FALSE]
        }
        nv.perm=nv.perm[1:rp.rand,,drop=FALSE]
      }
    }else{
      nv.perm=nv.perm[1:rp.rand,,drop=FALSE]
    }
    
    rp.cv=sapply(1:rp.rand,
                 function(i)
                 {
                   nv=sort(nv.perm[i,])
                   nt=sort((1:obn)[-nv])
                   
                   idv=which((x3col[,1] %in% rownames(xmatrix)[nv])|(x3col[,2] %in% rownames(xmatrix)[nv]))
                   data.t=data.frame(x=dist.3col(xmatrix[nt,nt,drop=FALSE])[,3],
                                     y=dist.3col(ymatrix[nt,nt,drop=FALSE])[,3],
                                     stringsAsFactors = FALSE)
                   data.v=data.frame(x=x3col[idv,3],
                                     y=y3col[idv,3],
                                     stringsAsFactors = FALSE)
                   
                   gm.i=lm(y~x,data=data.t)
                   y.i=data.frame(x=data.v$x,y.obs=data.v$y,y.predict=predict(gm.i,newdata = data.v),stringsAsFactors = FALSE)
                   cvi=sum((y.i$y.obs-y.i$y.predict)^2)/nrow(y.i)
                   toti=sum((y.i$y.obs-meany)^2)/nrow(y.i)
                   c(cvi,toti,y.i$x,y.i$y.obs,y.i$y.predict)
                 })
    cv.rp=mean(rp.cv[1,])
    r2.cv.rp=1-(cv.rp/(mean(rp.cv[2,])))
    names(cv.rp)<-names(r2.cv.rp)<-"RepeatLearning"
    nvnm=(nrow(rp.cv)-2)/3
    y.cvpred$RepeatLearning=data.frame(x=as.vector(rp.cv[3:(nvnm+2),]),
                                       y.obs=as.vector(rp.cv[(nvnm+3):(2*nvnm+2),]),
                                       y.predict=as.vector(rp.cv[(2*nvnm+3):(3*nvnm+2),]),
                                       stringsAsFactors = FALSE)
    cv.out=c(cv.out,cv.rp);r2cv.out=c(r2cv.out,r2.cv.rp)
    names(rp.nv)="RepeatLearning.nv";names(rp.rand)="RepeatLearning.rand.time";input.out=c(input.out,rp.nv,rp.rand)
  }
  list(coefficients=gm$coefficients,R2=r2.out,R2cv=r2cv.out,CVE=cv.out,inputs=input.out,
       pred.lm=data.frame(y=gm$y,y.pred=gm$fitted.values),pred.cv=y.cvpred)
}
