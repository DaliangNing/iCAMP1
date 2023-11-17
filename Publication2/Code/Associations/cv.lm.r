cv.lm<-function(form.lm,data,method=c("k-fold","monte","repeat"),
                k.m=5,mt.nv=0.1,mt.rand=100,rp.nv=mt.nv,rp.rand=100,
                nv.perm=NULL, perm.test.rand=999,norm.test=FALSE,
                indep.test=FALSE,out.list=TRUE)
{
  # cross-validated error and R square (coefficient of determination) 
  # Daliang Ning 2017.4.1
  obn=nrow(data)
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
  gm=lm(formula = form.lm,data = data,y=TRUE)
  gms=summary(gm)
  
  ttest.p=gms$coefficients[,4];names(ttest.p)=paste0(names(ttest.p),".p.ttest")
  aic.obs=AIC(gm);names(aic.obs)="lm.AIC"
  r2.obs=1-(sum((gm$y-gm$fitted.values)^2)/sum((gm$y-mean(gm$y))^2));names(r2.obs)="Obs"
  fn=length(gm$coefficients)-1
  r2.adj=r2.obs-((1-r2.obs)*fn/(obn-fn-1));names(r2.adj)="Adjusted"
  r2.out=c(r2.obs,r2.adj)
  
  if(!is.null(perm.test.rand))
  {
    perm.n=permute::shuffleSet(obn,perm.test.rand)
    perm.test.rand=nrow(perm.n)
    datap=model.frame(form.lm,data=data)
    r2.rand=sapply(1:perm.test.rand,
                    function(i)
                    {
                      datai=datap
                      datai[,1]=datap[perm.n[i,],1]
                      summary(lm(formula = form.lm,data=datai))$r.squared
                    })
    EPS <-.Machine$double.eps
    p.perm=(sum(r2.rand>=(r2.obs-EPS))+1)/(perm.test.rand+1)
  }else{p.perm=NA}
  names(p.perm)="lm.p.perm"
  
  if(is.null(gms$fstatistic)){model.p=NA}else{
    model.p=pf(gms$fstatistic[1],gms$fstatistic[2],gms$fstatistic[3],lower.tail = FALSE)
  }
  names(model.p)="lm.p.Ftest"

  ttest.p=c(aic.obs,model.p,p.perm,ttest.p)
  
  
  if(norm.test)
  {
    # normality test
    resid1=gm$residuals
    norm.p=shapiro.test(x=resid1)$p.value;names(norm.p)="normal.test";
    ttest.p=c(ttest.p,norm.p)
  }
  if(indep.test)
  {
    # independence test
    dw.p=car::durbinWatsonTest(gm,max.lag=(obn-1))$p
    names(dw.p)=paste0("Independ.test.lag.",1:length(dw.p))
    ttest.p=c(ttest.p,dw.p)
  }
  
  input.out<-cv.out<-r2cv.out<-vector()
  y.cvpred=list()
  if("k-fold" %in% method)
  {
    rand=sample(obn)%%k.m+1
    k.ms=sort(unique(rand))
    y.k=rep(0,obn)
    for(i in k.ms)
    {
      nv=which(rand==i)
      nt=which(rand!=i)
      gm.i=lm(formula = form.lm,data=data[nt,,drop=FALSE])
      y.k[nv]=predict(gm.i,newdata = data[nv,,drop=FALSE])
    }
    y.cvpred$kfold=data.frame(y=gm$y,cvpred=y.k)
    cv.k=(sum((gm$y-y.k)^2))/obn
    r2.cv.k=1-(sum((gm$y-y.k)^2)/sum((gm$y-mean(gm$y))^2))
    names(cv.k)<-names(r2.cv.k)<-"kfold"
    cv.out=c(cv.out,cv.k);r2cv.out=c(r2cv.out,r2.cv.k)
    names(k.m)="kfolder.number";input.out=c(input.out,k.m)
  }
  
  if("monte" %in% method)
  {
    nvn=round(mt.nv*obn)
    ntn=obn-nvn
    m.cv=sapply(1:mt.rand,
                function(i)
                {
                  nt=sort(sample(1:obn,ntn))
                  nv=sort((1:obn)[-nt])
                  gm.i=lm(formula = form.lm,data=data[nt,,drop=FALSE])
                  y.i=predict(gm.i,newdata = data[nv,,drop=FALSE])
                  cvi=sum((gm$y[nv]-y.i)^2)/nvn
                  toti=sum((gm$y[nv]-mean(gm$y))^2)/nvn
                  c(cvi,toti,gm$y[nv],y.i)
                })
    cv.m=mean(m.cv[1,])
    r2.cv.m=1-(cv.m/(mean(m.cv[2,])));names(cv.m)<-names(r2.cv.m)<-"ManteCarlo"
    y.cvpred$MonteCarlo=data.frame(y=as.vector(m.cv[3:(nvn+2),]),cvpred=as.vector(m.cv[(nvn+3):(2*nvn+2),]))
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
                   gm.i=lm(formula = form.lm,data=data[nt,,drop=FALSE])
                   y.i=predict(gm.i,newdata = data[nv,,drop=FALSE])
                   cvi=sum((gm$y[nv]-y.i)^2)/nvn
                   toti=sum((gm$y[nv]-mean(gm$y))^2)/nvn
                   c(cvi,toti,gm$y[nv],y.i)
                 })
    cv.rp=mean(rp.cv[1,])
    r2.cv.rp=1-(cv.rp/(mean(rp.cv[2,])));names(cv.rp)<-names(r2.cv.rp)<-"RepeatLearning"
    y.cvpred$RepeatLearning=data.frame(y=as.vector(rp.cv[3:(nvn+2),]),cvpred=as.vector(rp.cv[(nvn+3):(2*nvn+2),]))
    cv.out=c(cv.out,cv.rp);r2cv.out=c(r2cv.out,r2.cv.rp)
    names(rp.nv)="RepeatLearning.nv";names(rp.rand)="RepeatLearning.rand.time";input.out=c(input.out,rp.nv,rp.rand)
  }
  
  res=list(coefficients=gm$coefficients,R2=r2.out,R2cv=r2cv.out,CVE=cv.out,p.values=ttest.p,inputs=input.out)
  
  if(!out.list)
  {
    res=unlist(res,use.names = TRUE)
  }else{
    res=c(res,list(pred.lm=data.frame(y=gm$y,y.pred=gm$fitted.values),pred.cv=y.cvpred))
  }
  res
}
