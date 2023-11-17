cv.PLS<-function(Y,X,nv=0.1,rand=100,cvmethod="monte",k.m=10)
{
  library(ropls)
  if(cvmethod=="monte")
  {
  n.test=round(length(Y)*nv)
  n.train=length(Y)-n.test
  trsq=seq(from=1,to=rand,by=10)
  datt=lapply(1:rand,
              function(j)
              {
                if(j %in% trsq) {message("j=",j,". ",date())}
                id.train=sort(sample(length(Y),size = n.train))
                id.test=which(!((1:length(Y)) %in% id.train))
                modj=try(ropls::opls(x=X[id.train,,drop=FALSE],y=Y[id.train],fig.pdfC='none',info.txtC='none',predI=NA,orthoI=0,permI = 20),silent = TRUE)
                if(class(modj)=="try-error"){modj=try(ropls::opls(x=X[id.train,,drop=FALSE],y=Y[id.train],fig.pdfC='none',info.txtC='none',predI=1,orthoI=0,permI = 20))}
                Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
              })
  }else if(cvmethod=="k-fold"){
    n.kme=round(length(Y)/k.m)
    rand=integer(length(Y))
    rand[order(Y)]=ceiling((1:length(Y))/n.kme)
    k.ms=sort(unique(rand))
    datt=lapply(k.ms,
                function(j)
                {
                  id.train=which(rand!=j)
                  id.test=which(rand==j)
                  modj=try(ropls::opls(x=X[id.train,,drop=FALSE],y=Y[id.train],fig.pdfC='none',info.txtC='none',predI=NA,orthoI=0,permI = 20),silent = TRUE)
                  if(class(modj)=="try-error"){modj=try(ropls::opls(x=X[id.train,,drop=FALSE],y=Y[id.train],fig.pdfC='none',info.txtC='none',predI=1,orthoI=0,permI = 20))}
                  Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                  data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
                })
  }
  datm=Reduce(rbind,datt)
  rownames(datm)=c()
  datm=datm[which(rowSums(is.na(datm))==0),,drop=FALSE]
  R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
  moda=try(ropls::opls(x=X,y=Y,fig.pdfC='none',info.txtC='none',predI=NA,orthoI=0,permI = 20),silent = TRUE)
  if(class(moda)=="try-error"){moda=try(ropls::opls(x=X,y=Y,fig.pdfC='none',info.txtC='none',predI=1,orthoI=0,permI = 20))}
  Y.moda=predict(moda)
  list(R2cv=R2fun(datm$obs,datm$test),R2=R2fun(Y,Y.moda),detail=datm,
       detail.mod=data.frame(ID=1:length(Y),obs=Y,test=Y.moda),Y=Y,X=X,mod=moda)
}

cv.RF<-function(Y,X,nv=0.1,rand=100,cvmethod="monte",k.m=10)
{
  library(randomForest)
  if(cvmethod=="monte")
  {
    n.test=round(length(Y)*nv)
    n.train=length(Y)-n.test
    trsq=seq(from=1,to=rand,by=10)
    datt=lapply(1:rand,
                function(j)
                {
                  if(j %in% trsq) {message("j=",j,". ",date())}
                  id.train=sort(sample(length(Y),size = n.train))
                  id.test=which(!((1:length(Y)) %in% id.train))
                  modj=randomForest::randomForest(x=X[id.train,,drop=FALSE],y=Y[id.train])
                  Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                  data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
                })
  }else if(cvmethod=="k-fold"){
    n.kme=round(length(Y)/k.m)
    rand=integer(length(Y))
    rand[order(Y)]=ceiling((1:length(Y))/n.kme)
    k.ms=sort(unique(rand))
    datt=lapply(k.ms,
                function(j)
                {
                  id.train=which(rand!=j)
                  id.test=which(rand==j)
                  modj=randomForest::randomForest(x=X[id.train,,drop=FALSE],y=Y[id.train])
                  Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                  data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
                })
  }
  datm=Reduce(rbind,datt)
  rownames(datm)=c()
  datm=datm[which(rowSums(is.na(datm))==0),,drop=FALSE]
  R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
  moda=randomForest::randomForest(x=X,y=Y)
  Y.moda=predict(moda)
  list(R2cv=R2fun(datm$obs,datm$test),R2=R2fun(Y,Y.moda),detail=datm,
       detail.mod=data.frame(ID=1:length(Y),obs=Y,test=Y.moda),Y=Y,X=X,mod=moda)
}

cv.LM<-function(Y,X,nv=0.1,rand=100,silent=TRUE,cvmethod="monte",k.m=10)
{
  library(randomForest)
  if(cvmethod=="monte")
  {
  n.test=round(length(Y)*nv)
  n.train=length(Y)-n.test
  trsq=seq(from=1,to=rand,by=10)
  datt=lapply(1:rand,
              function(j)
              {
                if((!silent) & (j %in% trsq)) {message("j=",j,". ",date())}
                id.train=sort(sample(length(Y),size = n.train))
                id.test=which(!((1:length(Y)) %in% id.train))
                modj=lm(Yt~.,data=data.frame(Yt=Y[id.train],X[id.train,,drop=FALSE],check.names = FALSE))
                Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
              })
  }else if(cvmethod=="k-fold"){
    n.kme=round(length(Y)/k.m)
    rand=integer(length(Y))
    rand[order(Y)]=ceiling((1:length(Y))/n.kme)
    k.ms=sort(unique(rand))
    datt=lapply(k.ms,
                function(j)
                {
                  id.train=which(rand!=j)
                  id.test=which(rand==j)
                  modj=lm(Yt~.,data=data.frame(Yt=Y[id.train],X[id.train,,drop=FALSE],check.names = FALSE))
                  Y.test=predict(modj,newdata=X[id.test,,drop=FALSE])
                  data.frame(ID=id.test,obs=Y[id.test],test=Y.test)
                })
  }
  datm=Reduce(rbind,datt)
  rownames(datm)=c()
  datm=datm[which(rowSums(is.na(datm))==0),,drop=FALSE]
  R2fun<-function(yv,fv){1-sum((fv-yv)^2)/sum((yv-mean(yv))^2)}
  moda=lm(Yt~.,data=data.frame(Yt=Y,X,check.names = FALSE))
  Y.moda=predict(moda)
  list(R2cv=R2fun(datm$obs,datm$test),R2=R2fun(Y,Y.moda),detail=datm,
       detail.mod=data.frame(ID=1:length(Y),obs=Y,test=Y.moda),Y=Y,X=X,mod=moda)
}
