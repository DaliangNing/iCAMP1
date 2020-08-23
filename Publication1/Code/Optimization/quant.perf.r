quant.perf<-function(expm,obsm,pron=c("HoS","HeS","HD","DL","DR"),RA=NULL)
{
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
  
  ccm=sapply(1:ncol(expm),
             function(i)
             {
               ccc(x=expm[,i],y=obsm[,i],out.vector = TRUE)[1:2]
             })
  cct=ccc(x=as.vector(as.matrix(expm)),y=as.vector(as.matrix(obsm)),out.vector = TRUE)[1:2]
  out=cbind(ccm,cct)
  colnames(out)=c(pron,"Total")
  rownames(out)=c("qACC","qPRC")
  if(!is.null(RA))
  {
    ccm2=sapply(1:ncol(expm),
               function(i)
               {
                 ccc(x=expm[,i]*RA,y=obsm[,i]*RA,out.vector = TRUE)[1:2]
               })
    cct2=ccc(x=as.vector(as.matrix(expm)*RA),y=as.vector(as.matrix(obsm)*RA),out.vector = TRUE)[1:2]
    out2=cbind(ccm2,cct2)
    colnames(out2)=c(pron,"Total")
    rownames(out2)=c("qACC.RA","qPRC.RA")
    out=rbind(out,out2)
  }
  out
}