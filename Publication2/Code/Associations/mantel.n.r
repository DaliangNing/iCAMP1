mantel.n<-function (xmatrix, ymatrix, method = "pearson", permutations = 999, strata = NULL,
          na.rm = FALSE, parallel = getOption("mc.cores"),tail=1,diag.in=FALSE) 
{
  if(!(tail %in% 1:2)){stop("Tail should be 1 or 2.")}
  EPS <- sqrt(.Machine$double.eps)
  xmatrix<-as.matrix(xmatrix)
  ymatrix<-as.matrix(ymatrix)
  
  if(diag.in)
  {
    xdis <- c(diag(xmatrix),as.vector(as.dist(xmatrix)))
    ydis <- c(diag(ymatrix),as.vector(as.dist(ymatrix)))
  }else{
    xdis<-as.vector(as.dist(xmatrix))
    ydis<-as.vector(as.dist(ymatrix))
  }
  
  if (na.rm)
  {
    use <- "complete.obs"
  }else{use <- "all.obs"}
  
  statistic <- cor(xdis, ydis, method = method, use = use)
  coeff=lm(xdis~ydis)$coefficients[[2]]
  variant <- match.arg(method, eval(formals(cor)$method))
  variant <- switch(variant, pearson = "Pearson's product-moment correlation", 
                    kendall = "Kendall's rank correlation tau", spearman = "Spearman's rank correlation rho", 
                    variant)
  N <- nrow(xmatrix)
  
  permat <-vegan:::getPermuteMatrix(permutations, N, strata = strata)
  
  if (ncol(permat) != N) stop(gettextf("'permutations' have %d columns, but data have %d observations", ncol(permat), N))
  
  permutations <- nrow(permat)
  
  if (permutations) {
    perm <- numeric(permutations)
    xmat <- xmatrix
    asdist <- row(xmat) > col(xmat)
    if(diag.in)
    {
      ptest <- function(take, ...) {
        xmat.p<-xmat[take, take]
        permvec <- c(diag(xmat.p),xmat.p[asdist])
        drop(cor(permvec, ydis, method = method, use = use))
      }
    }else{
      ptest <- function(take, ...) {
        xmat.p<-xmat[take, take]
        permvec <- xmat.p[asdist]
        drop(cor(permvec, ydis, method = method, use = use))
      }
    }
    
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
      if (.Platform$OS.type == "unix" && !hasClus) {
        perm <- do.call(rbind, parallel::mclapply(1:permutations, 
                                        function(i, ...) ptest(permat[i, ], ...), mc.cores = parallel))
      }else {
        if (!hasClus) {
          parallel <- parallel::makeCluster(parallel)
        }
        perm <- parallel::parRapply(parallel, permat, ptest)
        if (!hasClus) 
          parallel::stopCluster(parallel)
      }
    }else {
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
        signif <- (sum(perm <= (statistic + EPS)) + 1)/(permutations +1)
      }
    }else{signif=NA}
  }else {
    signif <- NA
    perm <- NULL
  }
  res <- list(call = match.call(), method = variant, statistic = statistic, coefficient=coeff,
              signif = signif, perm = perm, permutations = permutations, 
              control = attr(permat, "control"))
  res
}