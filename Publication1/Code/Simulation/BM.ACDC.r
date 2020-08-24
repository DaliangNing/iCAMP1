BM.ACDC<-function(tree,mean.trait=0,sig2=1,bounds=c(-Inf,Inf),
               g=NULL,rand.num=1,nworker=4,code.wd=getwd())
{
  requireNamespace("iCAMP")
  requireNamespace("phytools")
  if(!is.null(g))
  {
    drt=iCAMP::tree.droot(tree,nworker = nworker)
    drt.id=drt$node
    drt.d=drt$distRoot
    drt.dg=(1-g^(-drt.d))/(1-g^(-1))
    treen=tree
    treen$edge.length=drt.dg[match(tree$edge[,2],drt.id)]-drt.dg[match(tree$edge[,1],drt.id)]
  }else{treen=tree}
  
  BMt=phytools::fastBM(treen, a=mean.trait, mu=0, sig2=sig2, bounds=bounds, internal=FALSE, nsim=rand.num)
  out=data.frame(BMt)
  colnames(out)=paste0("Rand",1:rand.num)
  list(trait=out,tree.new=treen)
}