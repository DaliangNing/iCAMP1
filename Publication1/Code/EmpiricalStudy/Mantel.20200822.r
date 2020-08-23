
rm(list=ls())

###########################
# 1 # file/folder paths/names and package loading
# you may need to change them to the paths on your computer before testing the code.

wdm="./Data/EmpiricalData" # the folder saving input file
icres.file="OKCW5Y.T1R1000.iCAMP.bNRIiRCa.sum.csv"
treat.file="OKCW5Y.Treatment.csv"
env.file="OKCW5Y.geoenv.tempc.csv"
xy.file="OKCW5Y.Coordinates.csv"
prefix="OKCW5Y"
wd=wdm # the folder to save the results
tool.code.wd="./Code/Tools" # the folder saving tools.r

######
source(paste0(tool.code.wd,"/tools.r"))
library(vegan)

###########################
# 2 # load data
icres=lazyopen(icres.file)
head(icres)
treat=lazyopen(treat.file)
head(treat)
env.in=lazyopen(env.file)
xy=lazyopen(xy.file)

###########################
# 3 # transform and select

# geographic distance and PCNM
gdis=dist(xy)
gdis3c=dist.3col(gdis)

pcnmg=vegan::pcnm(gdis)
pcnm.perc=pcnmg$values[pcnmg$values>0]/sum(pcnmg$values[pcnmg$values>0])
sum(pcnm.perc[1:5])
pcnms=pcnmg$vectors[,1:5,drop=FALSE]

env=data.frame(env.in,xy,pcnms,stringsAsFactors = FALSE)

# log transform
logtr<-function(m)
{
  if(min(m)<=0)
  {
    newm=m-min(m)
    newm[which(newm==0)]=min(newm[newm>0])/20
    out=log(newm)
  }else{out=log(m)}
  out
}
envlog=sapply(1:ncol(env),function(i){logtr(env[,i])})
colnames(envlog)=colnames(env)
rownames(envlog)=rownames(env)
envlog[,"pH"]=env[,"pH"]

#################
# use envlog or env as env.use to set log-transform or not
env.use=envlog
# env.use=env

# use "H" or "N" to test in warming or control plots.
trt.sel="H"
pron.id=5

#################
# to focus on within-treatment within-year
yr1=treat$year[match(icres$sample1,rownames(treat))]
yr2=treat$year[match(icres$sample2,rownames(treat))]
tr1=treat$treatment[match(icres$sample1,rownames(treat))]
tr2=treat$treatment[match(icres$sample2,rownames(treat))]

# choose treatment and process and grouping information.
id.sel=which(tr1==trt.sel & tr2==trt.sel)
(pronn=colnames(icres)[pron.id+2])
y.3col=icres[id.sel,c(1:2,pron.id+2),drop=FALSE]
dim(y.3col)
y.ids=which(yr1[id.sel]==yr2[id.sel])
samps=unique(as.vector(as.matrix(y.3col[,1:2])))
treatu=treat[match(samps,rownames(treat)),,drop=FALSE]
head(treatu)
grp.rand=treatu[,c("plot","year"),drop=FALSE]


######################
menv=(env.use[match(icres$sample1,rownames(env.use)),]+
        env.use[match(icres$sample2,rownames(env.use)),])/2
colnames(menv)=paste0("m",colnames(env.use))
denv=abs(env.use[match(icres$sample1,rownames(env.use)),]-
           env.use[match(icres$sample2,rownames(env.use)),])
colnames(denv)=paste0("d",colnames(env.use))
mdenv=cbind(icres[,1:2],menv,denv)
head(mdenv)
dim(mdenv)
head(mdenv[id.sel,])


###########################
# 4 # Mantel test. mc here means to use multi-dimension constrained permutation.
# linear model
mcm=t(sapply(3:ncol(mdenv),
             function(i)
             {
               x.3col=mdenv[id.sel,c(1,2,i)]
               mci=mcMantel(y.3col = y.3col,y.ids = y.ids,
                            x.3col = x.3col, grp.rand = grp.rand,
                            grp.const = NULL,try.time = 5,
                            method = "pearson",permutations = 999)
               outi=c(r=mci$statistic,p=mci$signif)
             }))
rownames(mcm)=colnames(mdenv)[3:ncol(mdenv)]
colnames(mcm)=c("r","P")
save.file(pmcm,prefix = prefix,filename = paste0("mcMantel-LM.",trt.sel,".",pronn),folder =wd)

# generalized linear model. here i use the funciton mcMRM to get GLM-based Mantel test results.
mrm.glm=mcMRM(y.3col=y.3col,data.x=mdenv[id.sel,,drop=FALSE],
             y.ids=y.ids,grp.rand=grp.rand,grp.const=NULL,
             scale.yn=TRUE,method="glm",forward=TRUE,mod1.perm=TRUE,
             glm.family=quasibinomial,lmm.formula=NULL,lmm.dat.rand=NULL,
             rand=1000,silent=FALSE)
mcmglm=mrm.glm$glm$mod1
save.file(pmcm,prefix = prefix,filename = paste0("mcMantel-GLM.",trt.sel,".",pronn),folder =wd)

###########################
# 5 # Partial Mantel test. 

pmcm=sapply(3:ncol(mdenv),
            function(i)
            {
              x.3col=mdenv[id.sel,c(1,2,i)]
              outi=sapply(3:ncol(mdenv),
                          function(j)
                          {
                            message("----Now i=",i," j=",j,". ",date())
                            if(i!=j)
                            {
                              z.3col=mdenv[id.sel,c(1,2,j)]
                              mcij=mcMantel(y.3col = y.3col,y.ids = y.ids,
                                            x.3col = x.3col,z.3col = z.3col,
                                            grp.rand = grp.rand,grp.const = NULL,try.time = 5,
                                            method = "pearson",permutations = 999)
                              outij=c(r=mcij$statistic,p=mcij$signif)
                            }else{
                              mcij=mcMantel(y.3col = y.3col,y.ids = y.ids,
                                            x.3col = x.3col,z.3col = NULL,
                                            grp.rand = grp.rand,grp.const = NULL,try.time = 5,
                                            method = "pearson",permutations = 999)
                              outij=c(r=mcij$statistic,p=mcij$signif)
                            }
                            outij
                          })
              outii=as.vector(t(outi))
              names(outii)=paste0(rep(colnames(mdenv)[3:ncol(mdenv)],2),rep(c(".r",".p"),each=ncol(mdenv)-2))
              outii
            })
colnames(pmcm)=colnames(mdenv)[3:ncol(mdenv)]
head(pmcm)
save.file(pmcm,prefix = prefix,filename = paste0("pmcMantel.",trt.sel,".",pronn),folder =wd)

# End #