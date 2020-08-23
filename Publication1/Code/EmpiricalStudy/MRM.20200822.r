
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
library(MASS)

###########################
# 2 # load data
icres=lazyopen(icres.file)
head(icres)
treat=lazyopen(treat.file)
head(treat)
xy=lazyopen(xy.file)
gdis=dist(xy)
gdis3c=dist.3col(gdis)


env.in=lazyopen(env.file)
env=env.in
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
env.use=envlog
#env.use=env # to compare the results with and without log transformation.

yr1=treat$year[match(icres$sample1,rownames(treat))]
yr2=treat$year[match(icres$sample2,rownames(treat))]
tr1=treat$treatment[match(icres$sample1,rownames(treat))]
tr2=treat$treatment[match(icres$sample2,rownames(treat))]

trt.sel="H"
pron.id=5

# choose treatment and process and grouping information.
id.sel=which(tr1==trt.sel & tr2==trt.sel)
pronn=colnames(icres)[pron.id+2]
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

mrm.lm=mcMRM(y.3col=y.3col,data.x=mdenv[id.sel,,drop=FALSE],
             y.ids=y.ids,grp.rand=grp.rand,grp.const=NULL,
             scale.yn=TRUE,method="lm",forward=TRUE,mod1.perm=FALSE,
             glm.family=quasibinomial,lmm.formula=NULL,lmm.dat.rand=NULL,
             rand=1000,silent=FALSE)
save(mrm.lm,file = "")
(mrmsum=mrm.lm$lm$summary)
f.mrm=mrmsum$Add.Factor

mdenvu=mdenv[id.sel[y.ids],match(f.mrm,colnames(mdenv))]
mdenvus=data.frame(scale(mdenvu))

#datalm=data.frame(y=as.vector(scale(y.3col[y.ids,3])),mdenvus,stringsAsFactors = FALSE)
datalm=data.frame(y=y.3col[y.ids,3],mdenvus,stringsAsFactors = FALSE)
mod0<-lm(y~1,data=datalm)
full.model <- lm(y~.,data = datalm)
# Stepwise regression model
step.model <- stepAIC(mod0,scope = list(upper=full.model,lower=mod0), direction = "forward", 
                      trace = FALSE)
summary(step.model)
anova(step.model)

data.sel=step.model$model
mod.obs=lm(y~.,data=data.sel)
anova(mod.obs)
t.obs=summary(mod.obs)$coefficients[-1,3]
f.obs=anova(mod.obs)$`F value`[1:(ncol(data.sel)-1)]
R2.obs=summary(mod.obs)$r.squared

# significance test need to use permutaiton strategy defined by MRM
permut=mc.perm(grp.rand = grp.rand,grp.const = NULL,id.2col = y.3col[,1:2],rand = 1000)
dim(permut)

f.rand=sapply(1:nrow(permut),
              function(i)
              {
                yr=y.3col[permut[i,][y.ids],3]
                datar=data.sel
                datar[,1]=yr
                modr=lm(y~.,data=datar)
                anova(modr)$`F value`[1:(ncol(datar)-1)]
              })
data.frame(factor=colnames(data.sel)[-1],p.perm=rowSums(f.rand>matrix(f.obs,nrow = nrow(f.rand),ncol(f.rand)))/ncol(f.rand))

R2.rand=sapply(1:nrow(permut),
               function(i)
               {
                 yr=y.3col[permut[i,][y.ids],3]
                 datar=data.sel
                 datar[,1]=yr
                 modr=lm(y~.,data=datar)
                 summary(modr)$r.squared
               })
(p.mod=sum(R2.rand>=R2.obs)/length(R2.rand))


t.rand=sapply(1:nrow(permut),
              function(i)
              {
                yr=y.3col[permut[i,][y.ids],3]
                datar=data.sel
                datar[,1]=yr
                modr=lm(y~.,data=datar)
                summary(modr)$coefficients[-1,3]
              })
data.frame(p.perm=rowSums(abs(t.rand)>matrix(abs(t.obs),nrow = nrow(t.rand),ncol(t.rand)))/ncol(t.rand))

# End #