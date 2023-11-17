# 1 # specify files and folders
wd="/Code/MissingDataFill"
env.file="100Well.GeoENV.Raw.csv"
prefix="100Well"
Rtoolfile="/Code/handytool.r"

# 2 # load data
setwd(wd)
source(Rtoolfile)
envold=lazyopen(env.file)
dim(envold)
head(envold)

mis.num=colSums(is.na(envold))
plot(mis.num)
sum(mis.num>0)

env.ful=envold[,which(mis.num==0),drop=FALSE]
dim(env.ful);sum(is.na(env.ful))

id.mis=which(mis.num>0)
mis.nor=mis.num[order(mis.num)]
mis.norp=mis.nor[mis.nor>0]

r.rec=numeric(length(mis.norp))
f.name=character(length(mis.norp))
j=1
j.max=length(mis.norp)
lm.list=list()

for(j in 1:j.max)
{
  fi<-id.na<-lm.sel<-list()
  lm.r2=numeric(length(mis.norp))
  
  for(i in 1:length(mis.norp))
  {
    message("now j=",j,", i=",i," in ",length(mis.norp),". ",date())
    id=which(colnames(envold)==names(mis.norp)[i])
    fi[[i]]=envold[,id]
    id.na[[i]]=which(is.na(fi[[i]]))
    fi.t=fi[[i]][-id.na[[i]]]
    env.ful.t=env.ful[-id.na[[i]],,drop=FALSE]
    # 1 # begin with all-detected factors
    lm0=lm(fi.t~1,data=env.ful.t)
    lm1=lm(fi.t~.,data=env.ful.t)
    lm.sel[[i]]=step(lm0,scope = list(lower=lm0,upper=lm1),direction = "both",trace = 0)
    lm.sum=summary(lm.sel[[i]])
    lm.r2[i]=lm.sum$adj.r.squared
  }
  
  idm=match(max(lm.r2),lm.r2)
  r.rec[j]=lm.r2[idm]
  f.name[j]=names(mis.norp)[idm]
  
  lm.list[[j]]=lm.sel[[idm]]
  
  fact.sel=names(lm.sel[[idm]]$coefficients[-1])
  fact.na=env.ful[id.na[[idm]],match(fact.sel,colnames(env.ful)),drop=FALSE]
  pred.na=(as.matrix(fact.na) %*% matrix(lm.sel[[idm]]$coefficients[-1],nc=1))+lm.sel[[idm]]$coefficients[1]
  pred.na[pred.na<0]=0
  fi[[idm]][id.na[[idm]]]=pred.na
  env.ful=cbind(env.ful,fi[[idm]])
  colnames(env.ful)[ncol(env.ful)]=f.name[j]
  mis.norp=mis.norp[-idm]
}

r.rec
save.file(data.frame(adjust.r2=r.rec),prefix = prefix,filename = "r2.recode",folder = wd)

env.res=env.ful[,match(colnames(envold),colnames(env.ful))]
save.file(env.res,prefix = prefix,filename = "env.ful",folder = wd)

length(mis.norp)
mis.norp2=mis.nor[mis.nor>0]
length(mis.norp2)
names(lm.list)<-names(mis.norp2)
save(lm.list,file=paste(wd,"/",prefix,".env.fill_miss.lm.detail.rda",sep=""))
