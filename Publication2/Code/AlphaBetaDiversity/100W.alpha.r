rm(list=ls())
comf="/Data/100WSc.OTUtable.csv"
msif="/Code/StressLevel/100W.MSI.csv"
wd="/Code/AlphaBetaDiversity"

##################
# 1 # load data
##################
library(SpadeR)
source("/Code/handytool.r")
wd=iwd(wd)
comm=t(lazyopen(comf))
msi=lazyopen(msif)

##################
# 2 # iChao1
##################

ichao=t(sapply(1:nrow(comm),
               function(i)
               {
                 message("i=",i)
                 comi=data.frame(s1=comm[i,],stringsAsFactors = FALSE)
                 chaos=SpadeR::ChaoSpecies(data=comi,datatype="abundance",k=10,conf=0.95)
                 chaos$Species_table[5,]
               }))

rownames(ichao)=rownames(comm)
save.file(ichao,prefix="100WSc",filename = "iChao1")

cor.test(ichao[,1],log(msi[,1]))
plot(ichao[,1]~log(msi[,1]))

##################
# 3 # Shannon
##################
shannon=vegan::diversity(comm, index = "shannon")
save.file(data.frame(shannon),prefix="100WSc",filename = "Shannon")

cor.test(shannon,log(msi[,1]))
plot(shannon~log(msi[,1]))

# End