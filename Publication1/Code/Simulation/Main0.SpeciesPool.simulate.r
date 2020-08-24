rm(list=ls())


###########################
# 0 # file/folder paths and package loading
# you may need to change them to the paths on your computer before testing the code.

wd="./Data/SimulateData/SpeciesPool" # the folder saving the input and output files
code.wd="./Code/Simulation" # the folder save simulation functions
tool.code.wd="./Code/Tools" # the folder saving tools.r

source(paste0(tool.code.wd,"/tools.r"))
library(ape)
library(phytools)
library(vegan)
library(iCAMP)
wd=iwd(wd)

###########################################################
# 1 # use Stegen's tree (Stegen et al 2015 Frontiers in Microbiology) <doi:10.3389/fmicb.2015.00370> with 1139 OTUs
# 1.1 # tree
tree.file="Tree_sim_10002.nwk" # data file from Stegen et al 2015
tree=lazyopen(tree.file)
tree$tip.label=paste0("OTU",tree$tip.label) # set species IDs as OTU1 ... OTU1139
drt=iCAMP::tree.droot(tree,nworker = 4)
(drt.max=max(drt$distRoot)) # check max distance to root
tree$edge.length=tree$edge.length/drt.max # correct the distance to root, to range within 1.
tree
pd=iCAMP::pdist.big(tree = tree, wd = wd, output = TRUE, nworker = 4)
tree.save=list(tree=tree,pd=pd)
save(tree.save,file=paste0(wd,"/JS.tree.pd.rda"))
  
# 1.2 # OPEN, optimum environment, the key trait
# 1.2.1 # method 1: Stegen 2015, low phylogenetic signal across tree
opens=list()
enop.old=lazyopen("All_pops_sim_10002.csv") # data file from Stegen et al 2015
enopv.old=enop.old[,2]
names(enopv.old)=paste0("OTU",rownames(enop.old))
hist(enop.old[,2],breaks=100)
(Kvalue.old=phytools::phylosig(tree, x=enopv.old, method="K", test=FALSE, nsim=1000))
opens$JS=list(openv=enopv.old,K=Kvalue.old)

# 1.2.2 # method 2: Brownian, medium phylogenetic signal across tree
source(paste0(code.wd,"/BM.ACDC.r"))
BMt=BM.ACDC(tree = tree,mean.trait = 0.5,sig2 = 0.25^2,
            bounds = c(0,1),g=NULL,rand.num = 1,
            nworker = 4,code.wd = code.wd)

enop=BMt$trait
enopv=enop[,1];names(enopv)=rownames(enop)
treen=BMt$tree.new
hist(enop[,1],breaks = 50)
(Kvalue=phytools::phylosig(tree, x=enopv, method="K", test=FALSE, nsim=1000))
opens$BM=list(openv=enopv,K=Kvalue)

# 1.2.3 # method 3: ACDC model, high phylogenetic signal across tree
source(paste0(code.wd,"/BM.ACDC.r"))
g=2000
BMt=BM.ACDC(tree = tree,mean.trait = 0.5,sig2 = 0.25^2,
            bounds = c(0,1),g=g,rand.num = 1,
            nworker = 4,code.wd = code.wd)

enop=BMt$trait
enopv=enop[,1];names(enopv)=rownames(enop)
treen=BMt$tree.new
hist(enop[,1],breaks = 100)
(Kvalue=phytools::phylosig(tree, x=enopv, method="K", test=FALSE, nsim=1000))
opens$ACDC=list(openv=enopv,K=Kvalue)

save(opens,file=paste0(wd,"/JS.opens.rda"))

# End #
