# handy functions
#####################
# 1 # open file
#####################
open.file<-function(filename,folder=NA,header=TRUE,row.names=1,load.obj=1,stringsAsFactors=FALSE,check.names=FALSE,...)
{
  # version 20160730 by Daliang Ning (ningdaliang@gmail.com)
  # version 20171104 by Daliang Ning (ningdaliang@gmail.com)
  # version 20200612 by Daliang Ning (ningdaliang@gmail.com), add tsv
  if(is.na(folder)){folder=getwd()}
  suffix=regexpr("\\.[^\\.]*$", filename)
  filetype=substr(filename,suffix,nchar(filename))
  if(filetype %in% c(".nwk",".tre",".tree"))
  {
    #library(ape)
    file=ape::read.tree(file=paste(folder,"/",filename,sep = ""),...)
  }else if(filetype==".fna"|filetype==".fasta"){
    #library(seqinr)
    file=seqinr::read.fasta(file = paste(folder,"/",filename,sep = ""),seqtype = "DNA",as.string = FALSE,forceDNAtolower = FALSE,set.attributes = TRUE,legacy.mode = TRUE,seqonly = FALSE,strip.desc = FALSE,bfa = FALSE)
  }else if(filetype==".RData"){
    if(load.obj=="all")
    {
      old.ls=ls(envir = .GlobalEnv)
      sp=load(paste(folder,"/",filename,sep = ""),envir = .GlobalEnv)
      new.ls=sp
      message("All objects are loaded, including ",sp)
      message("Replaced objects include ",intersect(old.ls,new.ls))
      file=sp
    }else{
      sp=load(paste(folder,"/",filename,sep = ""))
      message("Only one object is loaded: ",sp[load.obj])
      file=get(sp[load.obj])
    }
  }else if(filetype==".rda")
  {
    sp=load(paste(folder,"/",filename,sep = ""))
    file=get(sp)
  }else{
    if(filetype %in% c(".txt",".tsv")){delim="\t"}else if(filetype==".csv"){delim=","}else if(filetype==".tabular"){delim=""}else{message("The file type of ",filename," is not included in this function.")}
    file=read.table(file = paste(folder,"/",filename,sep = ""),header=header,sep = delim,
                    row.names = row.names,stringsAsFactors=stringsAsFactors,check.names = check.names,
                    comment.char = "",...)
  }
  file
}

lazyopen<-function(filepath,load.obj = "all",...)
{
  if(is.null(filepath)){file=NULL}else{
    filepath=gsub("\\\\","/",filepath)
    if(regexpr("/",filepath)>0)
    {
      filepath=sub("file:///","",filepath)
      locate=as.vector(gregexpr("/",filepath)[[1]]);locate=locate[length(locate)]
      file=open.file(filename=substring(filepath,locate+1),folder = substr(filepath,1,locate-1),load.obj = load.obj,...)
    }else{
      file=open.file(filename=filepath,folder = getwd(),load.obj = load.obj,...)
    }
  }
  file
}
#####################
# 2 # save file
#####################
save.file<-function(obj,prefix=numeric(0),filename=NA,folder=NA,
                    type=c("csv","txt","RData","nwk","fasta"),
                    rowname.title="ID", replace=TRUE,row.names=TRUE,...)
{
  # version 20171109
  if(is.null(obj)){message(deparse(substitute(obj))," is NULL, nothing to save.");return(invisible())}
  if(is.na(folder)){folder=getwd()}
  if(is.na(filename)){filename=as.character(substitute(obj))}
  if(length(prefix)!=0){prefix=paste(prefix,".",sep = "")}
  tcode=format(Sys.time(),format = "%y%m%H%M%S")
  if(type[1]=="csv")
  {
    if(row.names)
    {
      obj=data.frame(obj,stringsAsFactors = FALSE)
      if(is.null(rownames(obj))) rownames(obj)=1:nrow(obj)
      obj=data.frame(rownames(obj),obj)
      colnames(obj)[1]=rowname.title
    }
    if(file.exists(paste(folder,"/",prefix,filename,".csv",sep = ""))&(!replace))
    {
      write.csv(obj,file = paste(folder,"/",prefix,filename,".",tcode,".csv",sep = ""),
                row.names=FALSE,...)
    }else{
      write.csv(obj,file = paste(folder,"/",prefix,filename,".csv",sep = ""),
                row.names=FALSE,...)
    }
  }else if(type[1]=="RData"){
    if(file.exists(paste(folder,"/",prefix,filename,".RData",sep = ""))&(!replace))
    {
      save(obj,file = paste(folder,"/",prefix,filename,".",tcode,".RData",sep = ""))
    }else{
      save(obj,file = paste(folder,"/",prefix,filename,".RData",sep = ""))
    }
  }else if(type[1]=="txt"){
    if(row.names)
    {
      obj=data.frame(obj,stringsAsFactors = FALSE)
      if(is.null(rownames(obj))) rownames(obj)=1:nrow(obj)
      obj=data.frame(rownames(obj),obj)
      colnames(obj)[1]=rowname.title
    }
    if(file.exists(paste(folder,"/",prefix,filename,".txt",sep = ""))&(!replace))
    {
      write.table(obj,file=paste(folder,"/",prefix,filename,".",tcode,".txt",sep = ""),
                  quote=FALSE,sep="\t",row.names=FALSE,...)
    }else{
      write.table(obj,file=paste(folder,"/",prefix,filename,".txt",sep = ""),
                  quote=FALSE,sep="\t",row.names=FALSE,...)
    }
  }else if(type[1]=="nwk"){
    #library(ape)
    if(file.exists(paste(folder,"/",prefix,filename,".nwk",sep = ""))&(!replace))
    {
      ape::write.tree(obj,file= paste(folder,"/",prefix,filename,".",tcode,".nwk",sep = ""))
    }else{
      ape::write.tree(obj,file= paste(folder,"/",prefix,filename,".nwk",sep = ""))
    }
  }else if(type[1]=="fasta"){
    #library(seqinr)
    if(file.exists(paste(folder,"/",prefix,filename,".fasta",sep = ""))&(!replace))
    {
      seqinr::write.fasta(obj,names = names(obj),file.out = paste(folder,"/",prefix,filename,".",tcode,".fasta",sep = ""))
    }else{
      seqinr::write.fasta(obj,names = names(obj),file.out = paste(folder,"/",prefix,filename,".fasta",sep = ""))
    }
  }else{
    warning("The file ",filename," was not saved, because type is not included in this function.")
  }
  invisible()
}

#####################
# 3 # set working directory
#####################
iwd<-function(wd=getwd())
{
  wd=sub("file:///","",wd)
  wd=gsub("\\\\","/",wd)
  if(length(list.dirs(wd))==0)
  {
    wd=substr(wd,1,regexpr("\\/[^\\/]*$",wd)-1)
  }
  setwd(wd)
  message("set wd as ",getwd(),"\n",date())
  getwd()
}

#####################
# 4 # match sample or taxa IDs in different tables
#####################
match.name<-function(name.check=integer(0),rn.list=list(integer(0)),cn.list=list(integer(0)),both.list=list(integer(0)),
                     v.list=list(integer(0)),lf.list=list(integer(0)),tree.list=list(integer(0)),group=integer(0),rerank=TRUE,silent=FALSE)
{
  if(length(name.check)==0)
  {
    if(length(tree.list[[1]])!=0){tname1=tree.list[[1]]$tip.label}else{tname1=integer(0)}
    name.check=unique(c(rownames(rn.list[[1]]),colnames(cn.list[[1]]),colnames(both.list[[1]]),names(v.list[[1]]),names(lf.list[[1]]),tname1))
    #union(rownames(rn.list[[1]]),union(colnames(cn.list[[1]]),colnames(both.list[[1]])))
  }
  
  namelist<-function(aa,prefix=integer(0))
  {
    if(length(names(aa))==0)
    {
      names(aa)=paste(prefix,1:length(aa),sep=".")
    }else{
      for(i in 1:length(aa))
      {
        if(nchar(names(aa)[i])==0){names(aa)[i]=paste(prefix,i,sep=".")}
      }
    }
    aa
  }
  
  if(length(rn.list[[1]])!=0)
  {
    rn.list=namelist(rn.list,prefix="rn.list")
    rname.list=lapply(1:length(rn.list), function(i,rl){rownames(rl[[i]])},rl=rn.list)
    names(rname.list)=paste(names(rn.list),".rowname",sep = "")
  }else{
    rname.list<-rn.list<-integer(0)
  }
  
  if(length(cn.list[[1]])!=0)
  {
    cn.list=namelist(cn.list,prefix="cn.list")
    cname.list=lapply(1:length(cn.list), function(i,rl){colnames(rl[[i]])},rl=cn.list)
    names(cname.list)=paste(names(cn.list),".colname",sep = "")
  }else{
    cname.list<-cn.list<-integer(0)
  }
  
  if(length(both.list[[1]])!=0)
  {
    both.list=namelist(both.list,prefix="both.list")
    bname.r.list=lapply(1:length(both.list), function(i,rl){rownames(rl[[i]])},rl=both.list)
    names(bname.r.list)=paste(names(both.list),".rowname",sep = "")
    bname.c.list=lapply(1:length(both.list), function(i,rl){colnames(rl[[i]])},rl=both.list)
    names(bname.c.list)=paste(names(both.list),".colname",sep = "")
    bname.list=c(bname.r.list,bname.c.list)
  }else{
    bname.list<-both.list<-integer(0)
  }
  
  if(length(v.list[[1]])!=0)
  {
    v.list=namelist(v.list,prefix="v.list")
    vname.list=lapply(1:length(v.list),function(i,rl){names(rl[[i]])},rl=v.list)
    names(vname.list)=paste0(names(v.list),".vname")
  }else{
    vname.list<-v.list<-integer(0)
  }
  
  if(length(lf.list[[1]])!=0)
  {
    lf.list=namelist(lf.list,prefix="lf.list")
    fname.list=lapply(1:length(lf.list),function(i,rl){names(rl[[i]])},rl=lf.list)
    names(fname.list)=paste0(names(lf.list),".name")
  }else{
    fname.list<-lf.list<-integer(0)
  }
  
  if(length(tree.list[[1]])!=0)
  {
    tree.list=namelist(tree.list,prefix = "tree.list")
    tname.list=lapply(1:length(tree.list),function(i,rl){rl[[i]]$tip.label},rl=tree.list)
    names(tname.list)=paste0(names(tree.list),".tipname")
  }else{
    tname.list<-tree.list<-integer(0)
  }
  
  if(length(group)!=0)
  {
    if(is.data.frame(group)|is.matrix(group))
    {
      if(ncol(group)>1){warning("Group has more than one column, only the first column was used. ",date())}
      group=as.matrix(group)[,1,drop=FALSE]
      if(is.null(rownames(group)))
      {
        if(nrow(group)!=length(name.check)){stop("Group does not have the same length as expected. ",date())}else{
          warning("Group has no names. The checked name list was assigned to Group.")
          rownames(group)=name.check
        }
      }
      gname.list=list(group.name=rownames(group))
    }else if(is.vector(group)|is.factor(group)){
      if(is.null(names(group)))
      {
        if(length(group)!=length(name.check)){stop("Group does not have the same length as expected. ",date())}else{
          warning("Group has no names. The checked name list was assigned to Group.")
          group=matrix(group,nc=1)
          rownames(group)=name.check
        }
      }else{
        gname.list=list(group.name=names(group))
        grp=matrix(as.vector(group),nc=1)
        rownames(grp)=names(group)
        group=grp
      }
    }else{stop("The format of Group can only be vector, factor, data.frame, or matrix. others are not acceptable. ", date)}
  }else{gname.list=integer(0)}
  
  lname.list=c(rname.list,cname.list,bname.list,vname.list,fname.list,tname.list,gname.list)
  lntn.list=c(rname.list,cname.list,bname.list,vname.list,fname.list,gname.list)
  if(length(lntn.list)>0)
  {
    tot.mis=sum(sapply(lntn.list, function(oldname,name.check){if(length(oldname)==length(name.check)){sum(oldname!=name.check)}else{1}},name.check=name.check))
  }
  if(length(tname.list)>0)
  {
    tot.mis=tot.mis+sum(sapply(tname.list, function(oldname,name.check){sum(!(name.check %in% oldname))+sum(!(oldname %in% name.check))},name.check=name.check))
  }
  
  if(tot.mis==0)
  {
    if(!silent) message("All match very well.")
    output=c(rn.list,cn.list,both.list,v.list,lf.list,tree.list,group)
  }else{
    check.name<-function(i,lname.list,name.check)
    {
      cname=lname.list[[i]]
      los.id=(!(name.check %in% cname))
      add.id=(!(cname %in% name.check))
      mis=sum(add.id)+sum(los.id)
      if(mis>0)
      {
        if(!silent) message("Mismatch warning: ",names(lname.list)[i]," has ",mis," mismatched names. ")
        if(!silent) if(sum(los.id)>0) message("lost name: ", paste(name.check[los.id]," ",sep = ""))
        if(!silent) if(sum(add.id)>0) message("unexpected name: ",paste(cname[add.id]," ",sep = ""))
      }
      mis
    }
    checked=sapply(1:length(lname.list), check.name,lname.list=lname.list,name.check=name.check)
    output=checked
    if(rerank)
    {
      if(sum(checked)>0)
      {
        lname.list=c(list(name.check),lname.list)
        iname=Reduce(intersect,lname.list)
      }else{iname=name.check}
      
      if(is.list(rn.list))
      {
        rn.list.new=lapply(rn.list, function(old,iname){old[match(iname,rownames(old)),,drop=FALSE]},iname=iname)
        names(rn.list.new)=names(rn.list)
      }else{rn.list.new=integer(0)}
      
      if(is.list(cn.list))
      {
        cn.list.new=lapply(cn.list, function(old,iname){old[,match(iname,colnames(old)),drop=FALSE]},iname=iname)
        names(cn.list.new)=names(cn.list)
      }else{cn.list.new=integer(0)}
      
      if(is.list(both.list))
      {
        both.list.new=lapply(both.list, function(old,iname){old[match(iname,rownames(old)),match(iname,colnames(old)),drop=FALSE]},iname=iname)
        names(both.list.new)=names(both.list)
      }else(both.list.new=integer(0))
      
      if(is.list(v.list))
      {
        v.list.new=lapply(v.list,function(old,iname){old[match(iname,names(old))]},iname=iname)
        names(v.list.new)=names(v.list)
      }else{v.list.new=integer(0)}
      
      if(is.list(lf.list))
      {
        f.list.new=lapply(lf.list,function(old,iname){old[match(iname,names(old))]},iname=iname)
        names(f.list.new)=names(lf.list)
      }else{f.list.new=integer(0)}
      
      if(is.list(tree.list))
      {
        tree.list.new=lapply(tree.list,function(tr,iname){rm.tip=tr$tip.label[!(tr$tip.label %in% iname)];ape::drop.tip(tr,rm.tip)},iname=iname)
        names(tree.list.new)=names(tree.list)
      }else{tree.list.new=integer(0)}
      
      if(length(group)!=0)
      {
        group.new=list(group=group[match(iname,rownames(group)),,drop=FALSE])
      }else{group.new=integer(0)}
      
      output=c(rn.list.new,cn.list.new,both.list.new,v.list.new,f.list.new,tree.list.new,group.new)
      if(!silent) message("The names are re-ranked.")
    }
  }
  output
}

