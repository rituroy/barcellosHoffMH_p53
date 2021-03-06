## USE R 2.4.0 - results different in R 2.2.1
## Run section 2 before sections of higher order


## ---------------------------------

computerFlag="cluster"
computerFlag=""


if (computerFlag=="cluster") {
    setwd("/home/royr/project/barcellosHoffMH/p53")
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"barcellosHoffMH/p53",sep=""))
}

## ----------------------------------------------

library(oligo)

##############################################
##############################################
## Section 1
## Process data

dataset="A0501_A0508"
dataset="A1020"

fName=paste("_",dataset,sep="")

datadir="data/"
phen=read.table(file=paste(datadir,"phen",fName,".txt", sep=""), header=T, sep="\t", quote="", comment.char="", as.is=T)
#phen=phen[1:4,]

fileList=paste(datadir,"cel/",dataset,"/",phen$celFileName,sep="")

if (F) {
    datadir="tmp/"
    fileList=dir(datadir,pattern="CEL")
    phen=data.frame(id=sapply(fileList,function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F),celFileName=fileList,type=c(0,0,1,1),stringsAsFactors=F)
    fileList=paste(datadir,fileList,sep="")
}

affyGeneFS=read.celfiles(filenames=fileList)
phen=phen[match(colnames(affyGeneFS),phen$celFileName),]
sampleNames(affyGeneFS)=phen$id

if (F) {
    R> library(oligoData)
    R> data(affyGeneFS)
    R> affyGeneFS
}

## ---------------------------
fit=fitProbeLevelModel(affyGeneFS)

if (F) {
png(paste("image",fName,"_%03d.png",sep=""),width=3*240, height=2*240)
image(fit)
dev.off()
}

png(paste("nusePlot",fName,".png",sep=""))
res=NUSE(fit,las=2)
dev.off()
summary(apply(res,2,median,na.rm=T))
resN=res

png(paste("rlePlot",fName,".png",sep=""))
res=RLE(fit,las=2)
dev.off()
summary(apply(res,2,median,na.rm=T))
resR=res
## ---------------------------

genePS=rma(affyGeneFS, target="probeset")
geneCore=rma(affyGeneFS, target="core")

if (F) {
    R> featureData(exonPS)=getNetAffx(exonPS, "probeset")
    R> featureData(exonCore)=getNetAffx(exonCore, "transcript")
    R> featureData(exonFull)=getNetAffx(exonFull, "transcript")
    R> featureData(exonExtd)=getNetAffx(exonExtd, "transcript")
    R> featureData(geneCore)=getNetAffx(geneCore, "transcript")

    R> exonCore
    R> featureData(exonCore)
    R> varLabels(featureData(exonCore))
    pData(featureData(exonCore))[1:2, "geneassignment"]
}
featureData(genePS)=getNetAffx(genePS, "probeset")

## ---------------------------
png(paste("boxplot",fName,".png",sep=""))
res=boxplot(genePS,las=2)
dev.off()
summary(res$stat[3,])
resB=res

png(paste("histogram",fName,".png",sep=""))
hist(genePS)
dev.off()

if (F) {
png(paste("MAplot",fName,"_%03d.png",sep=""),width=3*240, height=2*240)
par(mfrow=c(2,3))
MAplot(genePS)
dev.off()
}
## ---------------------------

expr=exprs(genePS)
anno=pData(featureData(genePS))
anno2=anno
names(anno2)[match(c("probesettype"),names(anno2))]="category"
anno2$geneSym=sapply(anno2$geneassignment,function(x) {
    y=strsplit(x," /// ")[[1]]
    strsplit(y," // ")[[1]][2]
},USE.NAMES=F)

save(expr,phen,file="data.RData")
save(anno,file="data.RData")
save(anno2,file="data.RData")

#######################################################
#######################################################

## --------------------------------------

if (F) {
    expr=data.frame(probesetid)
    rm(probesetid)

    anno=read.table(file=paste(dirSrc,"Affymetrix/MoGene1_0st/MoGene-1_0-st-v1.na31.mm9.transcript.txt",sep=""), header=T, sep="\t", quote="", comment.char="", as.is=T)
    table(anno$transcript_cluster_id==anno$probesetid,exclude=NULL)

    table(duplicated(expr$probesetid))
    table(duplicated(anno$probesetid))
    table(expr$probesetid%in%anno$probesetid) ## why so many unannotated probesets. There weren''t so many with na24.hg18 version - check Ingrid''s data analysis
    #TRUE
    #35556
    table(anno$probesetid%in%expr$probesetid) ## No probsets in annotation file which are not in arrays?
    anno[!(anno$probesetid%in%expr$probesetid),]

    id=match(expr$probesetid,anno$probesetid)
    if (any(is.na(id))) {
        tmp=anno[1:sum(is.na(id)),]
        tmpN=as.numeric(rep(NA,nrow(tmp)))
        for (j in 1:ncol(tmp)) {
            if (is.numeric(tmp[,j])) {tmp[,j]=tmpN
            } else {tmp[,j]=""}
        }
        tmp$probesetid=expr$probesetid[is.na(id)]
        tmp=rbind(anno[id[!is.na(id)],],tmp)
        anno=tmp[match(expr$probesetid,tmp$probesetid),]
    } else {
        anno=anno[match(expr$probesetid,anno$probesetid),]
    }

    save(anno,file="anno.RData")
    rm(anno)

    ## --------------------------------------
    anno2=read.table(file=paste(dirSrc,"Affymetrix/MoGene1_0st/anno_MoGene-1_0-st-v1.na31.mm9.transcript.txt",sep=""), header=T, sep="\t", quote="", comment.char="", as.is=T)

    id=match(expr$probesetid,anno2$probesetid)
    if (any(is.na(id))) {
        tmp2=anno2[1:sum(is.na(id)),]
        tmpN=as.numeric(rep(NA,nrow(tmp2)))
        for (j in 1:ncol(tmp2)) {
            if (is.numeric(tmp2[,j])) {tmp2[,j]=tmpN
            } else {tmp2[,j]=""}
        }
        tmp2$probesetid=expr$probesetid[is.na(id)]
        tmp2=rbind(anno2[id[!is.na(id)],],tmp2)
        anno2=tmp2[match(expr$probesetid,tmp2$probesetid),]
    } else {
        anno2=anno2[match(expr$probesetid,anno2$probesetid),]
    }

    save(anno2,file="anno2.RData")

    rm(anno2)

    ## --------------------------------------

    notRequired=function() {
        id=match(anno1$probesetid,anno2$probesetid)
        id1=which(!is.na(id))
        id2=id[id1]
        summary(abs(expr1[id1,]-expr2[id2,]))
    }

}

#######################################################
## Filter out unexpressed probes
#######################################################

#source(paste(dirSrc,"functions/affyFuncs.3.R",sep=""))

datadir=paste("results/",dataset,"",sep="")
datadir=""

load(paste(datadir,"data.RData",sep=""))

## No missing values
probeNA=rep(NA,nrow(expr))
for (i in 1:nrow(expr)) {
    probeNA[i]=mean(is.na(expr[i,]),na.rm=T)
}
table(probeNA,exclude=NULL)
save(probeNA,file=paste("probeNA.RData",sep=""))

for (pThres in c(.1,.25)) {
    thres=quantile(expr,pThres,na.rm=T)
    probeMean=rep(NA,nrow(expr))
    for (i in 1:nrow(expr)) {
        probeMean[i]=mean(expr[i,]<thres,na.rm=T)
    }
    save(probeMean,file=paste("probeMean_",pThres*100,"Perc.RData",sep=""))
}

probeMad=rep(NA,nrow(expr))
for (i in 1:nrow(expr)) {
    probeMad[i]=mad(expr[i,],na.rm=T)
}
save(probeMad,file=paste("probeMad.RData",sep=""))

## ----------------------------

load(file=paste(datadir,"anno2.RData",sep=""))
load(file=paste(datadir,"anno.RData",sep=""))
#anno=anno[,c("probesetid","crosshybtype","number_cross_hyb_probes")]
anno=anno[,c("probesetid","crosshybtype")]

datadir=paste("results/",dataset,"",sep="")
datadir=""
thres=.10
thres=.25
load(file=paste(datadir,"probeMean_",100*thres,"Perc.RData",sep=""))
load(file=paste(datadir,"probeMad.RData",sep=""))
probeMadThres=quantile(probeMad,thres,na.rm=T)

table(probeMean<1,exclude=NULL)
table(probeMad>probeMadThres,exclude=NULL)
table(probeMean<1,probeMad>probeMadThres,exclude=NULL)
table(anno$crosshybtype,exclude=NULL)
table(anno2$assignXhyb_mRNA==1,exclude=NULL)
table(anno$crosshybtype,anno2$assignXhyb_mRNA==1,exclude=NULL)

## USE THIS
#clId_probeMean.RData
clId=probeMean<1 & anno2$category%in%c("main")
save(clId,file=paste("clId.RData",sep=""))

clId=probeMean<1 & probeMad>probeMadThres & anno2$category%in%c("main")

#clId_xhyb.RData
clId=probeMean<1 & anno2$category%in%c("main") & !is.na(anno$crosshybtype) & anno$crosshybtype==1

#clId_assignScore100.RData
clId=probeMean<1 & anno2$category%in%c("main") & !is.na(anno$crosshybtype) & anno$crosshybtype==1 & !is.na(anno2$assignScore_mRNA) & anno2$assignScore_mRNA==100

#clId_probeMad_assignScore100.RData
clId=probeMean<1 & probeMad>probeMadThres & anno2$category%in%c("main") & !is.na(anno$crosshybtype) & anno$crosshybtype==1 & !is.na(anno2$assignScore_mRNA) & anno2$assignScore_mRNA==100


#######################################################

eset=list(expr=expr,ann=anno2,phen=phen)
save(eset,file="eset.RData")

save.image("tmp.RData")



#######################################################
#######################################################
## Section 2
## Load processed data

datadir="results/"
load(paste(datadir,"eset.RData",sep=""))
load(paste(datadir,"clId.RData",sep=""))
phen=eset$phen
datadir="docs/"
phen2=read.table(paste(datadir,"Suervised cluster RNAseq and Microarray data - microarray.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
names(phen2)[match(c("Sample.ID","Fast.Slopes","X","Sample.ID.1","Slow.Slopes"),names(phen2))]=c("id","slope","X","id2","slope2")
phen3=phen2
phen3$id=phen3$id2
phen3$slope=phen3$slope2
phen3$slopeCat="slow"
phen2$slopeCat="fast"
phen2=rbind(phen2,phen3)
phen2$id=gsub("-","_",gsub(" ","",phen2$id))
j=match(phen$id,phen2$id); j1=which(!is.na(j)); j2=j[j1]
phen$slope=NA
phen$slopeCat=""
phen[j1,c("slope","slopeCat")]=phen2[j2,c("slope","slopeCat")]
phen$slopeCat[which(phen$slopeCat=="")]=NA
eset$phen=phen
rm(phen,phen2,phen3)

#######################################################
#######################################################
## Section 3
## Association of heatmap clusters with clinical variables

datadir=""

fileList=dir(,pattern="clusterInfoFeature")
for (fId in 1:length(fileList)) {
    fName=fileList[fId]
    tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}


varList=c("treatment","group","experimentNo")
var1="clustId"
fileList=dir(,pattern="clusterInfoSample")
tmp=rep(NA,length(fileList))
tmpC=rep("",length(fileList))
n=length(fileList)
out=data.frame(fileList,tmp,tmp,tmp,stringsAsFactors=F)
names(out)=c("file",varList)
for (fId in 1:length(fileList)) {
    fName=fileList[fId]
    out$file[fId]=fName
    tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    for (vId in 1:length(varList)) {
        out[fId,varList[vId]]=fisher.test(tbl[,var1],tbl[,varList[vId]])$p.value
    }
}
out



#######################################################
#######################################################
## DELETE


datadir=""
tbl1=read.table(file=paste(datadir,"tmp1.txt", sep=""), header=F, sep="\t", quote="", comment.char="", as.is=T)
tbl2=read.table(file=paste(datadir,"tmp2.txt", sep=""), header=F, sep="\t", quote="", comment.char="", as.is=T)

x1=tbl1[,1]
x2=tbl2[,1]


write.table(paste("qRscript tmp_mmu.R ",x1[!x1%in%x2],sep=""), file="tmp3.txt",col.names=F,row.names=F, sep="\t",quote=F)

/home/royr/project/barcellosHoffMH/p53/htseqCountData/tmp,tophat2,gtf/Mus_musculus
/home/royr/project/barcellosHoffMH/p53/tophat2Data/tmp,tophat2,gtf/Mus_musculus

/cbc/data/fastqData/tmp/Mus_musculus


#######################################################
#######################################################
## Section 4
## DE genes
## Run section 2 first

library(limma)
library(qvalue)
library(sva)
source(paste(dirSrc,"functions/TTest.9.1.6.R",sep=""))

fName1=""

varInfo=data.frame(formula=c("group+experimentNo","group","experimentNo"),variable=c("treatExpt","group","experimentNo"),name=c("_treatAndExptNoAdj","_treatAdj","_exptNoAdj"),stringsAsFactors=F)
adjFlag="_treatAndExptNoAdj"

varInfo=data.frame(formula=c("treatExpt","group","experimentNo"),variable=c("treatExpt","group","experimentNo"),name=c("_treatExptNoAdj","_treatAdj","_exptNoAdj"),stringsAsFactors=F)
adjFlag="_exptNoAdj"
adjFlag="_treatAdj"
adjFlag="_treatExptNoAdj"
adjFlag=""

## Model with treatment adjusted
varInfo=data.frame(formula=c("group"),variable=c("group"),name=c("_treatAdj"),stringsAsFactors=F)
adjFlag=""

## Model with no adjustment
varInfo=NULL
adjFlag=""

grp=as.factor(eset$phen$group)
getLog2Mean=function(x) {
    apply(x,1,mean,na.rm=T)
}
log2FC=as.matrix(cbind(
    gammaVsham=getLog2Mean(eset$expr[clId,which(grp=="_")])-getLog2Mean(eset$expr[clId,which(grp=="Sham")]),
    hzeVsham=getLog2Mean(eset$expr[clId,which(grp=="HZE")])-getLog2Mean(eset$expr[clId,which(grp=="Sham")]),
    hzeVgamma=getLog2Mean(eset$expr[clId,which(grp=="HZE")])-getLog2Mean(eset$expr[clId,which(grp=="_")])
))

grp=as.integer(eset$phen$group!="Sham")
grp[which(eset$phen$group=="HZE")]=2
grp=as.factor(grp)
grp=as.factor(eset$phen$group)
contrasts(grp)=contr.sum(sum(!duplicated(grp)))
#contrasts(grp)=contr.treatment(sum(!duplicated(grp)))
design=model.matrix(~0+grp,ref="Sham")
colnames(design)=sub("(Intercept)","intercept",colnames(design),fixed=T)
fit=lmFit(eset$expr[clId,],design)
contMat=makeContrasts(gammaVsham=grp_-grpSham,hzeVsham=grpHZE-grpSham,hzeVgamma=grpHZE-grp_,levels=design)
#contMat=makeContrasts(gammaVsham=grp1,hzeVsham=grp2,hzeVgamma=grp2-grp1,levels=design)
#contMat=makeContrasts(gammaVsham=grp2,hzeVsham=grp3,hzeVgamma=grp3-grp2,levels=design)
fit2=contrasts.fit(fit, contMat)
fit2=eBayes(fit2)

phen0=eset$phen
phen0$group[which(phen0$group=="Sham")]="0Sham"
phen0$treatExpt=paste(phen0$group,phen0$experimentNo)
phen0$slopeCat3=as.integer(phen0$slopeCat=="fast")
phen0$slopeCat3[!phen0$slopeCat%in%c("slow","fast")]=NA

#grp=as.integer(eset$phen$group!="0Sham")
#grp[which(eset$phen$group=="HZE")]=2
compList=c("gammaVsham","hzeVsham","hzeVgamma")
subsetList=""

compList="slope"
subsetList=c("",paste("_",c("gamma","hze","sham"),sep=""))

compList="fastVsSlow"
subsetList=""
subsetList=c("",paste("_",c("gamma"),sep=""))
subsetList=c("",paste("_",c("gamma","hze","sham"),sep=""))

colId=2
tmp=matrix(nrow=nrow(fit2$coef),ncol=length(compList)*length(subsetList),dimnames=list(rownames(fit2$coef),paste(compList,rep(subsetList,each=length(compList)),sep="")))
fit3=list(coef=tmp,p.value=tmp)
fit4=list(coef=tmp,p.value=tmp)
for (subsetFlag in subsetList) {
    for (compId in 1:length(compList)) {
        compFlag=compList[compId]
        varList2="experimentNo"
        varType="categorical"
        switch(compFlag,
        "gammaVsham"={
            varList="group"; grpUniq=c("0Sham","_"); grpName=c("sham","gamma")
        },
        "hzeVsham"={
            varList="group"; grpUniq=c("0Sham","HZE"); grpName=c("sham","hze")
        },
        "hzeVgamma"={
            varList="group"; grpUniq=c("_","HZE"); grpName=c("gamma","hze")
        },
        "slope"={
            varType="continuous"
            varList="slope"; grpUniq=grpName=NULL
        },
        "fastVsSlow"={
            varList="slopeCat3"; grpUniq=c(0,1); grpName=c("slow","fast")
        }
        )
        samId=1:nrow(phen0)
        if (is.null(grpUniq)) {
            samId=samId[!is.na(phen0[samId,varList])]
        } else {
            samId=samId[which(phen0[samId,varList]%in%grpUniq)]
        }
        switch(subsetFlag,
            "_gamma"={
                samId=samId[which(phen0$group[samId]%in%c("_"))]
            },
            "_hze"={
                samId=samId[which(phen0$group[samId]%in%c("HZE"))]
            },
            "_sham"={
                samId=samId[which(phen0$group[samId]%in%c("0Sham"))]
            }
        )
        phen=phen0[samId,]
        covFlag=""
        if (is.null(varInfo)) {
            varList2=NULL
        } else {
            if (varList=="slope") {
                var2Info=varInfo
                if (adjFlag%in%c("_treatExptNoAdj","_treatAndExptNoAdj") & sum(!duplicated(phen$group))==1) {
                    adj2Flag="_exptNoAdj"
                } else {
                    adj2Flag=adjFlag
                }
                k=which(var2Info$name==adj2Flag)
                if (length(k)==1 && k>1) var2Info=var2Info[-(1:(length(k)-1)),]
                covFlag=""
                varList2=NULL
                for (k in 1:nrow(var2Info)) {
                    x=table(phen[,var2Info$variable[k]])
                    x=x[which(x!=1)]
                    if (length(x)>1) {
                        varList2=strsplit(var2Info$formula[k],"+",fixed=T)[[1]]
                        covFlag=var2Info$name[k]
                        break
                    }
                }
                if (F) {
                    varList2="group"
                    x=table(phen[,varList2])
                    x=x[which(x!=1)]
                    if (adj2Flag!="_treatAdj") x=1
                    if (length(x)>1) {
                        covFlag="_treatAdj"
                    } else {
                        varList2="experimentNo"
                        x=table(phen[,varList2])
                        x=x[which(x!=1)]
                        if (adj2Flag!="_exptNoAdj") x=1
                        if (length(x)>1) {
                            covFlag="_exptNoAdj"
                        } else {
                            covFlag=""
                            varList2=NULL
                        }
                    }
                }
                if (F) {
                    if (covFlag!="") {
                        for (varThis in varList2) {
                            x=table(phen[,varThis])
                            if (any(x==1)) {
                                cat("\nExcluded singular groups ",names(x)[x==1],"\n",sep="")
                                x=x[which(x!=1)]
                                samId=samId[which(phen[,varThis]%in%names(x))]
                                phen=phen0[samId,]
                            }
                        }
                    }
                }
            } else {
                var2Info=varInfo
                if (adjFlag%in%c("_treatExptNoAdj","_treatAndExptNoAdj") & sum(!duplicated(phen$group))==1) {
                    adj2Flag="_exptNoAdj"
                } else {
                    adj2Flag=adjFlag
                }
                k=which(var2Info$name==adj2Flag)
                if (length(k)==1 && k>1) var2Info=var2Info[-(1:(length(k)-1)),]
                covFlag=""
                varList2=NULL
                for (k in 1:nrow(var2Info)) {
                    x=table(phen[,varList],phen[,var2Info$variable[k]])
                    y=x!=0
                    if (nrow(x)>1 && ncol(x)>1 && any(y[1,]==y[2,])) {
                        varList2=strsplit(var2Info$formula[k],"+",fixed=T)[[1]]
                        covFlag=var2Info$name[k]
                        break
                    }
                }
            }
        }
        cat("\n\n===================== ",compFlag,covFlag,subsetFlag,"\n",sep=": ")
        cat("No. of samples: ",length(samId),"\n",sep="")
        if (covFlag=="") {
            modelThis=paste("~",varList,sep="")
        } else {
            modelThis=paste("~",varList,paste("+as.factor(",varList2,")",collapse=""),sep="")
        }
        if (varType=="categorical") {
            modelThis=sub("varList","as.factor(varList)",modelThis)
        }
        modelThis=as.formula(modelThis)
        design=model.matrix(modelThis,data=phen)
        fit=lmFit(eset$expr[clId,samId],design)
        fit=eBayes(fit)
        fit3$coef[,paste(compFlag,subsetFlag,sep="")]=fit$coef[,colId]
        fit3$p.value[,paste(compFlag,subsetFlag,sep="")]=fit$p.value[,colId]
        out=data.frame(probesetid=rownames(fit$coef),coef=fit$coef[,colId],t=fit$t[,colId],pv=fit$p.value[,colId])
        write.table(out,file=paste("stat_",compFlag,covFlag,subsetFlag,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
        
        mod=design
        if (covFlag=="") {
            mod0=matrix(rep(1,nrow(phen)),ncol=1); colnames(mod0)="(Intercept)"
        } else {
            mod0=model.matrix(as.formula(paste("~",paste("as.factor(",varList2,")",collapse="+"),sep="")), data=phen)
        }
        #cat(colnames(mod),"\n")
        #cat(colnames(mod0),"\n")
        svObj=try(sva(eset$expr[clId,samId],mod,mod0))
        nm=NULL
        if (class(svObj)!="try-error") {
            if (is.matrix(svObj$sv)) {
                nm=c(colnames(mod),paste("sv",1:ncol(svObj$sv),sep=""))
            } else if (is.numeric(svObj$sv)) {
                if (length(svObj$sv)>1) {
                    nm=c(colnames(mod),"sv1")
                } else {
                    nm=NULL
                }
            }
        }
        if (!is.null(nm)) {
            #dat=svObj$sv; colnames(dat)=paste("sv",1:ncol(dat),sep="")
            design=cbind(mod,svObj$sv)
            colnames(design)=nm
            fit=lmFit(eset$expr[clId,samId],design)
            fit=eBayes(fit)
            fit4$coef[,paste(compFlag,subsetFlag,sep="")]=fit$coef[,colId]
            fit4$p.value[,paste(compFlag,subsetFlag,sep="")]=fit$p.value[,colId]
            out=data.frame(probesetid=rownames(fit$coef),coef=fit$coef[,colId],t=fit$t[,colId],pv=fit$p.value[,colId])
            write.table(out,file=paste("stat_",compFlag,covFlag,subsetFlag,"_sva",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
        } else {
            cat("Could not run SVA !!!\n",sep="")
        }
    }
}

## ----------------------------------------------
fitThis=fit3
fitThis=fit2
png("tmp.png")
par(mfrow=c(2,2))
for (k in 1:ncol(fitThis$coef)) {
    lim=range(log2FC[,k],fitThis$coef[,k],na.rm=T)
    plot(log2FC[,k],fitThis$coef[,k],xlim=lim,ylim=lim,main=colnames(fitThis$coef)[k],xlab="Observed log2FC",ylab="Estimated log2FC")
    abline(c(0,1))
}
dev.off()

colIdPV="pv"
pThres=0.05
for (k in 1:ncol(fitThis$coef)) {
    if (k==0) {
        stat2=data.frame(F=fitThis$F,pv=fitThis$F.p.value,qv=rep(NA,nrow(fitThis$coef)))
        heading="ANOVA"
    } else {
        stat2=data.frame(pv=fitThis$p.value[,k],qv=rep(NA,nrow(fitThis$coef)))
        heading=colnames(fitThis$coef)[k]
    }
    cat("\n\n===========",heading,"=====\n")
    i=which(!is.na(stat2$pv))
    #stat2$qv[i]=getAdjustedpv(stat2[i,colIdPV],method="qv",strict=T)
    stat2$qv[i]=qvalue(stat2$pv[i])$qvalues
    print(table(stat2$qv<pThres))
    print(summary(stat2$pv))
}

colId=c("probesetid","seqname","start","stop","geneSym")
pThres2=0.05
k=which(colnames(fitThis$coef)=="hzeVgamma")
stat2=data.frame(pv=fitThis$p.value[,k],qv=rep(NA,nrow(fitThis$coef)))
heading=colnames(fitThis$coef)[k]
table(gammaVshamPV0.05=fitThis$p.value[,"gammaVsham"]<pThres2,hzeVshamPV0.05=fitThis$p.value[,"hzeVsham"]<pThres2)
"
                hzeVshamPV0.05
gammaVshamPV0.05  FALSE   TRUE
            FALSE 219663   5535
            TRUE    9187   2061
"
i=which(fitThis$p.value[,"gammaVsham"]<pThres2 & fitThis$p.value[,"hzeVsham"]<pThres2 & !is.na(stat2$pv))
i=which(fitThis$p.value[,"gammaVsham"]>=pThres2 & fitThis$p.value[,"hzeVsham"]>=pThres2 & !is.na(stat2$pv))
i=which(fitThis$p.value[,"gammaVsham"]<pThres2 & fitThis$p.value[,"hzeVsham"]>=pThres2 & !is.na(stat2$pv))
i=which(fitThis$p.value[,"gammaVsham"]>=pThres2 & fitThis$p.value[,"hzeVsham"]<pThres2 & !is.na(stat2$pv))
i=which(fitThis$p.value[,"gammaVsham"]<pThres2)
i=which(fitThis$p.value[,"hzeVsham"]<pThres2)
stat2$qv=NA
res=try(qvalue(stat2$pv[i]))
if (class(res)=="try-error") {
    stat2$qv[i]=p.adjust(stat2$pv[i],method="BH")
} else {
    stat2$qv[i]=res$qvalues
}
print(table(stat2$qv<pThres))
i=which(stat2$qv<pThres)
eset$ann[clId,colId][i,]
"
probesetid seqname    start     stop geneSym
17385631   17385631    chr2 60444171 60444195  Pla2r1
17385636   17385636    chr2 60452429 60452453  Pla2r1
"
i=which(fitThis$p.value[,"gammaVsham"]>=pThres2 & fitThis$p.value[,"hzeVsham"]>=pThres2 & !is.na(stat2$pv))
stat2$qv=NA
stat2$qv[i]=qvalue(stat2$pv[i])$qvalues
print(table(stat2$qv<pThres))
i=which(stat2$qv<pThres)
eset$ann[clId,colId][i,]

i=which(fitThis$p.value[,"gammaVsham"]<pThres2 & fitThis$p.value[,"hzeVsham"]<pThres2 & !is.na(stat2$pv))
stat2$qv=NA
stat2$qv[i]=qvalue(stat2$pv[i])$qvalues
ii=which(stat2$qv<pThres)
png(paste("scatterPlots","_",eset$ann$geneSym[clId][ii[1]],fName1,".png",sep=""),width=4*240,height=2*240)
par(mfcol=c(1,2))
lim=range(c(eset$expr[clId,][ii,]),na.rm=T)
for (i2 in ii) {
    i1=which(clId)[i2]
    boxplot(eset$expr[i1,]~grp,names=c("Sham","Gamma","HZE"),ylim=lim,main=paste(eset$ann$probesetid[i1],", ",eset$ann$geneSym[i1],sep=""))
    for (k in 1:sum(!duplicated(grp))) {
        j=which(grp==(k-1))
        points(rep(k,length(j)),eset$expr[i1,j])
        lines(k+c(-.5,.5),rep(mean(eset$expr[i1,j],na.rm=T),2),col="green")
    }
}
legend(3,lim[2]-.1,"mean",col="green",lty="solid")
dev.off()

png(paste("densityPlot_slope.png"))
plot(density(phen0$slope,na.rm=T),main="slope")
dev.off()

## ----------------------------------------------

if (F) {
    # fit30 & fit40 - models with no adjustment
    # fit31 & fit41 - models with treatment adjustment
    fit3=fit30
    fit3$coef=cbind(fit31$coef[,1],fit30$coef)
    fit3$p.value=cbind(fit31$p.value[,1],fit30$p.value)
    colnames(fit3$coef)=colnames(fit3$p.value)=c("fastVsSlow_treatAdj",colnames(fit30$coef))
    fit4=fit40
    fit4$coef=cbind(fit41$coef[,1],fit40$coef)
    fit4$p.value=cbind(fit41$p.value[,1],fit40$p.value)
    colnames(fit4$coef)=colnames(fit4$p.value)=c("fastVsSlow_treatAdj",colnames(fit40$coef))
}

## ----------------------------------------------
fName1="_moGene2.0"
fName1="_sva_moGene2.0"

colIdPV="pv"; colNamePV="PV"; pThres=10^-6
colIdPV="qv"; colNamePV="QV"; pThres=0.05

pThres2=0.05

onePlotFlag=F
onePlotFlag=T

cexThis=.8
cexThis=1.2

compFlag3=""

compList1=c("gammaVsham","hzeVsham","hzeVgamma")
compList1=c("gammaVsham","hzeVsham")
subsetList=""

compList1="slope"
subsetList=c("",paste("_",c("gamma","hze","sham"),sep=""))

compList1="fastVsSlow"
subsetList=c("",paste("_",c("gamma","hze","sham"),sep=""))
subsetList=c("_treatAdj","",paste("_",c("gamma","hze","sham"),sep=""))

fName2=ifelse(length(compList1)==1,paste("_",compList1,sep=""),subsetList)
if (fName2=="_") fName2=""

switch(fName1,
    "_moGene2.0"={
        fitThis=fit3
    },
    "_sva_moGene2.0"={
        fitThis=fit4
    }
)

out1=NULL
if (compFlag3=="") {
    if (onePlotFlag) {
        if (compFlag3=="") {
            #png(paste("dePlots",fName2,fName1,".png",sep=""),width=6*240,height=3*240)
            png(paste("dePlots",fName2,fName1,".png",sep=""),width=length(compList1)*length(subsetList)*2*240,height=3*240)
            par(mfcol=c(2,length(compList1)*length(subsetList)))
        } else {
            png(paste("dePlots",fName2,fName1,".png",sep=""),width=6*240,height=2*240)
            par(mfcol=c(2,6))
        }
    }
}
for (subsetFlag in subsetList) {
    for (compFlag in compList1) {
        if (compFlag3=="") compList=compFlag
        if (compFlag3!="") {
            if (onePlotFlag) {
                if (compFlag3=="") {
                    png(paste("dePlots",compFlag3,"_",compFlag,subsetFlag,fName1,".png",sep=""),width=4*240,height=2*240)
                    par(mfcol=c(2,4))
                } else {
                    png(paste("dePlots",compFlag3,"_",compFlag,subsetFlag,fName1,".png",sep=""),width=length(compList)*240,height=2*240)
                    par(mfcol=c(2,length(compList)))
                }
            }
        }
        subsetName="All samples: "
        switch(subsetFlag,
            "_gamma"={
                subsetName="Gamma: "
            },
            "_hze"={
                subsetName="HZE: "
            },
            "_sham"={
                subsetName="Sham: "
            }
        )
        for (compThis in compList) {
            compFlag2=compFlag
            colGeneId="probesetid"
            stat_1=data.frame(probesetid=rownames(fitThis$coef),coef=fitThis$coef[,paste(compFlag,subsetFlag,sep="")],pv=fitThis$p.value[,paste(compFlag,subsetFlag,sep="")],qv=rep(NA,nrow(fitThis$coef)))
            i=which(!is.na(stat_1$pv))
            if (length(i)==0) next
            stat_1$qv[i]=qvalue(stat_1$pv[i])$qvalues
            ann_1=eset$ann[clId,]
            switch(compFlag,
                "gammaVsham"={
                    compName="Gamma vs. sham"
                },
                "hzeVsham"={
                    compName="HZE vs. sham"
                },
                "hzeVgamma"={
                    compName="HZE vs. gamma"
                },
                "slope"={
                    compName="Slope"
                },
                "fastVsSlow"={
                    compName="Fast vs. slow"
                }
            )
            switch(subsetFlag,
                "_treatAdj"={
                    compName=paste(compName,", treat adjusted",sep="")
                }
            )
            if (compFlag3!="") compName=paste(compThis,", ",compName,sep="")
            i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
            ann_1=ann_1[i1,]
            cat("\n\n",subsetName,compName,"\n",sep="")
            x=table(stat_1$coef>0,stat_1[,colIdPV]<pThres)
            rownames(x)=c("Down","Up")[match(rownames(x),c("FALSE","TRUE"))]
            colnames(x)=paste(colNamePV,c(">=","<"),pThres,sep="")[match(colnames(x),c("FALSE","TRUE"))]
            print(x)
            
            cexMain=cexLab=1.5
            cexMain=cexLab=2.5
            cexMain=2.5; cexLab=1.5
            
            if (!onePlotFlag) png(paste("histogram_",compFlag,subsetFlag,".png",sep=""))
            hist(stat_1$pv,xlab="P-value",main=paste(subsetName,compName,sep=""),pch=19,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
            if (!onePlotFlag) dev.off()
            
            if (!onePlotFlag) png(paste("volcanoPlot_",compFlag,subsetFlag,"_",tolower(colNamePV),pThres,".png",sep=""))
            i=which(stat_1[,colIdPV]<pThres)
            plot(stat_1$coef,-log10(stat_1$pv),main=paste(subsetName,compName,"\nNo. with ",tolower(colNamePV)," < ",pThres,": ",length(i),sep=""),xlab="Log2 fold change",ylab="-Log10(p-value)",pch=19,cex=.8,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
            if (length(i)!=0) points(stat_1$coef[i],-log10(stat_1$pv[i]),col="red",pch=19,cex=cexThis)
            if (!onePlotFlag) dev.off()

            out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("coef","pv","qv")]
            rownames(out2)=ann_1[,colGeneId]
            chr=sapply(ann_1$seqname,function(x) {
                y=strsplit(x,"_")[[1]]
                chr=sub("chr","",y[1])
                chr1=chr
                if (length(chr)!=1) print(x)
                switch(chr1,
                "X"={chr="23"},
                "Y"={chr="24"},
                "M"={chr="99"},
                "Un"={chr="100"}
                )
                chr=as.integer(chr)
                if (length(y)!=1) chr=1000+chr
                chr
            },USE.NAMES=F)
            i=1:nrow(out2)
            i=order(chr,ann_1$start)
            i=order(out2$pv)
            out2=out2[i,]
            names(out2)=paste(names(out2),"_",compFlag2,subsetFlag,sep="")
            names(out2)=sub("coef_","log2FC_",names(out2))
            if (is.null(out1)) {
                out1=ann_1
            }
            out1=cbind(out1,out2[match(out1[,colGeneId],rownames(out2)),])
            out2=cbind(ann_1[i,which(!names(ann_1)%in%colGeneId)],out2)
            #write.table(out2, file=paste("stat_",compFlag2,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)

            if (F) {
                out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("coef","pv","qv")]
                i=order(out2$pv); i=i[which(out2$pv[i]<pThres2)]
                out2=out2[i,]
                names(out2)=paste(names(out2),"_",compFlag2,sep="")
                out2=cbind(ann_1[i,which(!names(ann_1)%in%colGeneId)],out2)
                write.table(out2, file=paste("stat_",compFlag2,"_pv",pThres2,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
            }

            if (F) {
                write.table(unique(ann_1$geneId), file=paste("ensGeneId_",compFlag,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
                i=which(stat_1$pv<pThres)
                if (length(i)!=0) write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_pv",pThres,"_",compFlag,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
                i=order(stat_1$pv)[1:200]
                write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_top200_",compFlag,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
            }
        }
        if (onePlotFlag & compFlag3!="") dev.off()
    }
}
write.table(out1, file=paste("stat",fName2,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
if (onePlotFlag & compFlag3=="") dev.off()

"
Gamma vs. sham

QV>=0.05
Down   118560
Up     117886


HZE vs. sham

QV>=0.05
Down   119872
Up     116574


HZE vs. gamma

QV>=0.05
Down   118942
Up     117504
"

## Check distribution of top genes
png(paste("scatterPlot_topGenes",fName2,fName1,"_%1d.png",sep=""))
varList="slope"
samId=1:nrow(phen0)
samId=samId[!is.na(phen0[samId,varList])]
for (i in which(stat_1[,colIdPV]<pThres)) {
    plot(phen0$slope[samId],eset$expr[clId,samId][i,],main=paste(ann_1$probesetid[i]," ",ann_1$geneSym[i],": pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2 expression")
}
dev.off()


## ----------------------------------------------

##############################################
library(sva)

load("results/eset.RData")

modCom=model.matrix(~as.factor(group), data=eset$phen)

exprCom=ComBat(dat=eset$expr, batch=eset$phen$experimentNo, mod=modCom, par.prior=TRUE,prior.plots=FALSE)
#exprComNP=ComBat(dat=eset$expr, batch=eset$phen$experimentNo, mod=modCom, par.prior=FALSE,prior.plots=FALSE)

mod=model.matrix(~as.factor(group)+as.factor(experimentNo), data=eset$phen)
mod0=model.matrix(~as.factor(experimentNo), data=eset$phen)


numSV=num.sv(eset$expr,mod=mod,method="leek")
numSV
svObj=sva(eset$expr,mod,mod0,n.sv=numSV)

svObj0=sva(eset$expr,mod,mod0)


save(exprCom,exprComNP,svObj,svObj0,file="sva.RData")

library(coin)

colId=c("group","experimentNo")
pvMat=matrix(nrow=ncol(svObj0$sv),ncol=length(colId))
colnames(pvMat)=colId
for (k1 in colId) {
    for (k2 in 1:ncol(svObj0$sv)) {
        pvMat[k2,k1]=pvalue(kruskal_test(svObj0$sv[,k2]~as.factor(eset$phen[,k1])))
    }
}
# None of the surrogate variables are associated with phenotypic variables
pvMat

##############################################
## NOT USED

chr=sub("chr","",ann_1$seqname)
chr[which(chr=="X")]="22"
chr[which(chr=="Y")]="23"
chr=as.integer(chr)
x=paste(chr*10^9+ann_1$start)
x=paste(chr*10^9+ann_1$start)
summary(order(x))

##############################################
##############################################
