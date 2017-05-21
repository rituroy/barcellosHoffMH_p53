## USE R 2.4.0 - results different in R 2.2.1

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
dataset="A0501_A0508"
dataset="A1020"

fName=paste("_",dataset,sep="")

datadir <- "data/"
phen <- read.table(file = paste(datadir,"phen",fName,".txt", sep=""), header = T, sep = "\t", quote="", comment.char="", as.is=T)
#phen=phen[1:4,]

fileList=paste(datadir,"cel/",dataset,"/",phen$celFileName,sep="")

if (F) {
datadir="tmp/"
fileList=dir(datadir,pattern="CEL")
phen=data.frame(id=sapply(fileList,function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F),celFileName=fileList,type=c(0,0,1,1),stringsAsFactors=F)
fileList=paste(datadir,fileList,sep="")
}

affyGeneFS <- read.celfiles(filenames=fileList)
phen=phen[match(colnames(affyGeneFS),phen$celFileName),]
sampleNames(affyGeneFS)=phen$id

if (F) {
    R> library(oligoData)
    R> data(affyGeneFS)
    R> affyGeneFS
}

## ---------------------------
fit <- fitProbeLevelModel(affyGeneFS)

if (F) {
png(paste("image",fName,"_%03d.png",sep=""),width = 3*240, height = 2*240)
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

genePS <- rma(affyGeneFS, target = "probeset")
geneCore <- rma(affyGeneFS, target = "core")

if (F) {
    R> featureData(exonPS) <- getNetAffx(exonPS, "probeset")
    R> featureData(exonCore) <- getNetAffx(exonCore, "transcript")
    R> featureData(exonFull) <- getNetAffx(exonFull, "transcript")
    R> featureData(exonExtd) <- getNetAffx(exonExtd, "transcript")
    R> featureData(geneCore) <- getNetAffx(geneCore, "transcript")

    R> exonCore
    R> featureData(exonCore)
    R> varLabels(featureData(exonCore))
    pData(featureData(exonCore))[1:2, "geneassignment"]
}
featureData(genePS) <- getNetAffx(genePS, "probeset")

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
png(paste("MAplot",fName,"_%03d.png",sep=""),width = 3*240, height = 2*240)
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

    anno <- read.table(file=paste(dirSrc,"Affymetrix/MoGene1_0st/MoGene-1_0-st-v1.na31.mm9.transcript.txt",sep=""), header=T, sep="\t", quote="", comment.char="", as.is=T)
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
    anno2 <- read.table(file=paste(dirSrc,"Affymetrix/MoGene1_0st/anno_MoGene-1_0-st-v1.na31.mm9.transcript.txt",sep=""), header=T, sep="\t", quote="", comment.char="", as.is=T)

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


datadir <- ""
tbl1 <- read.table(file = paste(datadir,"tmp1.txt", sep=""), header = F, sep = "\t", quote="", comment.char="", as.is=T)
tbl2 <- read.table(file = paste(datadir,"tmp2.txt", sep=""), header = F, sep = "\t", quote="", comment.char="", as.is=T)

x1=tbl1[,1]
x2=tbl2[,1]


write.table(paste("qRscript tmp_mmu.R ",x1[!x1%in%x2],sep=""), file="tmp3.txt",col.names=F,row.names=F, sep="\t",quote=F)

/home/royr/project/barcellosHoffMH/p53/htseqCountData/tmp,tophat2,gtf/Mus_musculus
/home/royr/project/barcellosHoffMH/p53/tophat2Data/tmp,tophat2,gtf/Mus_musculus

/cbc/data/fastqData/tmp/Mus_musculus


#######################################################
#######################################################

library(limma)
library(qvalue)
library(sva)
source(paste(dirSrc,"functions/TTest.9.1.6.R",sep=""))


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
contrasts(grp) <- contr.sum(sum(!duplicated(grp)))
#contrasts(grp) <- contr.treatment(sum(!duplicated(grp)))
design=model.matrix(~0+grp,ref="Sham")
colnames(design)=sub("(Intercept)","intercept",colnames(design),fixed=T)
fit <- lmFit(eset$expr[clId,],design)
contMat <- makeContrasts(gammaVsham=grp_-grpSham,hzeVsham=grpHZE-grpSham,hzeVgamma=grpHZE-grp_,levels=design)
#contMat <- makeContrasts(gammaVsham=grp1,hzeVsham=grp2,hzeVgamma=grp2-grp1,levels=design)
#contMat <- makeContrasts(gammaVsham=grp2,hzeVsham=grp3,hzeVgamma=grp3-grp2,levels=design)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)

grp=as.integer(eset$phen$group!="Sham")
grp[which(eset$phen$group=="HZE")]=2
compList=c("gammaVsham","hzeVsham","hzeVgamma")
tmp=matrix(nrow=nrow(fit2$coef),ncol=length(compList),dimnames=list(rownames(fit2$coef),compList))
fit3=list(coef=tmp,p.value=tmp)
fit4=list(coef=tmp,p.value=tmp)
for (k in 1:length(compList)) {
    compFlag=compList[k]
    switch(compFlag,
    "gammaVsham"={
        samId=which(grp%in%c(0,1))
    },
    "hzeVsham"={
        samId=which(grp%in%c(0,2))
    },
    "hzeVgamma"={
        samId=which(grp%in%c(1,2))
    }
    )
    phen=cbind(eset$phen[samId,],grp=grp[samId])
    design=model.matrix(~grp+experimentNo,data=phen)
    fit <- lmFit(eset$expr[clId,samId],design)
    fit <- eBayes(fit)
    fit3$coef[,k]=fit$coef[,2]
    fit3$p.value[,k]=fit$p.value[,2]
    
    mod = model.matrix(~as.factor(grp)+as.factor(experimentNo), data=phen)
    mod0 = model.matrix(~as.factor(experimentNo), data=phen)
    svObj=sva(eset$expr[clId,samId],mod,mod0)
    dat=svObj$sv; colnames(dat)=paste("sv",1:ncol(dat),sep="")
    design = cbind(mod,dat)
    fit <- lmFit(eset$expr[clId,samId],design)
    fit <- eBayes(fit)
    fit4$coef[,k]=fit$coef[,2]
    fit4$p.value[,k]=fit$p.value[,2]
}

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

## ----------------------------------------------
fName1="_moGene2.0"
fName1="_sva_moGene2.0"

colIdPV="pv"; colNamePV="PV"; pThres=10^-6
colIdPV="qv"; colNamePV="QV"; pThres=0.05

pThres2=0.05

onePlotFlag=F
onePlotFlag=T

compFlag3=""

compList1=c("gammaVsham","hzeVsham","hzeVgamma")
compList1=c("gammaVsham","hzeVsham")

switch(fName1,
    "_moGene2.0"={
        fitThis=fit3
    },
    "_sva_moGene2.0"={
        fitThis=fit4
    }
)

out=NULL
if (compFlag3=="") {
    if (onePlotFlag) {
        if (compFlag3=="") {
            #png(paste("dePlots",fName1,".png",sep=""),width=6*240,height=3*240)
            png(paste("dePlots",fName1,".png",sep=""),width=length(compList1)*2*240,height=3*240)
            par(mfcol=c(2,length(compList1)))
        } else {
            png(paste("dePlots",fName1,".png",sep=""),width=6*240,height=2*240)
            par(mfcol=c(2,6))
        }
    }
}
for (compFlag in compList1) {
    if (compFlag3=="") compList=compFlag
    if (compFlag3!="") {
        if (onePlotFlag) {
            if (compFlag3=="") {
                png(paste("dePlots",compFlag3,"_",compFlag,fName1,".png",sep=""),width=4*240,height=2*240)
                par(mfcol=c(2,4))
            } else {
                png(paste("dePlots",compFlag3,"_",compFlag,fName1,".png",sep=""),width=length(compList)*240,height=2*240)
                par(mfcol=c(2,length(compList)))
            }
        }
    }
    for (compThis in compList) {
        compFlag2=compFlag
        colGeneId="probesetid"
        stat_1=data.frame(probesetid=rownames(fitThis$coef),coef=fitThis$coef[,compFlag],pv=fitThis$p.value[,compFlag],qv=rep(NA,nrow(fitThis$coef)))
        i=which(!is.na(stat_1$pv))
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
        }
        )
        if (compFlag3!="") compName=paste(compThis,", ",compName,sep="")
        i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
        ann_1=ann_1[i1,]
        cat("\n\n",compName,"\n",sep="")
        x=table(stat_1$coef>0,stat_1[,colIdPV]<pThres)
        rownames(x)=c("Down","Up")[match(rownames(x),c("FALSE","TRUE"))]
        colnames(x)=paste(colNamePV,c(">=","<"),pThres,sep="")[match(colnames(x),c("FALSE","TRUE"))]
        print(x)
        
        cexMain=cexLab=1.5
        cexMain=cexLab=2.5
        cexMain=2.5; cexLab=1.5
        
        if (!onePlotFlag) png(paste("histogram_",compFlag,".png",sep=""))
        hist(stat_1$pv,xlab="P-value",main=compName,pch=19,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
        if (!onePlotFlag) dev.off()
        
        if (!onePlotFlag) png(paste("volcanoPlot_",compFlag,"_",tolower(colNamePV),pThres,".png",sep=""))
        i=which(stat_1[,colIdPV]<pThres)
        plot(stat_1$coef,-log10(stat_1$pv),xlab="Log2 fold change",ylab="-Log10(p-value)",main=paste(compName,"\nNo. with ",tolower(colNamePV)," < ",pThres,": ",length(i),sep=""),pch=19,cex=.8,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
        if (length(i)!=0) points(stat_1$coef[i],-log10(stat_1$pv[i]),col="red",pch=19,cex=.8)
        if (!onePlotFlag) dev.off()

        out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("coef","pv","qv")]
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
        names(out2)=paste(names(out2),"_",compFlag2,sep="")
        names(out2)=sub("coef_","log2FC",names(out2))
        out2=cbind(ann_1[i,which(!names(ann_1)%in%colGeneId)],out2)
        write.table(out2, file=paste("stat_",compFlag2,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)

        if (F) {
            out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("coef","pv","qv")]
            i=order(out2$pv); i=i[which(out2$pv[i]<pThres2)]
            out2=out2[i,]
            names(out2)=paste(names(out2),"_",compFlag2,sep="")
            out2=cbind(ann_1[i,which(!names(ann_1)%in%colGeneId)],out2)
            write.table(out2, file=paste("stat_",compFlag2,"_pv",pThres2,fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
        }

        if (F) {
            write.table(unique(ann_1$geneId), file=paste("ensGeneId_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
            i=which(stat_1$pv<pThres)
            if (length(i)!=0) write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_pv",pThres,"_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
            i=order(stat_1$pv)[1:200]
            write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_top200_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
        }
    }
    if (onePlotFlag & compFlag3!="") dev.off()
}
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

## ----------------------------------------------

##############################################
library(sva)

load("results/eset.RData")

modCom = model.matrix(~as.factor(group), data=eset$phen)

exprCom=ComBat(dat=eset$expr, batch=eset$phen$experimentNo, mod=modCom, par.prior = TRUE,prior.plots = FALSE)
#exprComNP=ComBat(dat=eset$expr, batch=eset$phen$experimentNo, mod=modCom, par.prior = FALSE,prior.plots = FALSE)

mod = model.matrix(~as.factor(group)+as.factor(experimentNo), data=eset$phen)
mod0 = model.matrix(~as.factor(experimentNo), data=eset$phen)


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
