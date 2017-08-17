## Run deGenes_rnaseq.R section 1 first
## Use Voom+sva with minCnt=1

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

##############################################
## Section 1

datadir="/Users/royr/Downloads/tmp/"
datadir=""
datadir="results/funcAnno/"
#tbl=read.table(paste(datadir,"ann_TransgenicMouseVsWildtypeMouse_fromDavid.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
datadir="results/rnaSeq/"
datadir="results/rnaSeq/misc/"

#phen=read.table(paste("docs/sampleInfo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
cnt10=read.table(paste(datadir,"count_raw_Mus_musculus.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,row.names="geneId")
phen10=read.table(paste(datadir,"sample_Mus_musculus.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

cnt_1=t(t(cnt10[,match(phen10$id,colnames(cnt10))])*phen10$norm.factors)
write.table(cbind(geneId=rownames(cnt_1),cnt_1), file="count_norm_Homo_sapiens.txt",col.names=T,row.names=F, sep="\t",quote=F)

if (F) {
    source("biomartify.R")
    source("funcs1.R")
    top=stat1_36[1:10,]
    rownames(top)=top$geneId
    top1=biomartify(top, organism = "HomoSapiens")
}

library(biomaRt)

if (F) {
    marts=listMarts()
    index<-grep("ensembl",marts)
    mart=useMart(marts[index])
    mart=useMart("ENSEMBL_MART_MOUSE")
    mart="hsapiens_gene_ensembl"
    mart="mmusculus_gene_ensembl"
    listDatasets(mart = mart)
    #martDisconnect(mart = mart)
}

tmpC=rep("",nrow(cnt10))
ann10=data.frame(geneId=rownames(cnt10),geneSym=tmpC,stringsAsFactors=F)
dataset="mmusculus_gene_ensembl"
ensembl=useMart("ensembl", dataset = dataset)
attrs=c("ensembl_gene_id", "mgi_symbol", "mgi_id", "entrezgene", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description")
res=getBM(attributes=attrs, filters = "ensembl_gene_id", values = rownames(cnt10), mart=ensembl, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE)
i=match(rownames(cnt10),res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
if (length(i12)!=0) {
    tmp=res[1:length(i12),]
    for (k in 1:ncol(tmp)) {
        tmp[,k]=NA
    }
    tmp$ensembl_gene_id=rownames(cnt10)[i12]
    res=rbind(res,tmp)
}
id=unique(res[,1][duplicated(res[,1])])
for (k in 1:length(id)) {
    i=which(res[,1]==id[k])
    if (sum(!duplicated(res$chromosome_name[i]))!=1 | sum(!duplicated(res$start_position[i]))!=1 | sum(!duplicated(res$end_position[i]))!=1) cat(k,"\n")
}
## No replicate gene IDs with different position
i=match(rownames(cnt10),res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]
res=res[i2,]
ann10=res
names(ann10)[match(c("ensembl_gene_id"),names(ann10))]="geneId"
write.table(ann10, file=paste("ann_Mus_musculus.txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)

#alcpmT=aveLogCPM(dgeT, normalized.lib.sizes=TRUE)


##############################################
## Section 2
## DE genes
## Set params datadir

#library(clusterProfiler)
library(qvalue)
source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))

dir("results/rnaSeq/deGene",pattern="stat_")

plotFlag="plotAveLogCPMvsVoomWilcoxPV"
plotFlag="scatterPlot_topGenes"
plotFlag="dePlots"

adjFlag="treatExptNoAdj"
adjFlag="treatAndExptNoAdj"

fName1="_rnaSeq_p53"
datadir="results/rnaSeq/deGene/"
datadir="./"
datadir="results/rnaSeq/deGene/slope/"
datadir="results/rnaSeq/slopeCat/"
datadir="results/rnaSeq/slope/"
datadir="results/rnaSeq/deGene/treatment/"
datadir="results/rnaSeq/treatAndExptNoAdj/slope/"
datadir="results/rnaSeq/deGene/slopeCat/"

colIdPV="pv"; colNamePV="PV"; pThres=10^-6
colIdPV="qv"; colNamePV="QV"; pThres=0.05

numPlotFlag="single"
numPlotFlag="subset"
numPlotFlag="all"

minCntFlag=c(10)
compFlag3=""

minCntFlag=c(0,2,5,10,20)
minCntFlag=c(0,5,10,20,50)
minCntFlag=c(0,1,5,10,20)
minCntFlag=c(1)
compFlag3="_edgerVsVoom"
compFlag3="_edgeR"
compFlag3="_wilcox_voom"
compFlag3="_voom"

## --------------
minCntFlag=c(0,1,2,5,10,20,50)
minCntFlag=c(10)
minCntFlag=c(1)
compFlag3="_voom_sva"
compFlag3="_voom"
compFlag3="_sva"

#subsetList=""; compFlag1=c("_A1026VsA1014"); compThis="wilcox"
#subsetList="_A1014A1402"; compFlag1=c("_gammaVsSham"); compThis="wilcox"

compThis0="wilcox"

varType="categorical"
compList0=paste("_",c("fast"),"VsSlow",sep=""); subsetList=c("",paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe","sham"),sep=""))

## --------------
minCntFlag=c(0,1,2,5,10,20,50)
minCntFlag=c(10)
minCntFlag=c(1)
compFlag3="_voom_sva"
compFlag3="_voom"
compFlag3="_sva"

#subsetList=""; compFlag1=c("_A1026VsA1014"); compThis="wilcox"
#subsetList="_A1014A1402"; compFlag1=c("_gammaVsSham"); compThis="wilcox"

compThis0="wilcox"

varType="categorical"
compList0="_A1026VsA1014"; subsetList=""
compList0=paste("_",c("capeGamma"),"VsSham",sep=""); subsetList=""
compList0="_gammaVsSham"; subsetList="_A1014A1402"
compList0=paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe"),"VsSham",sep=""); subsetList=""
compList0=paste("_",c("capeGamma"),"VsGamma",sep=""); subsetList=""
compList0=paste("_",c("capeGamma"),"VsGamma",sep=""); subsetList="_A1402"

## --------------
minCntFlag=c(0,1,2,5,10,20,50)
minCntFlag=c(10)
minCntFlag=c(1)
compFlag3="_voom_sva"
compFlag3="_voom"
compFlag3="_sva"
compThis0=NULL
varType="continuous"
compList0="_slope"
subsetList=c("",paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe","sham"),sep=""))
#subsetList=""

## --------------
minCntFlag=c(0,1,2,5,10,20,50)
minCntFlag=c(10)
minCntFlag=c(1)
compFlag3=""
compFlag3="_voom_sva"
compFlag3="_sva"
compFlag3="_voom"

compThis0="wilcox"

varType="categorical"
compList0="_fastVsSlow"
subsetList=c("",paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe","sham"),sep=""))

## --------------

## --------------
#for (minCntFlag in c(0,1,2,5,10,20,50)) {
#for (minCntFlag in c(1)) {
    
    fName2=ifelse(length(compList0)==1,compList0,subsetList)
    if (fName2=="_") fName2=""
    fName2=paste(fName2,compFlag3,"_minCnt",minCntFlag,fName1,sep="")
    #fName2=fName1
    
    if (F) {
        if (numPlotFlag=="all" & plotFlag=="dePlots") {
            png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,".png",sep=""),width=length(compList)*240,height=2*240)
            par(mfcol=c(2,length(compList)))
        }
    }
    if (numPlotFlag=="subset" & plotFlag=="dePlots") {
        png(paste(plotFlag,compFlag3,"_minCnt",minCntFlag,fName1,".png",sep=""),width=length(subsetList)*240,height=2*240)
        par(mfcol=c(2,length(subsetList)))
    }
    for (subsetFlag in subsetList) {
        for (compFlag0 in compList0) {
            if (compFlag3!="") {
                x=dir(datadir,pattern=paste("stat",compFlag0,sep=""))
                x=x[grep(paste(subsetFlag,fName1,sep=""),x)]
                x=x[grep(compFlag3,x)]
                if (length(x)==0) {
                    next
                }
            }
            if (!is.null(compThis0)) {
                x=dir(datadir,pattern=paste("stat",compFlag0,sep=""))
                x=x[grep(paste(subsetFlag,fName1,sep=""),x)]
                x=x[grep("wilcox",x)]
                if (length(x)==0) {
                    next
                }
                statW=read.table(paste(datadir,"stat",compFlag0,subsetFlag,fName1,"_",compThis0,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                names(statW)=c("geneId","log2FC","pv")
                statW$qv=NA
                i=which(!is.na(statW$pv))
                statW$qv[i]=qvalue(statW$pv[i])$qvalues
            }
            
            compList1=compFlag0
            
            out=NULL
            if (compFlag3=="") {
                if (numPlotFlag=="all") {
                    if (compFlag3=="") {
                        if (plotFlag=="plotAveLogCPMvsVoomWilcoxPV") {
                            png(paste(plotFlag,subsetFlag,fName2,".png",sep=""),width=3*240,height=2*240)
                            par(mfcol=c(2,3))
                        } else {
                            png(paste(plotFlag,subsetFlag,fName2,".png",sep=""),width=4*240,height=2*240)
                            par(mfcol=c(2,4))
                        }
                    } else {
                        if (plotFlag=="plotAveLogCPMvsVoomWilcoxPV") {
                            png(paste(plotFlag,subsetFlag,fName2,".png",sep=""),width=3*240,height=2*240)
                            par(mfcol=c(2,3))
                        } else {
                            png(paste(plotFlag,subsetFlag,fName2,".png",sep=""),width=6*240,height=2*240)
                            par(mfcol=c(2,6))
                        }
                    }
                }
            }
            for (compFlag1 in compList1) {
                if (compFlag3=="") compList=compFlag1
                if (compFlag3=="_edgerVsVoom") compList=c("edgeR",paste("voom_minCnt",minCntFlag,sep=""))
                if (compFlag3=="_voom") compList=paste("voom_minCnt",minCntFlag,sep="")
                if (compFlag3=="_edgeR") compList=paste("edgeR_minCnt",minCntFlag,sep="")
                if (compFlag3=="_wilcox_voom") compList=c("wilcox",paste("voom_minCnt",minCntFlag,sep=""))
                if (compFlag3=="_voom_sva") compList=paste("voom",c("","_sva"),"_minCnt",minCntFlag,sep="")
                if (compFlag3=="_sva") compList=paste("voom_sva_minCnt",minCntFlag,sep="")
                if (compFlag3!="") {
                    if (numPlotFlag=="all") {
                        if (compFlag3=="") {
                            if (plotFlag=="plotAveLogCPMvsVoomWilcoxPV") {
                                png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,".png",sep=""),width=3*240,height=2*240)
                                par(mfcol=c(2,3))
                            } else if (plotFlag=="dePlots") {
                                png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,fName2,".png",sep=""),width=4*240,height=2*240)
                                par(mfcol=c(2,4))
                            }
                        } else {
                            if (plotFlag=="plotAveLogCPMvsVoomWilcoxPV") {
                                ##png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,".png",sep=""),width=length(compList)*240,height=2*240)
                                ##par(mfcol=c(2,ceiling(length(compList)/2)))
                                #png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,".png",sep=""),width=(length(compList))*2*240,height=2*240)
                                #par(mfcol=c(1,length(compList)))
                                png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,fName2,".png",sep=""),width=5*240,height=1*240)
                                par(mfrow=c(1,5))
                            } else if (plotFlag=="dePlots") {
                                png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,fName2,".png",sep=""),width=length(compList)*240,height=2*240)
                                par(mfcol=c(2,length(compList)))
                                #png(paste(plotFlag,compFlag3,compFlag1,subsetFlag,".png",sep=""),width=5*240,height=2*240)
                                #par(mfcol=c(2,5))
                            }
                        }
                    }
                }
                for (compThis in compList) {
                    compFlag2=compFlag1
                    colGeneId="geneId"
                    #stat_1=read.table(paste(datadir,"stat",compFlag1,subsetFlag,fName1,"_",compThis,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    x=dir(datadir,pattern=paste("stat",compFlag1,sep=""))
                    x=x[grep(paste(subsetFlag,fName1,"_",compThis,".txt",sep=""),x)]
                    #x=x[grep(paste(subsetFlag,fName1,compThis,".txt",sep=""),x)]
                    if (length(x)>1) {
                        x=x[grep(paste("_",c(sub("_","",adjFlag),"exptNoAdj","investAdj","treatAdj"),subsetFlag,fName1,sep="",collapse="|"),x)]
                    }
                    if (length(x)==0) {
                        next
                    }
                    cat("\n\n===========",x,"==============\n")
                    #stat_1=read.table(paste(datadir,"stat",compFlag1,subsetFlag,fName1,"_",compThis,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    stat_1=read.table(paste(datadir,x,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    names(stat_1)=c("geneId","log2FC","pv")
                    ann_1=ann10
                    i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
                    samId=1:nrow(phen10)
                    subsetName=ifelse(compFlag1=="_slope","All samples: ","")
                    switch(subsetFlag,
                    "_A1402"={
                        samId=which(phen10$experimentNo%in%c("A1402"))
                        subsetName="A1402: "
                    },
                    "_A1014A1402"={
                        samId=which(phen10$experimentNo%in%c("A1014","A1402"))
                        subsetName="A1014, A1402: "
                    },
                    "_cape"={
                        samId=which(phen10$treatment2%in%c("CAPE"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_capeGamma"={
                        samId=which(phen10$treatment2%in%c("CAPE+ gamma-rays"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_gamma"={
                        samId=which(phen10$treatment2%in%c("gamma-rays"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_gammaTe"={
                        samId=which(phen10$treatment2%in%c("gamma-rays-TE"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_hze"={
                        samId=which(phen10$treatment2%in%c("HZE"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_hzeTe"={
                        samId=which(phen10$treatment2%in%c("HZE-TE"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    },
                    "_sham"={
                        samId=which(phen10$treatment2%in%c("Sham"))
                        subsetName=paste(unique(phen10$treatment2[samId]),": ",sep="")
                    }
                    )
                    grpName=NULL
                    switch(compFlag1,
                    "_A1026VsA1014"={
                        compName="A1026 vs. A1014"
                        grpUniq=c("A1014","A1026")
                        varList="experimentNo"
                    },
                    "_capeVsSham"={
                        compName="CAPE vs. sham"
                        grpUniq=c("Sham","gamma-rays")
                        varList="treatment2"
                    },
                    "_capeGammaVsSham"={
                        compName="CAPE+gamma vs. sham"
                        grpUniq=c("Sham","CAPE+ gamma-rays")
                        varList="treatment2"
                    },
                    "_gammaVsSham"={
                        compName="Gamma vs. sham"
                        grpUniq=c("Sham","gamma-rays")
                        varList="treatment2"
                    },
                    "_gammaTeVsSham"={
                        compName="Gamma-TE vs. sham"
                        grpUniq=c("Sham","gamma-rays")
                        varList="treatment2"
                    },
                    "_hzeVsSham"={
                        compName="HZE vs. sham"
                        grpUniq=c("Sham","gamma-rays")
                        varList="treatment2"
                    },
                    "_hzeTeVsSham"={
                        compName="HZE-TE vs. sham"
                        grpUniq=c("Sham","gamma-rays")
                        varList="treatment2"
                    },
                    "_capeGammaVsGamma"={
                        compName="CAPE+gamma vs. gamma"
                        grpUniq=c("gamma-rays","CAPE+ gamma-rays")
                        varList="treatment2"
                    },
                    "_slope"={
                        compName="Slope"
                        grpUniq=NULL
                        varList="slope"
                    },
                    "_fastVsSlow"={
                        compName="Fast vs. slow"
                        grpName=c("slow","fast")
                        grpUniq=c(0,1)
                        varList="slopeCat2"
                        varList="slopeCat3"
                    }
                    )
                    if (is.null(grpUniq)) {
                        samId=samId[!is.na(phen10[samId,varList])]
                    } else {
                        if (is.null(grpName)) grpName=grpUniq
                        samId=samId[which(phen10[samId,varList]%in%grpUniq)]
                    }
                    stat_1$qv=NA
                    i=which(!is.na(stat_1$pv))
                    stat_1$qv[i]=qvalue(stat_1$pv[i])$qvalues
                    #if (compFlag3!="") compName=paste(compThis,", ",compName,sep="")
                    header=paste(subsetName,compName,sep="")
                    header=paste(subsetName,ifelse(subsetFlag=="","","\n"),compName,sep="")
                    fName3=paste(compFlag2,subsetFlag,fName2,"_",compThis,sep="")
                    fName3=paste(compFlag2,subsetFlag,fName2,sep="")
                    if (!is.null(grpUniq)) {
                        if (all(is.na(stat_1$log2FC))) {
                            stat_1$log2FC=apply(lcpmT[,which(phen10[,varList]==grpUniq[2])],1,mean,na.rm=T)-apply(lcpmT[,which(phen10[,varList]==grpUniq[1])],1,mean,na.rm=T)
                        }
                    }
                    if (plotFlag=="plotAveLogCPMvsVoomWilcoxPV") {
                        lim=c(-3,3)
                        #plot(alcpmT[i1],log10(stat_1$pv)-log10(statW$pv[i1]),ylim=lim,main=paste(header,"\nN = ",nrow(stat_1),sep=""),xlab="Avg logCPM",ylab="log10(PV-limma)-log10(PV-wilcox)",cex=.2,pch=20,cex.main=2.5,cex.axis=2,cex.lab=1.5)
                        plot(alcpmT[i1],log10(stat_1$pv)-log10(statW$pv[i1]),ylim=lim,main=paste(header,"\nN = ",nrow(stat_1),sep=""),xlab="Avg logCPM",ylab="log10(PV-limma)-log10(PV-wilcox)",cex=.2,pch=20,cex.main=1,cex.axis=1.5,cex.lab=1.5)
                        abline(h=0,col="red")
                    } else if (plotFlag=="scatterPlot_topGenes") {
                        for (signifFlag in c("_topGene","_signifWilcoxNotVoom","_signifVoomNotWilcox")) {
                            if ((is.null(compThis0) | compThis=="wilcox") & signifFlag%in%c("_signifWilcoxNotVoom","_signifVoomNotWilcox")) next
                            switch(signifFlag,
                            "_topGene"={
                                #ii=order(stat_1$pv)[1:64]
                                ii=order(stat_1$pv,decreasing=F)
                                ii=ii[which(stat_1$qv[ii]<pThres)]
                                if (length(ii)!=0) ii=ii[1:min(64,length(ii))]
                            },
                            "_signifWilcoxNotVoom"={
                                ii=order(stat_1$pv,decreasing=T)
                                #ii=ii[which(statW[i1[ii],colNamePV]<pThres)]
                                ii=ii[which(statW$pv[i1[ii]]<pThres)]
                                ii=ii[1:min(64,length(ii))]
                            },
                            "_signifVoomNotWilcox"={
                                ii=order(statW$pv[i1],decreasing=T)
                                #ii=ii[which(stat_1[ii,colNamePV]<pThres)]
                                ii=ii[which(stat_1$pv[ii]<pThres)]
                                ii=ii[1:min(64,length(ii))]
                            }
                            )
                            if (length(ii)==0) next
                            subDir=""
                            if (T) {
                                subDir=paste(compThis,signifFlag,sep="")
                                if (!file.exists(subDir)){
                                    dir.create(file.path(subDir))
                                }
                                subDir=paste(subDir,"/",sep="")
                            }
                            #png(paste(subDir,plotFlag,"_",compThis,compFlag1,subsetFlag,"_%1d.png",sep=""),width=8*2*240,height=4*2*240)
                            png(paste(subDir,plotFlag,"_",compThis,compFlag1,subsetFlag,"_%1d.png",sep=""),width=4*240,height=2*240)
                            par(mfcol=c(2,4))
                            for (i in ii) {
                                lim=NULL
                                if (varType=="categorical") {
                                    ##boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
                                    #boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=1.5,cex.lab=1.5,cex.axis=1)
                                    boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],names=grpName,ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]]," ",ann10$mgi_symbol[i1[i]],"\nlog2FC ",round(stat_1$log2FC[i],2),", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=1.5,cex.lab=1.5,cex.axis=1)
                                    for (grpId in 1:length(grpUniq)) {
                                        j=which(phen10[samId,varList]==grpUniq[grpId])
                                        #points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=5)
                                        points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=1)
                                        lines(grpId+c(-.4,.4),rep(mean(lcpmT[i1[i],samId[j]]),2),col="green")
                                    }
                                } else {
                                    ##plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
                                    #plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",pch=19,cex=.8,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
                                    plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]]," ",ann10$mgi_symbol[i1[i]],"\nlog2FC ",round(stat_1$log2FC[i],2),", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",pch=19,cex=.8,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
                                }
                            }
                            dev.off()
                        }
                    } else {
                        #i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
                        ann_1=ann_1[i1,]
                        cat("\n\n",header,"\n",sep="")
                        x=table(stat_1$log2FC>0,stat_1[,colIdPV]<pThres)
                        rownames(x)=c("Down","Up")[match(rownames(x),c("FALSE","TRUE"))]
                        colnames(x)=paste(colNamePV,c(">=","<"),pThres,sep="")[match(colnames(x),c("FALSE","TRUE"))]
                        print(x)
                        
                        if (F) {
                            lim=c(0,1)
                            png(paste("tmp_",fName,".png",sep=""))
                            plot(stat$pvBeta[i],stat$qvBeta[i],xlim=lim,ylim=lim,xlab="P-value",ylab="Q-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
                            abline(0,1)
                            dev.off()
                            
                            png(paste("qqplot_",fName,".png",sep=""))
                            pvs=sort(na.exclude(stat$pvBeta[i]))
                            qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
                            abline(0,1)
                            dev.off()
                        }
                        
                        if (numPlotFlag=="single") png(paste("histogram_",compFlag1,subsetFlag,".png",sep=""))
                        hist(stat_1$pv,xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
                        if (numPlotFlag=="single") dev.off()
                        
                        if (numPlotFlag=="single") png(paste("volcanoPlot_",compFlag1,subsetFlag,"_",tolower(colNamePV),pThres,".png",sep=""))
                        i=which(stat_1[,colIdPV]<pThres)
                        plot(stat_1$log2FC,-log10(stat_1$pv),xlab="Log2 fold change",ylab="-Log10(p-value)",main=paste(header,"\nN=",nrow(stat_1),", no. with ",tolower(colNamePV),"<",pThres,": ",length(i),sep=""),pch=19,cex=.1,cex.axis=1.5,cex.lab=1.5)
                        if (length(i)!=0) points(stat_1$log2FC[i],-log10(stat_1$pv[i]),col="red",pch=19,cex=.2)
                        if (numPlotFlag=="single") dev.off()
                        
                        #out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("foldChange","log2FC","logCPM","pv","FDR")]
                        out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("log2FC","pv","qv")]
                        i=order(out2$pv)
                        #i=order(out2$pv); i=i[out2[i,colIdPV]<pThres]
                        out2=out2[i,]
                        names(out2)=paste(names(out2),compFlag2,sep="")
                        out2=cbind(ann_1[i,which(!names(ann_1)%in%c("gene_id"))],out2)
                        for (k in 1:ncol(out2)) {
                            if (is.character(out2[,k])) {
                                out2[,k]=gsub(", ...","",out2[,k])
                            }
                        }
                        write.table(out2, file=paste("stat",fName3,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
                        
                        if (F) {
                            write.table(unique(ann_1$geneId), file=paste("ensGeneId_",compFlag1,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
                            i=which(stat_1[,colIdPV]<pThres)
                            if (length(i)!=0) write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_pv",pThres,"_",compFlag1,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
                            i=order(stat_1[,colIdPV])[1:200]
                            write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_top200_",compFlag1,subsetFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
                        }
                    }
                }
                if (numPlotFlag=="all" & compFlag3!="" & plotFlag!="scatterPlot_topGenes") dev.off()
            }
            if (numPlotFlag=="all" & compFlag3=="" & plotFlag!="scatterPlot_topGenes") dev.off()
        }
    }
    if (numPlotFlag=="subset" & plotFlag=="dePlots")  dev.off()
        
#}

## Check log2 fold change direction
j1=which(phen12$type=="untreated")
j2=which(phen12$type=="TGFbeta")
log2Expr=cbind(wt=apply(log2(cnt12[i1,j1]),1,mean,na.rm=T),tr=apply(log2(cnt12[i1,j2]),1,mean,na.rm=T))
i=which(stat_1$pv<pThres)
cbind(stat_1$log2FC[i],log2Expr[i,])[which(stat_1$log2FC[i]>0),][1:10,]
cbind(stat_1$log2FC[i],log2Expr[i,])[which(stat_1$log2FC[i]<0),][1:10,]

## Check distribution of top genes
## Top 2 genes from cape+gamma vs. gamma for all samples, adjusted for experiment no.
## Not used because cape+gamma was done in A1402 only while gamma was done in A1014 &
## A1402, leading to +ve association with cape+gamma for ENSMUSG00000001305 when the
## association is -ve when considering A1402 samples only
png(paste("scatterPlot_2TopGenes_capeGammaVsGamma_exprNoAdj.png"))
ii=ii[1:2]
x1=sub("+","\n",paste(phen10$treatment2[samId],phen10$experimentNo[samId],sep="\n"),fixed=T)
grp2Uniq=sort(unique(x1))
par(mfcol=c(1,2))
for (i in ii) {
    lim=NULL
    if (varType=="categorical") {
        boxplot(lcpmT[i1[i],samId]~x1,names=grpName,ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],"\nlog2FC ",round(stat_1$log2FC[i],2),", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=1,cex.lab=1.5,cex.axis=.7,las=3)
        for (grpId in 1:length(grp2Uniq)) {
            j=which(x1==grp2Uniq[grpId])
            points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=1)
            lines(grpId+c(-.4,.4),rep(mean(lcpmT[i1[i],samId[j]]),2),col="green")
        }
    }
}
dev.off()
png(paste("scatterPlotLegend.png"))
sampleColorLegend(tls="mean",col="green",lty="solid")
dev.off()

png(paste("densityPlot_slope.png"))
plot(density(dgeT$samples$slope,na.rm=T),main="slope")
dev.off()

## NOT USED
par(mfcol=c(1,2))
for (i in ii) {
    lim=NULL
    if (varType=="categorical") {
        ##boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
        #boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=1.5,cex.lab=1.5,cex.axis=1)
        boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]]," ",ann10$mgi_symbol[i1[i]],"\nlog2FC ",round(stat_1$log2FC[i],2),", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=1.5,cex.lab=1.5,cex.axis=1)
        for (grpId in 1:length(grpUniq)) {
            j=which(phen10[samId,varList]==grpUniq[grpId])
            #points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=5)
            points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=1)
            lines(grpId+c(-.4,.4),rep(mean(lcpmT[i1[i],samId[j]]),2),col="green")
        }
    } else {
        ##plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
        #plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",pch=19,cex=.8,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
        plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]]," ",ann10$mgi_symbol[i1[i]],"\nlog2FC ",round(stat_1$log2FC[i],2),", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",pch=19,cex=.8,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
    }
}

## Scatter plot of log2FC for cape+gamma vs. gamma for all samples vs. A1402 samples
stat_12=stat_1
table(stat_11[,1]==stat_12[,1])
plot(stat_11$log2FC,stat_12$log2FC)
i=which(stat_11$qv<pThres)
points(stat_11$log2FC[i],stat_12$log2FC[i],col="red")


###########################################################
###########################################################


## ----------------------------------------------
library(clusterProfiler)
library(qvalue)

outFormat="png"
outFormat="pdf"

outFormat2="pdf"

fName1="_rnaSeq_p53"
minCntFlag=c(10)
minCntFlag=c(1)


compFlag3="_voom"
compFlag3="_sva"

organismThis=c("human","hsa")
orgDB="org.Hs.eg.db"
organismThis=c("mouse","mm")
orgDB="org.Mm.eg.db"

colIdEst="logFC"; colIdPV=c("PValue","FDR"); pThres=0.05
colIdEst="log2FC"; colIdPV=c("pv","qv"); pThres=0.05

compList0=paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe"),"VsSham",sep=""); subsetFlag=""
compList0="_A1026VsA1014"; subsetFlag=""
compList0="_gammaVsSham"; subsetFlag="_A1014A1402"
compList0=paste("_",c("capeGamma"),"VsSham",sep=""); subsetFlag=""
compList0=paste("_",c("capeGamma"),"VsSham",sep=""); subsetFlag="_A1402"

for (compFlag0 in compList0) {
    colInfo=data.frame(key=c("PValue","FDR"),value=c("pv","qv"),stringsAsFactors=F)
    compList1=compFlag0
    for (pThres in c(0.05)) {
        #for (pThres in c(NA)) {
        for (compFlag1 in compList1) {
            compFlag2=compFlag1
            if (compFlag3=="_voom") compList=paste("voom_minCnt",minCntFlag,sep="")
            if (compFlag3=="_sva") compList=paste("voom_sva_minCnt",minCntFlag,sep="")
            compThis=compList
            colGeneId="geneId"
            x=dir(datadir,pattern=paste("stat",compFlag1,sep=""))
            x=x[grep(paste(subsetFlag,fName1,"_",compThis,".txt",sep=""),x)]
            if (length(x)>1) {
                if (subsetFlag=="") {
                    x=x[grep("_exptNoAdj|investAdj",x)]
                } else {
                    x=x[grep(subsetFlag,x)]
                }
            }
            cat("===========",x,"==============\n")
            #stat_1=read.table(paste(datadir,"stat",compFlag1,subsetFlag,fName1,"_",compThis,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            stat_1=read.table(paste(datadir,x,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            names(stat_1)=c("geneId","log2FC","pv")
            ann_1=ann10
            i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
            samId=1:nrow(phen10)
            subsetName=""
            switch(subsetFlag,
                "_A1402"={
                    samId=which(phen10$experimentNo%in%c("A1402"))
                    subsetName="A1402: "
                },
                "_A1014A1402"={
                    samId=which(phen10$experimentNo%in%c("A1014","A1402"))
                    subsetName="A1014, A1402: "
                }
            )
            switch(compFlag1,
                "_A1026VsA1014"={
                    compName="A1026 vs. A1014"
                    grpUniq=c("A1014","A1026")
                    varList="experimentNo"
                },
                "_capeVsSham"={
                    compName="CAPE vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_capeGammaVsSham"={
                    compName="CAPE+gamma vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_gammaVsSham"={
                    compName="Gamma vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_gammaTeVsSham"={
                    compName="Gamma-TE vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_hzeVsSham"={
                    compName="HZE vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_hzeTeVsSham"={
                    compName="HZE-TE vs. sham"
                    grpUniq=c("Sham","gamma-rays")
                    varList="treatment2"
                },
                "_slope"={
                    compName="Slope"
                    grpUniq=NULL
                    varList="slope"
                }
            )
            if (!is.null(grpUniq)) {
                samId=samId[which(phen10[samId,varList]%in%grpUniq)]
            }
            stat_1$qv=NA
            i=which(!is.na(stat_1$pv))
            stat_1$qv[i]=qvalue(stat_1$pv[i])$qvalues
            if (compFlag3!="") compName=paste(compThis,", ",compName,sep="")
            cat("\n\n================== ",compFlag1,", ",pThres," ==================\n",sep="")
            for (pThres1 in c(0.05,0.01,0.001,0.0005)) {
                x=table(stat_1[,colIdPV[2]]<pThres1)
                names(x)=paste("QV",c(">=","<"),pThres1,sep="")[match(names(x),c("FALSE","TRUE"))]
                print(x)
            }
            
            
            iA2=match(stat_1$geneId,ann_1$geneId)
            
            if (F) {
                idType(OrgDb = orgDB)
                out=bitr(unique(ann_1$mgi_symbol[!is.na(ann_1$mgi_symbol)]), fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb=orgDB)
                table(is.na(as.integer(out$ENTREZID))) ## ALL FALSE
                ann_1$entrezId=""
                i=match(ann_1$mgi_symbol,out$SYMBOL); i1=which(!is.na(i)); i2=i[i1]
                ann_1$entrezId[i1]=out$ENTREZID[i2]
            }
            ann_1$entrezId=as.character(ann_1$entrezgene)
            ann_1$entrezId[is.na(ann_1$entrezId)]=""
            table(ann_1$entrezId=="")
            
            ann2=ann_1[iA2,]
            if (is.na(pThres)) {
                fName2=paste("_",compFlag1,subsetFlag,"_top200genes",sep="")
                prId=order(stat_1[,colIdPV[2]])[1:200]
            } else {
                fName2=paste("_",compFlag1,subsetFlag,"_",colInfo$value[which(colInfo$key==colIdPV[2])],pThres,sep="")
                prId=which(stat_1[,colIdPV[2]]<pThres)
            }
            if (length(prId)==0) next
            
            #for (enrichFlag in c("go","kegg","do","david")) {
            #for (enrichFlag in c("go","kegg")) {
            for (enrichFlag in c("go","kegg")) {
                ontList=""
                switch(enrichFlag,
                "go"={
                    ontList=c("BP","CC","MF")
                    pThres2=0.05; qThres2=0.05
                },
                "kegg"={
                    pThres2=0.4; qThres2=0.4
                    pThres2=0.05; qThres2=0.05
                },
                "do"={
                    pThres2=0.05; qThres2=0.05
                },
                "david"={
                }
                )
                pThres2=0.005; qThres2=0.05
                pThres2=0.05; qThres2=0.05
                
                for (ontThis in ontList) {
                    cat("\n\n================== ",compFlag1,", ",pThres,", ",enrichFlag,", ",ontThis," ==================\n",sep="")
                    i0=which(!duplicated(ann2$entrezId) & ann2$entrezId!="")
                    i=prId[which(!duplicated(ann2$entrezId[prId]) & ann2$entrezId[prId]!="")]
                    res=NULL
                    switch(enrichFlag,
                    "go"={
                        cat("------------- groupGO\n\n",sep="")
                        fName3=paste("groupGO",fName2,ifelse(ontThis=="","","_"),ontThis,sep="")
                        res=NULL
                        #res=groupGO(gene=ann2$entrezId[prId],organism=organismThis[1],ont=ontThis,level=3,readable=TRUE)
                        res=groupGO(gene=ann2$entrezId[i],OrgDb=orgDB,ont=ontThis,level=3,readable=TRUE)
                        #print(head(summary(res)))
                        if (outFormat=="png") {
                            png(paste(fName3,".png",sep=""))
                        } else if (outFormat=="pdf") {
                            pdf(paste(fName3,".pdf",sep=""))
                        }
                        #plot(1:5,1:5)
                        barplot(res, drop=TRUE, showCategory=12)
                        dev.off()
                        if (F) {
                            png(paste("tmp",fName2,ifelse(ontThis=="","","_"),ontThis,".png",sep=""))
                            plot(1:5,1:5)
                            dev.off()
                        }
                        res=NULL
                        res=try(enrichGO(gene=ann2$entrezId[i],keytype="ENTREZID",universe=ann2$entrezId[i0],OrgDb=orgDB,ont=ontThis,pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,readable=TRUE))
                        if (class(res) %in% c("try-error","NULL")) {
                            print(paste(ontThis,": Error",sep=""))
                        } else {
                            if (nrow(res@result)==0) {
                                print(paste(ontThis,": Nothing significant",sep=""))
                            } else {
                                res=simplify(res)
                            }
                        }
                    },
                    "kegg"={
                        res=enrichKEGG(gene=ann2$entrezId[i],organism=organismThis[2],universe=ann2$entrezId[i0],pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,use_internal_data=TRUE)
                    },
                    "do"={
                        res=enrichDO(gene=ann2$entrezId[i],ont="DO",universe=ann2$entrezId[i0],pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,minGSSize=5,readable=FALSE)
                    },
                    "david"={
                        res=enrichDAVID(gene=ann2$entrezId[i],idType="ENTREZ_GENE_ID",listType="Gene",annotation="KEGG_PATHWAY",david.user="clusterProfiler@hku.hk")
                    }
                    )
                    head(summary(res))
                    if (class(res) %in% c("try-error","NULL")) {
                        print(paste(ontThis,": Error",sep=""))
                    } else {
                        if (nrow(res@result)==0) {
                            print(paste(ontThis,": Nothing significant",sep=""))
                        } else {
                            fName3=paste(enrichFlag,"Anno",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                            write.table(res@result, file=paste(fName3,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
                            cat("PV<",signif(max(res@result$pvalue,na.rm=T),2),", QV<",signif(max(res@result$qvalue,na.rm=T),2),"\n")
                            fName3=paste(enrichFlag,"EnrichMap",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                            if (outFormat=="png") {
                                png(paste(fName3,".png",sep=""))
                            } else if (outFormat=="pdf") {
                                pdf(paste(fName3,".pdf",sep=""), pointsize=4)
                            }
                            try(enrichMap(res,vertex.label.font=.4))
                            dev.off()
                            fName3=paste(enrichFlag,"CgetPlot",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                            if (outFormat=="png") {
                                #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=12)
                                png(paste(fName3,".png",sep=""),width=1.5*480,height=1.5*480,pointsize=12)
                            } else if (outFormat=="pdf") {
                                pdf(paste(fName3,".pdf",sep=""))
                            }
                            try(cnetplot(res, categorySize="pvalue", foldChange=2^stat_1$coef[i0]))
                            dev.off()
                            if (enrichFlag=="go") {
                                fName3=paste(enrichFlag,"PlotGraph",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                                if (outFormat2=="png") {
                                    #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=12)
                                    #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=24)
                                    png(paste(fName3,".png",sep=""),width=5*480,height=5*480,pointsize=24)
                                } else if (outFormat2=="pdf") {
                                    pdf(paste(fName3,".pdf",sep=""))
                                }
                                try(plotGOgraph(res))
                                dev.off()
                            }
                        }
                    }
                    switch(enrichFlag,
                    "go"={
                    },
                    "kegg"={
                        res=enrichKEGG(gene=ann2$entrezId[i],organism=organismThis[2],universe=ann2$entrezId[i0],pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,use_internal_data=TRUE)
                        gseKEGG(geneList, organism = "hsa", keyType = "kegg", exponent = 1,
                        nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE,
                        seed = FALSE)                },
                    "do"={
                    },
                    "david"={
                    }
                    )
                    
                }
            }
        }
    }
}

###########################################################
###########################################################


#png("tmp.png",width=8*2*240,height=4*2*240)
#par(mfcol=c(2,4))
png("tmp.png",width=4*240,height=2*240)
par(mfcol=c(2,4))
#par(c(5, 4, 4, 2) + 0.1)
par(c(5, 8, 8, 2) + 0.1)
for (i in ii) {
    lim=NULL
    if (varType=="categorical") {
        boxplot(lcpmT[i1[i],samId]~phen10[samId,varList],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
        for (grpId in 1:length(grpUniq)) {
            j=which(phen10[samId,varList]==grpUniq[grpId])
            points(rep(grpId,length(j)),lcpmT[i1[i],samId[j]],col="black",cex=5)
            lines(grpId+c(-.4,.4),rep(mean(lcpmT[i1[i],samId[j]]),2),col="green")
        }
    } else {
        cexMain=cexLab=1.5
        cexMain=cexLab=2.5
        cexMain=2.5; cexLab=1.5
        #plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",cex.main=2.5,cex.axis=4,cex.lab=1.5)
        #plot(stat_1$coef,-log10(stat_1$pv),main=paste(subsetName,compName,"\nNo. with ",tolower(colNamePV)," < ",pThres,": ",length(i),sep=""),xlab="Log2 fold change",ylab="-Log10(p-value)",pch=19,cex=.8,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
        plot(phen10[samId,varList],lcpmT[i1[i],samId],ylim=lim,main=paste(header,"\n",ann10$geneId[i1[i]],", pv ",signif(stat_1$pv[i],2),sep=""),xlab=varList,ylab="log2(CPM)",pch=19,cex=.8,cex.main=cexMain,cex.lab=cexLab,cex.axis=1.5)
    }
}
dev.off()


###########################################################
###########################################################

datadir="results/rnaSeq/deGene/slopeCat/"
stat1=read.table(paste(datadir,x,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
