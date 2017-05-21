
###########################################################
library(marray)
source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
#source(paste(dirSrc,"functions/heatmap.5.7.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))

load("results/clId.RData")

outFormat="pdf"
outFormat="png"

datadir=""

centrFlag="_noCentering"
centrFlag=""

subsetFlag=""

numPr=500
numPr=2000

pThres=10^-8
pThres=10^-6
pThres=0.05

compList=paste("_rnd",numPr,sep="")
compList=paste("_topVar",numPr,sep="")

datFlag="_combatAdj"
datFlag=""

colGeneId="probesetid"; colIdPV="FDR"; colNamePV="QV"

fName1="_moGene2.0"

tblCC=NULL
for (compFlag in compList) {
    if (length(grep("_rnd",compFlag))==1) {
        rndVec=paste("_rnd",1:4,sep="")
        #		rndVec=paste("_rnd",1:20,sep="")
    } else {
        rndVec=""
    }
    
    for (rndId in rndVec) {
        limFCmmu=c(-6,6)
        if (compFlag%in%c("_topSignif")) {
            load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
            switch(compFlag,
            "_topSignif"={
                stat_1=stat1_4
            }
            )
            timeThis=as.integer(sub("hrs","",strsplit(compFlag,"_")[[1]][2]))
            compName1=paste(timeThis,"hr: TGFbeta vs. untreated",sep="")
            cnt_1=eset$expr
            switch(datFlag,
                "_combatAdj"={cnt_1=exprCom
                }
            )
            phen_1=eset$phen
        } else {
            if (substr(compFlag,1,nchar("_topVar"))=="_topVar") {
                compName1=paste("Top ",sub("_topVar","",compFlag)," most variable probesets",sep="")
            }
            if (substr(compFlag,1,nchar("_rnd"))=="_rnd") {
                compName1=paste("Random ",sub("_rnd","",compFlag)," probesets",sep="")
            }
            cnt_1=eset$expr
            switch(datFlag,
            "_combatAdj"={cnt_1=exprCom
            }
            )
            ann_1=eset$ann
            phen_1=eset$phen
            stat_1=NULL
        }
        cnt_1=cnt_1[clId,]
        ann_1=ann_1[clId,]
        phen_1$id2=phen_1$id
        phen_1$id2=sub("Donor_","",sub("TGFbeta","tgfb",sub("hrs","h",sub("_hrs","hrs",phen_1$id))))
        for (transFlag in c("")) {
            if (F) {
                if (transFlag=="") {
                    subsetFlag=subsetName=""
                } else {
                    subsetFlag=paste("_",tolower(transFlag),sep="")
                    subsetName=paste(", ",transFlag,sep="")
                }
            }
            for (subsetFlag in c("","_sham","_gamma","_hze")) {
                compName2=compName1
                if (subsetFlag=="") {
                    subsetName=""
                    prId=NULL
                    samId=1:nrow(phen_1)
                    sampleBar="cluster"
                    geneBar="clusterPr"
                    nClust=c(2,3)
                } else {
                    grpUniq=sub("_","",subsetFlag)
                    subsetName=paste(", ",grpUniq,sep="")
                    if (grpUniq=="gamma") grpUniq="_"
                    fName2=paste(datadir,"clusterInfoFeature",fName1,compFlag,centrFlag,rndId,".txt",sep="")
                    fName2=paste(datadir,"clusterInfoFeature",fName1,centrFlag,rndId,".txt",sep="")
                    prId=read.table(file=fName2, header=T, sep="\t", quote="", comment.char="", as.is=T)
                    prId=prId[,"probesetid"]
                    samId=which(tolower(phen_1$group)==grpUniq)
                    sampleBar=""
                    sampleBar="cluster"
                    geneBar=""
                    nClust=c(NA,NA)
                }
                if (compFlag%in%c("_topSignif")) {
                    i1=which(stat_1[,colIdPV]<pThres)
                    if (length(i1)==0) next
                    fNameOut=paste(fName1,compFlag,subsetFlag,centrFlag,"_",colNamePV,pThres,datFlag,sep="")
                    header=paste(compName2,subsetName,", ",colNamePV,"<",pThres,sep="")
                    dat0=eset$expr
                    switch(datFlag,
                    "_combatAdj"={dat0=exprCom
                    }
                    )
                    dat0=dat0[match(stat_1[,colGeneId][i1],rownames(dat0)),]
                } else {
                    #fNameOut=paste(fName1,compFlag,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    fNameOut=paste(fName1,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    header=paste(compName2,subsetName,sep="")
                    dat0=cnt_1
                }
                expr=dat0[,samId]
                annRow=ann_1[match(rownames(expr),ann_1$probesetid),]
                phen=phen_1[samId,]
                
                i2=1:nrow(expr)
                if (length(grep("_rnd",compFlag))==1) {
                    if (compFlag==paste("_rnd",numPr,sep="")) {
                        i1=1:numPr
                    }
                    header=paste(header,": ",length(i1)," random probesets",sep="")
                    geneBar="clusterPr"
                    #set.seed(5453)
                    i2=sample(1:nrow(expr),length(i1),replace=F)
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("_topVar",compFlag))==1) {
                    if (is.null(prId)) {
                        varGene=apply(expr,1,var,na.rm=T)
                        i2=order(varGene,decreasing=T)[1:numPr]
                    } else {
                        compName2=paste(compName2," based on all samples",sep="")
                        i2=match(prId,annRow[,"probesetid"])
                    }
                    header=paste(compName2,subsetName,sep="")
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("_top",compFlag))==1) {
                    header=paste(header,", n=",nrow(expr),sep="")
                    geneBar=""
                    geneBar="clusterPr"
                    #i=order(annRow$logFC)
                    i=1:nrow(expr)
                } else {
                    geneBar=""
                    expr=expr[i2,]
                    annRow=cbind(annRow[i2,],logFC=stat_1$logFC[match(annRow[i2,colGeneId],stat_1[,colGeneId])])
                    i=order(annRow$logFC)
                }
                
                if (transFlag=="") {
                    j=1:ncol(expr)
                } else {
                    j=which(phen$translocation==transFlag)
                }
                
                
                arrayData=expr[i,j]
                annRow=annRow[i,]
                annCol=phen[j,]
                
                annColAll=eset$phen
                annColAll$id2=annColAll$id
                
                if (centrFlag=="") {
                    centr=apply(arrayData,1,median,na.rm=T)
                    arrayData=arrayData-centr
                }
                
                varList=c("treatment","group","experimentNo")
                varName=paste(varList," ",sep="")
                k=which(varList%in%names(annCol))
                varListAll=varList
                varNameAll=varName
                varList=varList[k]
                varName=varName[k]
                
                colList=c("skyblue","blue","yellow","purple","red")
                colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
                colList2=c("skyblue","blue")
                colHM=c("red","blue","grey")
                
                distMethod="pearson"
                linkMethod="ward.D2"
                
                #cloneName=annRow$geneSymbol
                cloneName=rep("",nrow(annRow))
                if (length(grep("_rnd|_topVar",compFlag))==1) {
                    cloneCol=NULL
                } else {
                    cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
                    k1=1; kk=which(names(annRow)=="logFC")
                    x=round(annRow[,kk]); x=x-min(x,na.rm=T)+1
                    grpUniq=sort(unique(x[!is.na(x)]))
                    x=round(annRow[,kk]); x[x<limFCmmu[1]]=limFCmmu[1]; x[x>limFCmmu[2]]=limFCmmu[2]; x=x+limFCmmu[2]+1
                    grpUniq=limFCmmu[1]:limFCmmu[2]
                    cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                    cloneCol[k1,]=cloneColUniq[x]
                    rownames(cloneCol)="log2FC "
                }
                
                if (F) {
                    if (subsetFlag=="") {
                        samName=rep("",ncol(arrayData))
                    } else {
                        samName=annCol$id2
                    }
                }
                samName=annCol$id2
                samCol=NULL
                samCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("lib.size")) {
                        j=match(annCol$id,annColAll$id)
                        x=round(annColAll[,varList[varId]])
                        lim=range(x,na.rm=T)
                        #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        grpUniq=lim[1]:lim[2]
                        samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        samCol[varId,]=samColUniq[x[j]]
                    } else {
                        if (varList[varId]%in%c("time")) {
                            x=annColAll[,varList[varId]]
                        } else {
                            x=as.character(annColAll[,varList[varId]])
                        }
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annCol$id,annColAll$id)]
                        if (length(grpUniq)<=length(colList2)) {
                            samCol[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            samCol[varId,]=colList[x]
                        } else {
                            samCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(samCol)=varName
                
                print("summary(range(c(arrayData),na.rm=T))")
                print(summary(range(c(arrayData),na.rm=T)))
                if (centrFlag=="") {
                    limit=c(-120000,120000)
                    limit=c(-10000,10000)
                    limit=c(-8,8)
                    limit=c(-1,1)
                    limit=c(-3,3)
                } else {
                    limit=c(8,13)
                }
                main=NULL
                main=header
                
                switch(distMethod,
                "pearson"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    if (sampleBar=="cluster") {
                        clustC=hclust(distMat, method=linkMethod)
                    } else {
                        clustC=NA
                        nClust[2]=NA
                    }
                    if (geneBar=="clusterPr") {
                        distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                        clustR=hclust(distMat, method=linkMethod)
                    } else {
                        clustR=NA
                        nClust[1]=NA
                    }
                },
                "spearman"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                },
                "euclidean"={distMat=dist(t(arrayData), method=distMethod)
                    distMat=dist(arrayData, method=distMethod)
                }
                )
                
                if (F) {
                    subDir <- paste(compFlag,sep="")
                    if (!file.exists(subDir)){
                        dir.create(file.path(subDir))
                    }
                    subDir=paste(subDir,"/",sep="")
                }
                subDir=""
                if (outFormat=="png") {
                    margins=c(6,1)
                    margins=c(10,20)
                    png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
                } else {
                    margins=c(12,5)
                    pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
                }
                totalC=nrow(eset$phen)
                totalC=ncol(arrayData)
                #hcc=heatmap3(x=arrayData, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
                hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3],totalC=totalC)
                dev.off()
                
                if (!is.na(nClust[1])) {
                    if (F) {
                        #png(paste("clusterVariables",fNameOut,".png",sep=""))
                        pdf(paste("clusterVariables",fNameOut,".pdf",sep=""))
                        plot(clustR,main=paste("Variable clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
                        dev.off()
                    }
                    
                    clustId=cutree(clustR,k=nClust[1])[clustR$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    
                    #tbl=as.data.frame(as.matrix(arrayData[clustR$order,]),stringsAsFactors=F)
                    tbl=data.frame(variable=clustR$labels[clustR$order],clustId,order=1:nrow(arrayData),stringsAsFactors=F)
                    tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                    write.table(tbl, paste("clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
                
                if (is.na(nClust[2])) {
                    clustId=paste("cluster",1,sep="")
                    tbl=cbind(annCol[clustC$order,which(names(annCol)%in%c("order"))],clustId,order=1:nrow(annCol))
                    write.table(tbl, paste("clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                } else {
                    if (F) {
                        #png(paste("clusterSamples",fNameOut,".png",sep=""))
                        pdf(paste("clusterSamples",fNameOut,".pdf",sep=""))
                        plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                        dev.off()
                    }
                    
                    clustId=cutree(clustC,k=nClust[2])[clustC$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    
                    tbl=cbind(annCol[clustC$order,which(!names(annCol)%in%c("order"))],clustId,order=1:nrow(annCol))
                    write.table(tbl, paste("clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
            }
        }
    }
}
if (!is.null(samCol)) {
    for (varId in 1:length(varListAll)) {
        if (outFormat=="png") {
            png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""))
        } else {
            pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
        }
        if (varListAll[varId]%in%c("age","wbc")) {
            x=round(annColAll[,varListAll[varId]])
            lim=range(x,na.rm=T)
            #lim=quantile(x,probs=c(.1,.9),na.rm=T)
            grpUniq=lim[1]:lim[2]
            samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
            heatmapColorBar(limit=lim,cols=c(samColUniq[c(length(samColUniq),1)],median(samColUniq)))
        } else {
            if (varList[varId]%in%c("time")) {
                x=annColAll[,varListAll[varId]]
            } else {
                x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
            }
            grpUniq=table(x)
            #		grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
            grpUniq=names(grpUniq)
            k=1:length(grpUniq)
            if (length(grpUniq)<=length(colList2)) {
                sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varNameAll[varId])
            } else if (length(grpUniq)<=length(colList)) {
                sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId])
            } else {
                sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId])
            }
        }
        dev.off()
    }
}
if (F) {
    heatmapColorBar=function(limit,cols=c("green","red","black"),main=NULL) {
        try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
        maColorBar(try, scale=limit,k=3,main=main)
    }
}
if (outFormat=="png") {
    png(paste("heatmapColorRange.png",sep=""),width=480,height=140)
} else {
    pdf(paste("heatmapColorRange.pdf",sep=""))
}
heatmapColorBar(limit=limit,cols=colHM,main="Heatmap color range")
dev.off()

###########################################################
###########################################################


###########################################################
###########################################################
