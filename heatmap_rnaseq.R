## Run deGenes_rnaseq.R section 1 then analysis_rnaseq.R section 1 first

## ----------------------------------------------
## ----------------------------------------------

capWords=function(s, strict=FALSE) {
    cap=function(s) paste(toupper(substring(s, 1, 1)),
    {s=substring(s, 2); if(strict) tolower(s) else s},
    sep="", collapse=" " )
    sapply(strsplit(s, split=" "), cap, USE.NAMES=!is.null(names(s)))
}

tolowerWords=function(s, strict=FALSE) {
    cap=function(s) paste(tolower(substring(s, 1, 1)),
    {s=substring(s, 2); if(strict) tolower(s) else s},
    sep="", collapse=" " )
    sapply(strsplit(s, split=" "), cap, USE.NAMES=!is.null(names(s)))
}


########################################################################
########################################################################
########################################################################

library(psy)
library(marray)
library(qvalue)
library(sva)
#source(paste(dirSrc,"functions/colorCluster.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.5.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))


datType="_stdzCoralHaoxu"
datType=""
datType="_stdzExpt"
datType="_combatAdj"

datList=""
datList=c("","_stdzCoralHaoxu","_stdzExpt","_combatAdj")

for (datType in datList) {

    candGeneFlag="_candGene2"
    candGeneFlag=""

    varType="_categorical"
    varType=""

    outFormat="pdf"
    outFormat="png"

    absFlag=T
    absFlag=F

    transformFlag=""

    datadirG="results/waterland/pearson/clusterInfo/"
    datadirG="results/waterland/spearman/clusterInfo/"


    if (F) {
        #candGene=candGene0[which(candGene0$tbl=="12864_2006_811_MOESM1_ESM.txt" & candGene0$signif==1 & candGene0$dirn%in%c(-1,1)),]
        candGene=candGene0[which(candGene0$signif==1 & candGene0$dirn%in%c(-1,1)),]
        candGene=candGene[!duplicated(toupper(candGene$geneSym)),]
    }

    phenAll=dgeT$samples[,c("id","tubeLabel","sampleId","experimentNo","dataset","treatment","treatment2","lib.size","norm.factors")]
    names(phenAll)[match(c("treatment2"),names(phenAll))]="group"
    annAll=ann10
    varPred=lcpmT

    if (F) {
        stat_1=statV_4_0
        nm=data.frame(var1=c("PValue"),var2=c("pv"),stringsAsFactors=F)
        k=match(nm$var1,names(stat_1))
        if (sum(!is.na(k))!=0) names(stat_1)[k[!is.na(k)]]=nm$var2[!is.na(k)]
        stat_1$bh=p.adjust(stat_1$pv,method="BH")
        stat_1$qv=NA
        i=which(!is.na(stat_1$pv))
        res=qvalue(stat_1$pv[i])
        stat_1$qv[i]=res$qvalue
        i=match(stat_1$geneId,annAll$geneId)
        #i=i[which(maxCntVec[i]>=quantile(maxCntVec[i],probs=.25) & sdCntVec[i]>=quantile(sdCntVec[i],probs=.25))]
    }

    distMethod="pearson"; absFlag=F
    distMethod="spearman"; absFlag=F
    distMethod="euclidean"; absFlag=F

    phen=phenAll

    annColAll=phen
    annRowAll=annAll

    linkMethod="ward.D2"

    colHmap=c("blue", "red", "white")
    colHmap=c("red", "blue", "white")
    colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","peachpuff","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
    colList=c("skyblue","blue","yellow","purple","black","red","orange","green","cyan","darkgreen")

    centerFlag=F
    centerFlag=T

    centerList=c(F,T)
    centerList=F
    centerList=T

    subsetList=c("","_male","_female")
    subsetList=c("")

    subsetFList=c("_top500","_rnd500")

    subsetF2List=c("_maxCnt10_sd0.2")
    subsetF2List=c("_maxCnt50_sd0.2")
    subsetF2List=c("_maxCnt500_sd0.2")
    subsetF2List=c("_maxCnt500_sd0.5")
    subsetF2List=c("_maxCntSdMedian")
    subsetF2List=c("")

    subsetFList=c("","_top2000"); subsetF2List=c("_maxCnt1_sd0.2","_maxCnt10_sd0.2")
    subsetFList=c(""); subsetF2List=c("_maxCnt1_sd0.2","_maxCnt10_sd0.2")
    subsetFList=c("","_top2000"); subsetF2List=c("_maxCnt10_sd0.2")
    subsetFList=c(""); subsetF2List=c("_maxCnt10_sd0.2")
    subsetFList=c(""); subsetF2List=c("_maxCnt1_sd0.2")
    subsetFList=c("_top2000"); subsetF2List=c("_maxCnt1_sd0.2")


    header1=""
    if (datType!="") {
        switch(datType,
        "_combatAdj"={header1="Coral/Haoxu ComBat adjusted: "
        },
        "_stdzExpt"={header1="Standardized within expt: "
        },
        "_stdzCoralHaoxu"={header1="Standardized within Coral/Haoxu: "
        },
        )
        prId=1:nrow(varPred)
        x=strsplit(subsetF2Flag,"_")[[1]]
        prId=prId[which(maxCntVec[prId]>=as.integer(sub("maxCnt","",x[2])) & sdCntVec[prId]>=as.numeric(sub("sd","",x[3])))]
        if (datType=="_combatAdj") {
            out=matrix(nrow=nrow(varPred),ncol=ncol(varPred),dimnames=list(rownames(varPred),colnames(varPred)))
            group=phen$group!="Sham"
            modCom = model.matrix(~as.factor(group))
            exprCom=ComBat(dat=varPred[prId,], batch=phen$experimentNo, mod=modCom, par.prior = TRUE,prior.plots = FALSE)
            out[prId,]=exprCom
        } else {
            k=c()
            if (datType=="_stdzCoralHaoxu") {
                grp=rep("Haoxu",nrow(phen))
                grp[which(phen$experimentNo=="A1402")]="Coral"
            } else {
                grp=phen$experimentNo
            }
            grpUniq=unique(grp)
            out=matrix(nrow=nrow(varPred),ncol=ncol(varPred),dimnames=list(rownames(varPred),colnames(varPred)))
            for (gId in 1:length(grpUniq)) {
                j=which(grp==grpUniq[gId])
                x=apply(varPred[,j],1,mad,na.rm=T)
                i=(x!=0)
                out[i,j]=(varPred[i,j]-apply(varPred[i,j],1,median,na.rm=T))/x[i]
                if (sum(!i)!=0) {
                    k=c(k,which(!i))
                }
            }
            k=unique(k)
            out[k,]=NA
            centerList=F
        }
        varPred=out
    }

    for (centerFlag in centerList) {
        for (subsetF2Flag in subsetF2List) {
            for (subsetFFlag in subsetFList) {
                for (subsetFlag in subsetList) {
                    #fName=paste("_p53_rnaSeq",datType,varType,"_",distMethod,ifelse(absFlag,"_abs",""),ifelse(centerFlag,"_center",""),sep="")
                    fName=paste("_p53_rnaSeq",datType,subsetFFlag,sep="")
                    if (candGeneFlag=="") {
                        clustFList=NA
                    } else {
                        fName2=paste("clusterInfoFeature",fName,sep="")
                        x=dir(datadirG,pattern=fName2)
                        if (length(x)==0) {
                            clustFList=NA
                        } else {
                            f1=read.table(paste(datadirG,fName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                            clustFList=unique(f1$clustId)
                        }
                    }
                    for (clustFFlag in clustFList) {
                        header=paste(header1,"Coral & Haoxu RNA-seq data",sep="")
                        #fName=paste("_p53_rnaSeq",datType,subsetFFlag,subsetF2Flag,subsetFlag,varType,"_",distMethod,ifelse(absFlag,"_abs",""),ifelse(centerFlag,"_center",""),sep="")
                        #fName=paste("_p53_rnaSeq",datType,sep="")
                        fName=paste("_p53_rnaSeq",datType,subsetFFlag,sep="")
                        cat("\n\n",fName,"\n\n")
                        prId=1:nrow(varPred)
                        prId=which(apply(varPred,1,function(x) {mean(!is.na(x))})!=0)
                        #prId=prId[1:100]
                        if (!is.na(clustFFlag)) {
                            header=paste("CpG ",clustFFlag,", ",header,sep="")
                            fName=paste("_tgfbTarget",datType,subsetFFlag,subsetF2Flag,subsetFlag,"_cpg",capWords(clustFFlag),varType,"_",distMethod,ifelse(absFlag,"_abs",""),ifelse(centerFlag,"_center",""),sep="")
                            prId=prId[match(f1$variable[which(f1$clustId==clustFFlag)],rownames(varPred)[prId])]
                        }
                        samId=1:ncol(varPred)
                        #samId=1:10
                        if (length(prId)<2 | length(samId)<2) next

                        varFList=varFName=NULL

                        x=strsplit(subsetF2Flag,"_")[[1]]
                        if (length(x)==2) {
                            if (subsetF2Flag=="_maxCntSdMedian") {
                                thres=0.5
                                header=paste(header,", min cnt / sd >= median",sep="")
                            }
                            prId=prId[which(maxCntVec[prId]>=quantile(maxCntVec[prId],probs=as.numeric(thres)) & sdCntVec[prId]>=quantile(sdCntVec[prId],probs=as.numeric(thres)))]
                            cat("Max cnt median ",quantile(maxCntVec[prId],probs=as.numeric(thres)),"\n",sep="")
                            cat("SD median ",round(quantile(sdCntVec[prId],probs=as.numeric(thres)),2),"\n",sep="")
                        } else {
                            prId=prId[which(maxCntVec[prId]>=as.integer(sub("maxCnt","",x[2])) & sdCntVec[prId]>=as.numeric(sub("sd","",x[3])))]
                            header=paste(header,", min cnt >= ",as.integer(sub("maxCnt","",x[2])),", sd >= ",as.numeric(sub("sd","",x[3])),sep="")
                        }
                        if (length(prId)==0) {
                            cat("No genes !!!\n")
                            next
                        }
                        if (substr(subsetFFlag,1,nchar("_top"))=="_top") {
                            x=as.integer(sub("_top","",subsetFFlag))
                            prId=prId[order(sdCntVec[prId],decreasing=T)][1:x]
                            header=paste(header,": ",x," top genes",sep="")
                        } else if (substr(subsetFFlag,1,nchar("_rnd"))=="_rnd") {
                            x=as.integer(sub("_rnd","",subsetFFlag))
                            set.seed(5393)
                            prId=sample(prId,size=x,replace=F)
                            header=paste(header,": random ",x," genes",sep="")
                        }
                        #}
                        if (length(prId)==0) {
                            cat("No genes !!!\n")
                            next
                        }

                        switch(subsetFlag,
                            "_male"={
                                samId=samId[which(phen$patient.gender[samId]=="male")]
                                header=paste(header,", male",sep="")
                            },
                            "_female"={
                                samId=samId[which(phen$patient.gender[samId]=="female")]
                                header=paste(header,", female",sep="")
                            }
                        )
                        if (length(samId)==0) {
                            cat("No samples !!!\n")
                            next
                        }
                        varList=varName=NULL
                        if (F) {
                            varList=c("experimentNo","treatment","treatment2")
                            varName=paste(varList," ",sep="")
                            varName=sub("treatment2","group",varName)
                        }
                        varList=c("experimentNo","group")
                        varName=paste(varList," ",sep="")

                        varListAll=varList; varNameAll=varName
                        varFListAll=varFList; varFNameAll=varFName

                        annCol=phen[samId,]
                        annRow=annAll[prId,]
                        arrayData=varPred[prId,samId]

                        if (centerFlag) {
                            arrayData=arrayData-apply(arrayData,1,median,na.rm=T)
                        }

                        cloneName=rownames(arrayData)
                        cloneName=rep("",length(prId))
                        if (is.null(varFList)) {
                            rowCol=NULL
                        } else {
                            rowCol=matrix(nrow=length(varFList),ncol=nrow(annRow))
                            for (varId in 1:length(varFList)) {
                                if (sum(!duplicated(annRowAll[!is.na(annRowAll[,varFList[varId]]),varFList[varId]]))>10) {
                                    x=round(annRowAll[,varFList[varId]])+1
                                    x2=as.integer(as.factor(x))
                                    j=match(annRow$id,annRowAll$id); j1=which(!is.na(j)); j2=j[j1]
                                    lim=range(x,na.rm=T)
                                    grpUniq=lim[1]:lim[2]
                                    rowColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    rowCol[varId,j1]=rowColUniq[x2[j2]]
                                } else {
                                    x=as.character(annRowAll[,varFList[varId]])
                                    x[x==""]=NA; x=as.integer(as.factor(x))
                                    grpUniq=sort(unique(x))
                                    x=x[match(annRow$id,annRowAll$id)]
                                    if (length(grpUniq)<=length(colList)) {
                                        rowCol[varId,]=colList[x]
                                    } else {
                                        rowCol[varId,]=rainbow(length(grpUniq))[x]
                                    }
                                }
                            }
                            rownames(rowCol)=varFName
                        }
                        lineClone=seq(5,nrow(arrayData),by=5)
                        lineClone=NULL

                        samName=rep("",length(samId))
                        if (length(samId)<=10) samName=annCol$id[samId]
                        samName=annCol$id[samId]
                        lineSam=NULL

                        if (is.null(varList)) {
                            colCol=NULL
                        } else {
                            colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                            for (varId in 1:length(varList)) {
                                if (sum(!duplicated(annColAll[!is.na(annColAll[,varList[varId]]),varList[varId]]))>10) {
                                    x=round(annColAll[,varList[varId]])+1
                                    x2=as.integer(as.factor(x))
                                    j=match(annCol$id,annColAll$id); j1=which(!is.na(j)); j2=j[j1]
                                    lim=range(x,na.rm=T)
                                    grpUniq=lim[1]:lim[2]
                                    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    colCol[varId,j1]=colColUniq[x2[j2]]
                                } else {
                                    x=as.character(annColAll[,varList[varId]])
                                    x[x==""]=NA; x=as.integer(as.factor(x))
                                    grpUniq=sort(unique(x))
                                    x=x[match(annCol$id,annColAll$id)]
                                    if (length(grpUniq)<=length(colList)) {
                                        colCol[varId,]=colList[x]
                                    } else {
                                        colCol[varId,]=rainbow(length(grpUniq))[x]
                                    }
                                }
                            }
                            rownames(colCol)=varName
                        }

                        switch(distMethod,
                            "pearson" = {
                                x=cor(arrayData,method=distMethod,use="complete.obs")
                                if (absFlag) x=abs(x)
                                distMat <- as.dist(1-x)
                                clustC <- hclust(distMat, method=linkMethod)
                                x=cor(t(arrayData),method=distMethod,use="complete.obs")
                                if (absFlag) x=abs(x)
                                distMat <- as.dist(1-x)
                                clustR <- hclust(distMat, method=linkMethod)},
                            "spearman" = {
                                x=cor(arrayData,method=distMethod,use="complete.obs")
                                if (absFlag) x=abs(x)
                                distMat <- as.dist(1-x)
                                clustC <- hclust(distMat, method=linkMethod)
                                x=cor(t(arrayData),method=distMethod,use="complete.obs")
                                if (absFlag) x=abs(x)
                                distMat <- as.dist(1-x)
                                clustR <- hclust(distMat, method=linkMethod)},
                            "kappa" = {
                                distMat <- getKappaDist(arrayData,absolute=absFlag)
                                clustC <- hclust(distMat, method=linkMethod)
                                distMat <- getKappaDist(t(arrayData),absolute=absFlag)
                                clustR <- hclust(distMat, method=linkMethod)},
                            "euclidean" = {
                                distMat <- dist(t(arrayData), method=distMethod)
                                clustC <- hclust(distMat, method=linkMethod)
                                distMat <- dist(arrayData, method=distMethod)
                                clustR <- hclust(distMat, method=linkMethod)}
                        )

                        #cloneName=rep("",length(prId))
                        #cloneName[clustR$order][lineClone]=paste("",lineClone,sep="")
                        
                        print("dim(arrayData)")
                        print(dim(arrayData))

                        summary(c(arrayData))
                        dat2=arrayData
                        #limit=quantile(c(arrayData),probs=c(0.1,.9),na.rm=T)
                        limit=c(-1,1)
                        #if (datType=="_stdzd") limit=c()

                        margins=c(4,.2)
                        margins=c(4,4)
                        margins=c(10,1)
                        main=header

                        #nClust=c(3,NA)
                        #nClust[2]=4
                        #if (subsetFlag=="_hispanic") nClust[2]=6
                        nClust=c(3,3)
                        while (T) {
                            x=table(cutree(clustR,k=nClust[1]))
                            if(max(x)<5) {
                                nClust[1]=NA
                                break
                            }
                            if (sum(x>=5)==3) break
                            nClust[1]=nClust[1]+1
                        }
                        while (T) {
                            x=table(cutree(clustC,k=nClust[2]))
                            if (sum(x>=5)==3) break
                            nClust[2]=nClust[2]+1
                        }

                        if (outFormat=="png") {
                            png(paste("heatmap",fName,".png",sep=""),width=480*2,height=480*2); cexRow=cexCol=2
                        } else {
                            pdf(paste("heatmap",fName,".pdf",sep="")); cexRow=cexCol=1
                        }
                        par(cex.main=.7)
                        #hcc <- heatmap3(x=dat2, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=samName, labRow=cloneName, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=2, cexRow=2, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
                        #hcc <- heatmap3(x=dat2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=samName, labRow=cloneName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=cexCol, cexRow=cexRow, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
                        hcc <- heatmap3(x=dat2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=samName, labRow=cloneName, lineRow=lineClone, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=cexCol, cexRow=cexRow, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
                        dev.off()

                        fName2=""

                        if (!is.null(colCol)) {
                            for (varId in 1:length(varListAll)) {
                                if (sum(!duplicated(annColAll[!is.na(annColAll[,varListAll[varId]]),varListAll[varId]]))>10) {
                                    if (outFormat=="png") {
                                        png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fName2,".png",sep=""),width=480,height=140)
                                    } else {
                                        pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fName2,".pdf",sep=""))
                                    }
                                    x=round(annColAll[,varListAll[varId]])+1
                                    lim=range(x,na.rm=T)
                                    if (varList[varId]==c("mutPerc")) lim=limPerc+1
                                    #if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                                    grpUniq=lim[1]:lim[2]
                                    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    lim=lim-1
                                    heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                } else {
                                    x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                                    x2=as.character(annCol[,varList[varId]]); x2[x2==""]=NA
                                    if (all(is.na(x2))) next
                                    x2=table(x2)
                                    x2=names(x2)
                                    grp=table(x)
                                    grpUniq=names(grp)
                                    k=1:length(grpUniq)
                                    ttl=grpUniq[k]
                                    ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                                    if (varList[varId]=="hpv") {
                                        ttl=sub("0","HPV-",ttl)
                                        ttl=sub("1","HPV+",ttl)
                                    }
                                    #if (clusterFlag[1]=="_supervised") {
                                    k=match(x2,grpUniq)
                                    #} else {
                                    #k=1:length(grpUniq)
                                    #}
                                    if (!is.null(varFListAll) && varFListAll[varId]%in%c("sd","mutPerc")) {
                                        width = 480; height = 140
                                    } else {
                                        if (length(grpUniq)<6) {
                                            width = 480; height = 480
                                        } else {
                                            width = 560; height = 960
                                        }
                                    }
                                    if (outFormat=="png") {
                                        png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fName2,".png",sep=""),width=width,height=height)
                                    } else {
                                        pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fName2,".pdf",sep=""))
                                    }
                                    cexThis=NULL
                                    if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                                    if (outFormat=="pdf") cexThis=1
                                    if (length(grpUniq)<=length(colList)) {
                                        sampleColorLegend(tls=ttl[k],col=colList[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else {
                                        sampleColorLegend(tls=ttl[k],col=rainbow(length(grpUniq))[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    }
                                }
                                dev.off()
                            }
                        }
                        if (!is.null(rowCol)) {
                            for (varId in 1:length(varFListAll)) {
                                if (sum(!duplicated(annRowAll[!is.na(annRowAll[,varFListAll[varId]]),varFListAll[varId]]))>10) {
                                    if (outFormat=="png") {
                                        png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".png",sep=""),width=480,height=140)
                                    } else {
                                        pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".pdf",sep=""))
                                    }
                                    x=round(annRowAll[,varFListAll[varId]])+1
                                    lim=range(x,na.rm=T)
                                    grpUniq=lim[1]:lim[2]
                                    rowColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    lim=lim-1
                                    heatmapColorBar(limit=lim,cols=c(rowColUniq[c(length(rowColUniq),1,median(1:length(rowColUniq)))]))
                                } else {
                                    x=as.character(annRowAll[,varFListAll[varId]]); x[x==""]=NA
                                    x2=as.character(annRow[,varFList[varId]]); x2[x2==""]=NA
                                    if (all(is.na(x2))) next
                                    x2=table(x2)
                                    x2=names(x2)
                                    grp=table(x)
                                    grpUniq=names(grp)
                                    k=1:length(grpUniq)
                                    ttl=grpUniq[k]
                                    ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                                    k=match(x2,grpUniq)
                                    if (!is.null(varFListAll) && varFListAll[varId]%in%c("sd","mutPerc")) {
                                        width = 480; height = 140
                                    } else {
                                        if (length(grpUniq)<6) {
                                            width = 480; height = 480
                                        } else {
                                            width = 560; height = 960
                                        }
                                    }
                                    if (outFormat=="png") {
                                        png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".png",sep=""),width=width,height=height)
                                    } else {
                                        pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],fName2,".pdf",sep=""))
                                    }
                                    cexThis=NULL
                                    if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                                    if (outFormat=="pdf") cexThis=1
                                    if (length(grpUniq)<=length(colList)) {
                                        sampleColorLegend(tls=ttl[k],col=colList[k],legendTitle=varFNameAll[varId],cex=cexThis)
                                    } else {
                                        sampleColorLegend(tls=ttl[k],col=rainbow(length(grpUniq))[k],legendTitle=varFNameAll[varId],cex=cexThis)
                                    }
                                }
                                dev.off()
                            }
                        }

                        if (!is.na(nClust[1])) {
                            if (F) {
                                #png(paste("clusterVariables",fName,".png",sep=""))
                                pdf(paste("clusterVariables",fName,".pdf",sep=""))
                                plot(clustR,main=paste("Variable clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
                                dev.off()
                            }

                            clustId=cutree(clustR,k=nClust[1])[clustR$order]
                            k1=which(!duplicated(clustId))
                            for (k in 1:length(k1)) {
                                clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                            }

                            #tbl=as.data.frame(as.matrix(arrayData[clustR$order,]),stringsAsFactors=F)
                            #tbl=data.frame(geneId=clustR$labels[clustR$order],clustId,order=1:nrow(arrayData),stringsAsFactors=F)
                            tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                            write.table(tbl, paste("clusterInfoFeature",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                        }

                        if (!is.na(nClust[2])) {
                            if (F) {
                            #png(paste("clusterSamples",fName,".png",sep=""))
                            pdf(paste("clusterSamples",fName,".pdf",sep=""))
                            plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                            dev.off()
                            }

                            clustId=cutree(clustC,k=nClust[2])[clustC$order]
                            k1=which(!duplicated(clustId))
                            for (k in 1:length(k1)) {
                            clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                            }

                            tbl=cbind(annCol[clustC$order,],clustId,order=1:nrow(annCol))
                            write.table(tbl, paste("clusterInfoSample",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                        }
                    }
                }
            }
        }
        if (outFormat=="png") {
            png(paste("heatmapColorRange_",subsetF2Flag,".png",sep=""))
        } else {
            pdf(paste("heatmapColorRange_",subsetF2Flag,".pdf",sep=""))
        }
        heatmapColorBar(limit=limit,cols=colHmap)
        dev.off()
    }
}
