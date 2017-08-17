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

## ----------------------------------------------
## ----------------------------------------------


##############################################
## Section 1

verbose=T
verbose=F

dataset <- "p53"
dataset <- "tmp"

cohort="mmu"
organism <- "Mus_musculus"

if (computerFlag=="cluster") {
    R.utils::use("aroma.seq, edgeR, aroma.light, matrixStats")
    
    if (F) {
        datadir="htseqCountData/tmp,tophat2,gtf/Mus_musculus/"
        fName="coral_1CR_1FR.count"
        tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=5)
        
        datadir="htseqCountData/tmp,tophat2,gtf/Homo_sapiens/"
        fName="coral_1CR_1FR.count"
        tbl2=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=5)
    }
    
    counts <- HTSeqCountDataSet$byName(dataset, tags = "tophat2,gtf", organism = organism)
    counts <- setFullNamesTranslator(counts, function(name, ...) chartr(";", "/", name))
    counts
    data <- extractMatrix(counts, column = 2L, colClass = "integer")
    dge <- extractDGEList(counts)
    
    if (F) {
        db <- TabularTextFile("data/sampleInfo.txt")
        samples <- readDataFrame(db)
        samples=samples[match(colnames(data),samples$id),]
    }
    
    samples=data.frame(id=colnames(dge),dataset=sapply(colnames(dge),function(x) {
        y=strsplit(x,"_")[[1]]
        if (length(y)==2) {
            n=2
        } else {
            n=length(y)-2
        }
        y=paste(y[1:n],collapse="_")
    }),stringsAsFactors=F)
    rownames(samples)=NULL
} else {
    library("edgeR")
    library("matrixStats")
    
    datadir="data/count/"
    fileList=dir(datadir)
    tbl=read.table(paste(datadir,fileList[1],sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    nm=sub("_01","_1",sub(".count","",fileList,fixed=T))
    counts=matrix(nrow=nrow(tbl),ncol=length(fileList),dimnames=list(tbl[,1],nm))
    for (fId in 1:length(fileList)) {
        fName=fileList[fId]
        tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        if (all(rownames(counts)==tbl[,1])) {
            counts[,fId]=tbl[,2]
        }
    }
    phen=as.data.frame(t(sapply(fileList,function(x) {
        id=sub("_01","_1",sub(".count","",x,fixed=T))
        experimentNo=strsplit(id,"_")[[1]][1]
        c(id,experimentNo)
    },USE.NAMES=F)),stringsAsFactors=F)
    names(phen)=c("id","experimentNo")
    datadir="docs/"
    if (F) {
        phen1=read.table(paste(datadir,"RNA sequencing sample ID with treatments.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        phen2=read.table(paste(datadir,"RNA sequencing_Coral Omene_Barcellos-Hoff.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(phen2)[match(c("Labels.on.the.tube","Sample.ID","User.name","Date.and.Time","Nucleic.Acid","Unit","A260..Abs.","A280..Abs.","X260.280","X260.230","Volume"),names(phen2))]=c("id","sampleId","userName","dateTime","nucleicAcid","unit","A260Abs","A280Abs","ratio260_280","ratio260_230","volume")
        phen2$sampleId=gsub("-","_",gsub(" ","",phen2$sampleId))
        j=match(phen$id,phen2$sampleId); j1=which(!is.na(j)); j2=j[j1]
    }
    phen2=read.table(paste(datadir,"Lin-RNA sequencing sample IDs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(phen2)[match(c("dataset","Labels.on.the.tube","Sample.ID","Experiment..","Treatments","Groups"),names(phen2))]=c("dataset","id","sampleId","experimentNo","treatment","treatment2")
    phen2$id=gsub("-","_",gsub(" ","",phen2$id))
    phen2$sampleId=gsub("-","_",gsub(" ","",phen2$sampleId))
    phen2$treatment[which(phen2$treatment=="NTE-Sham")]="NTE_Sham"
    phen2$treatment2[which(phen2$treatment2=="sham")]="Sham"
    phen2$tubeLabel=phen2$id
    j=which(phen2$dataset=="Coral")
    phen2$id[j]=phen2$sampleId[j]
    j=match(phen$id,phen2$id); j1=which(!is.na(j)); j2=j[j1]
    #phen=phen[match(colnames(counts),phen$id),]
    phen=cbind(phen[match(colnames(counts),phen$id),],phen2[j2,which(!names(phen2)%in%names(phen))])
    phen2=read.table(paste(datadir,"Suervised cluster RNAseq and Microarray data - rnaSeq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
    names(phen2)[match(c("Sample.ID","Fast.Slopes","X","Sample.ID.1","Slow.Slopes"),names(phen2))]=c("id","slope","X","id2","slope2")
    phen3=phen2
    phen3$id=phen3$id2
    phen3$slope=phen3$slope2
    phen3$slopeCat="slow"
    phen2$slopeCat="fast"
    phen2=rbind(phen2,phen3)
    phen2$id=gsub(" ","",phen2$id)
    j=match(phen$sampleId,phen2$id); j1=which(!is.na(j)); j2=j[j1]
    phen$slope=NA
    phen$slopeCat=""
    phen[j1,c("slope","slopeCat")]=phen2[j2,c("slope","slopeCat")]
    rm(phen2,phen3)
    dge=DGEList(counts=counts,samples=phen,group=NULL,genes=NULL)
    samples=dge$samples
}
dge$samples$investigator=rep("Haoxu",nrow(dge$samples))
dge$samples$investigator[which(dge$samples$experimentNo=="A1402")]="Coral"


samples1=samples

data2=dge$counts
rho <- cor(data2, method = "spearman")
round(summary(as.vector(rho)),2)
j=2
round(rho[j,j],2)

dge$samples[1:3, -1]
samples=cbind(samples,dge$samples[,c("lib.size","norm.factors")])
rownames(samples)=NULL
if (verbose) {
    write.table(cbind(geneId=rownames(dge$counts),dge$counts), file=paste("count_raw_",organism,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    write.table(samples, file=paste("sample_",organism,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
}

dge$samples$group <- samples$experimentNo

dge1=dge
dge=dge1


kruskal.test(samples$lib.size~as.factor(samples$experimentNo))
#p-value = 2.854e-06

fName2="_p53"

dge <- calcNormFactors(dge, method = "TMM")
dgeT=dge


cpmT=cpm(dgeT,log=F)
lcpmT=cpm(dgeT,log=TRUE)
dgeT$samples$group=dgeT$samples$experimentNo

maxCntVec0=apply(dgeT$counts,1,max,na.rm=T)

maxCntVec=apply(dgeT$counts,1,max,na.rm=T)
png("maxCntPlot_raw.png")
hist(maxCntVec)
dev.off()
summary(maxCntVec)
maxCntVecRaw=maxCntVec

maxCntVec=apply(cpmT,1,max,na.rm=T)
png("maxCntPlot_cpm.png")
hist(maxCntVec)
dev.off()
summary(maxCntVec)
maxCntVecCpm=maxCntVec

sdCntVec=apply(lcpmT,1,sd,na.rm=T)
png("sdCntPlot.png")
hist(sdCntVec)
dev.off()
summary(sdCntVec)

alcpmT=aveLogCPM(dgeT, normalized.lib.sizes=TRUE)


if (F) {
    dge=dge1
    samples=samples1[match(rownames(dge$samples),samples1$id),]
    
    colSam=rainbow(ncol(dge$counts))
    
    dge <- calcNormFactors(dge, method = "TMM")
    
    png("densityPlot_normCount.png")
    par(mfcol=c(3,3))
    xlim=range(c(dge$counts),na.rm=T)
    xlim=c(0,50)
    j=1
    plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=paste("Densities of scale-normalized log2 gene counts across the ",ncol(dge$counts)," samples",sep=""),xlab="log(count)",col=colSam[j])
    for (j in 2:ncol(dge$counts)) {
        #	lines(density(dge$counts[,j],na.rm=T),col=colSam[j])
        plot(density(dge$counts[,j],na.rm=T),col=colSam[j],xlim=xlim)
    }
    dev.off()
    
    dge$samples$group <- samples$experimentNo
    
    dgeT <- dge[, !is.na(dge$samples$group)]
    
    dgeT <- estimateDisp(dgeT, trend = "none", robust = TRUE)
    print("table(is.infinite(dgeT$prior.df))")
    print(table(is.infinite(dgeT$prior.df)))
    "
    FALSE  TRUE
    30728 34488
    "
    
    save(dgeT,file="dgeT_tmp.RData")
}


## -----------------------------------------
## Section 2

subsetFlag=""

fName1=paste("_",organism,subsetFlag,sep="")
dge=dge1
samples=samples1[match(rownames(dge$samples),samples1$id),]

j=1:nrow(dge$samples)
subsetName=subsetFlag
dge=dge[,j]

colSam=rainbow(ncol(dge$counts))

dge <- calcNormFactors(dge, method = "TMM")
samples=samples1[match(rownames(dge$samples),samples1$id),]
samples=cbind(samples,dge$samples[,c("lib.size","norm.factors")])

datadir="results/comparison/"
datadir=""
ann <- read.table(paste(datadir,"ann_",organism,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
ann$length=sapply(ann$chromosome,function(x) {
    y=strsplit(x,":")[[1]][2]
    y=strsplit(y," ")[[1]][1]
    y=diff(as.integer(strsplit(y,"-")[[1]]))
    y
},USE.NAMES=F)
i=match(rownames(dge$counts),sub("*","",ann$gene_id,fixed=T))
x=rpkm(dge, gene.length=ann$length[i], normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
if (verbose) {
    write.table(cbind(geneId=rownames(dge$counts),x), file=paste("rpkm",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    #save(dge,file=paste("dge",fName1,".RData",sep=""))
    write.table(cbind(geneId=rownames(dge$counts),dge$counts), file=paste("count",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    write.table(samples, file=paste("sample",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
}

if (F) {
    png("densityPlot_normCount.png",width=3*240, height=3*240)
    par(mfcol=c(3,3))
    xlim=range(c(dge$counts),na.rm=T)
    xlim=c(0,50)
    xlim=NULL
    j=1
    plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=paste("Densities of scale-normalized log2 gene counts across the ",ncol(dge$counts)," samples",sep=""),xlab="log(count)",col=colSam[j])
    for (j in 2:ncol(dge$counts)) {
        lines(density(dge$counts[,j],na.rm=T),col=colSam[j])
        #	plot(density(dge$counts[,j],na.rm=T),col=colSam[j],xlim=xlim)
    }
    dev.off()
}

if (F) {
    png(paste("densityPlot_normCount",subsetName,".png",sep=""),width=3*240, height=3*240)
    par(mfcol=c(3,3))
    for (j in 1:ncol(dge$counts)) {
        xlim=c(0,quantile(dge$counts[,j],probs=.9,na.rm=T))
        plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=rownames(dge$samples)[j],col=colSam[j])
    }
    dev.off()
    
    png(paste("histogram_normCount",subsetName,".png",sep=""),width=3*240, height=3*240)
    par(mfcol=c(3,3))
    for (j in 1:ncol(dge$counts)) {
        xlim=c(0,quantile(dge$counts[,j],probs=.9,na.rm=T))
        xlim=c(0,quantile(dge$counts[,j],probs=.75,na.rm=T))
        hist(dge$counts[,j],xlim=xlim,main=rownames(dge$samples)[j])
    }
    dev.off()
}

dge$samples$group <- samples$experimentNo

grpUniq=unique(dge$samples$group)
grpName=grpUniq

## -----------------------------------------
## Section 3

## Voom
## Run section 1 first

library(limma)
library(edgeR)
library(sva)

fName2="_rnaSeq_p53"

dgeT$samples$slopeCat2=as.integer(dgeT$samples$slope>median(dgeT$samples$slope,na.rm=T))
dgeT$samples$slopeCat3=as.integer(dgeT$samples$slopeCat=="fast")
dgeT$samples$slopeCat3[!dgeT$samples$slopeCat%in%c("slow","fast")]=NA
dgeT$samples$treatment2[which(dgeT$samples$treatment2=="Sham")]="0Sham"
dgeT$samples$treatExpt=paste(dgeT$samples$treatment2,dgeT$samples$experimentNo)

## Will try to adjust for covariate adjFlag or a succeeding variable in varInfo
## Implemented for slope only !!!
varInfo=data.frame(formula=c("treatment2+experimentNo","treatment2","experimentNo"),variable=c("treatExpt","treatment2","experimentNo"),name=c("_treatAndExptNoAdj","_treatAdj","_exptNoAdj"),stringsAsFactors=F)
adjFlag="_treatAndExptNoAdj"

varInfo=data.frame(formula=c("treatExpt","treatment2","experimentNo"),variable=c("treatExpt","treatment2","experimentNo"),name=c("_treatExptNoAdj","_treatAdj","_exptNoAdj"),stringsAsFactors=F)
adjFlag="_exptNoAdj"
adjFlag=""
adjFlag="_treatAdj"
adjFlag="_treatExptNoAdj"

varInfo=data.frame(formula=c("treatment2"),variable=c("treatment2"),name=c("_treatAdj"),stringsAsFactors=F)
adjFlag=""

varInfo=NULL
adjFlag=""

compList=paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe"),"VsSham",sep=""); subsetList=""
compList="_A1026VsA1014"; subsetList=""
compList="_gammaVsSham"; subsetList="_A1014A1402"
compList="_capeGammaVsGamma"; subsetList=""
compList="_capeGammaVsGamma"; subsetList="_A1402"

compList="_slope"
subsetList=c("",paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe","sham"),sep=""))

compList="_slopeCat"
subsetList=c("",paste("_",c("cape","capeGamma","gamma","gammaTe","hze","hzeTe","sham"),sep=""))
subsetList=""

for (subsetFlag in subsetList) {
    for (compFlag in compList) {
        cat("\n\n=========== ",subsetFlag,", ",compFlag," ==============\n",sep="")
        varType="categorical"
        switch(compFlag,
            "_A1026VsA1014"={
                varList="experimentNo"; grpUniq=c("A1014","A1026"); grpName=grpUniq
            },
            "_capeVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","CAPE"); grpName=c("sham","cape")
            },
            "_capeGammaVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","CAPE+ gamma-rays"); grpName=c("sham","capeGamma")
            },
            "_gammaVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","gamma-rays"); grpName=c("sham","gamma")
            },
            "_gammaTeVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","gamma-rays-TE"); grpName=c("sham","gammaTe")
            },
            "_hzeVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","HZE"); grpName=c("sham","hze")
            },
            "_hzeTeVsSham"={
                varList="treatment2"; grpUniq=c("0Sham","HZE-TE"); grpName=c("sham","hzeTe")
            },
            "_capeGammaVsGamma"={
                varList="treatment2"; grpUniq=c("gamma-rays","CAPE+ gamma-rays"); grpName=c("gamma","capeGamma")
            },
            "_slope"={
                varType="continuous"
                varList="slope"; grpUniq=grpName=NULL
            },
            "_slopeCat"={
                varList="slopeCat2"; grpUniq=c(0,1); grpName=c("slow","fast")
                varList="slopeCat3"; grpUniq=c(0,1); grpName=c("slow","fast")
            }
        )

        #subsetFlag=""; varList="experimentNo"; grpUniq=c("A1014","A1026"); grpName=grpUniq
        #subsetFlag="_A1014A1402"; varList="treatment2"; grpUniq=c("0Sham","gamma-rays"); grpName=c("sham","gamma")


        maxCntVec=maxCntVecCpm
        samId=1:nrow(dgeT$samples)
        switch(subsetFlag,
            "_A1402"={
                samId=which(dgeT$samples$experimentNo%in%c("A1402"))
            },
            "_A1014A1402"={
                samId=which(dgeT$samples$experimentNo%in%c("A1014","A1402"))
            },
            "_cape"={
                samId=which(dgeT$samples$treatment2%in%c("CAPE"))
            },
            "_capeGamma"={
                samId=which(dgeT$samples$treatment2%in%c("CAPE+ gamma-rays"))
            },
            "_gamma"={
                samId=which(dgeT$samples$treatment2%in%c("gamma-rays"))
            },
            "_gammaTe"={
                samId=which(dgeT$samples$treatment2%in%c("gamma-rays-TE"))
            },
            "_hze"={
                samId=which(dgeT$samples$treatment2%in%c("HZE"))
            },
            "_hzeTe"={
                samId=which(dgeT$samples$treatment2%in%c("HZE-TE"))
            },
            "_sham"={
                samId=which(dgeT$samples$treatment2%in%c("0Sham"))
            }
        )
        if (is.null(grpUniq)) {
            samId=samId[!is.na(dgeT$samples[samId,varList])]
        } else {
            samId=samId[which(dgeT$samples[samId,varList]%in%grpUniq)]
        }
        if (varType=="categorical") {
            if (sum(!duplicated(dgeT$samples[samId,varList]))<2) {
                cat("Not Run. Less than 2 categories !!!\n")
                next
            }
        }
        
        i=1:nrow(dgeT$counts)
        dgeF=dgeT[i,samId]
        #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
        #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
        if (is.null(varInfo)) {
            covFlag=""
            varList2=NULL
        } else {
            var2Info=varInfo
            if (adjFlag%in%c("_treatExptNoAdj","_treatAndExptNoAdj") & sum(!duplicated(dgeF$samples$treatment2))==1) {
                adj2Flag="_exptNoAdj"
            } else {
                adj2Flag=adjFlag
            }
            if (varList=="slope") {
                k=which(var2Info$name==adj2Flag)
                if (length(k)==1 && k>1) var2Info=var2Info[-(1:(length(k)-1)),]
                covFlag=""
                varList2=NULL
                for (k in 1:nrow(var2Info)) {
                    x=table(dgeF$samples[!is.na(dgeF$samples[,varList]),var2Info$variable[k]])
                    x=x[which(x!=1)]
                    if (length(x)>1) {
                        varList2=strsplit(var2Info$formula[k],"+",fixed=T)[[1]]
                        covFlag=var2Info$name[k]
                        break
                    }
                }
                if (F) {
                    varList2="treatment2"
                    x=table(dgeF$samples[,varList2])
                    if (length(x)>1) {
                        covFlag="_treatAdj"
                    } else {
                        varList2="experimentNo"
                        x=table(dgeF$samples[,varList2])
                        if (length(x)>1) {
                            covFlag="_exptNoAdj"
                        } else {
                            covFlag=""
                            varList2=NULL
                        }
                    }
                }
            } else {
                k=which(var2Info$name==adj2Flag)
                if (length(k)==1 && k>1) var2Info=var2Info[-(1:(length(k)-1)),]
                covFlag=""
                varList2=NULL
                for (k in 1:nrow(var2Info)) {
                    x=table(dgeF$samples[,varList],dgeF$samples[,var2Info$variable[k]])
                    y=x!=0
                    if (nrow(x)>1 && ncol(x)>1 && any(y[1,]==y[2,])) {
                        varList2=strsplit(var2Info$formula[k],"+",fixed=T)[[1]]
                        covFlag=var2Info$name[k]
                        break
                    }
                }
                if (F) {
                    varList2="experimentNo"
                    x=table(dgeF$samples[,varList],dgeF$samples[,varList2])
                    y=x!=0
                    if (nrow(x)>1 & ncol(x)>1 & any(y[1,]==y[2,])) {
                        covFlag="_exptNoAdj"
                    } else {
                        varList2="investigator"
                        x=table(dgeF$samples[,varList],dgeF$samples[,varList2])
                        y=x!=0
                        if (nrow(x)>1 & ncol(x)>1 & any(y[1,]==y[2,])) {
                            covFlag="_investAdj"
                        } else {
                            covFlag=""
                            varList2=NULL
                        }
                    }
                }
            }
            }
        
        minCnt=10
        minCnt=1
        for (minCnt in c(0,1,2,5,10,20,50)) {
            #for (minCnt in c(1)) {
            cat("\n------------- minCnt ",minCnt," --------------\n",sep="")
            fName3=paste("_",ifelse(is.null(grpUniq),varList,paste(grpName[2],"Vs",capWords(grpName[1]),sep="")),covFlag,subsetFlag,fName2,"_voom_minCnt",minCnt,sep="")
            i=which(maxCntVec>=minCnt)
            #dgeF=list(counts=dgeT$counts[i,],genes=dgeT$genes[i,],samples=dgeT$samples)
            dgeF=dgeT[i,samId]
            #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
            #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
            ##dgeF$samples$group=as.factor(as.character(dgeF$samples[,varList]))
            if (covFlag=="") {
                modelThis=paste("~",varList,sep="")
            } else {
                #modelThis=paste("~",varList,"+as.factor(",varList2,")",sep="")
                modelThis=paste("~",varList,paste("+as.factor(",varList2,")",collapse=""),sep="")
            }
            if (varType=="categorical") {
                modelThis=sub("varList","as.factor(varList)",modelThis)
            }
            modelThis=as.formula(modelThis)
            design=model.matrix(modelThis,data=dgeF$samples)
            #design=model.matrix(~group,data=dgeF$samples)
            #cat(colnames(design),"\n",sep=", ")
            if (F) {
                #fit <- voom(dgeF$counts,design,plot=TRUE)
                #fit <- voom(dgeF$counts,design,save.plot=TRUE)
                png("tmp.png")
                plot(fit)
                dev.off()
            }
            #dat <- voom(dgeF$counts,design,save.plot=F)
            dat <- voom(dgeF,design,save.plot=F)
            #save(dat,design,file=paste("voom",fName3,".RData",sep=""))
            fit <- lmFit(dat,design)
            rm(dat)
            fit <- eBayes(fit)
            #topTable(fit,coef=ncol(design))
            colId=2
            top=cbind(geneId=rownames(fit$coef),logFC=fit$coef[,colId],PValue=fit$p.value[,colId])
            write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }

        ## -----------------------------------------
        ## Wilcox test

        if (varType=="categorical") {
            i=1:nrow(dgeT)
            dgeF=dgeT[i,samId]
            #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
            #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
            dat=lcpmT[i,samId]
            if (varType=="categorical") {
                dgeF$samples$group=as.factor(as.character(dgeF$samples[,varList]))
            } else {
                dgeF$samples$group=dgeF$samples[,varList]
            }
            grp=dgeF$samples$group
            fName3=paste("_",ifelse(is.null(grpUniq),varList,paste(grpName[2],"Vs",capWords(grpName[1]),sep="")),subsetFlag,fName2,"_wilcox",sep="")

            #save(dat,grp,fName3,file="tmp.RData")

            library(coin)
            #load("tmp.RData")
            distrib="exact"
            timeStamp=Sys.time()
            print(format(timeStamp, "%x %X"))
            tmp=rep(NA,nrow(dat))
            top=data.frame(geneId=rownames(dat),logFC=tmp,PValue=rep(NA,nrow(dat)),stringsAsFactors=F)
            i1=1:100
            i1=1:nrow(dat)
            print(nrow(dat))
            for (i in i1) {
                if (i%%1000==0) print(i)
                res=try(wilcox_test(dat[i,]~grp,distribution=distrib))
                x=try(pvalue(res))
                if (class(x)=="try-error") next
                top$PValue[i]=x
            }
            timeStamp=c(timeStamp,Sys.time())
            print(format(timeStamp[2], "%x %X"))
            print(diff(timeStamp))
            write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }

        ## -----------------------------------------
        ## SVA based

        #if (!is.null(varList2)) {
        if (T) {
            minCnt=1
            minCnt=50
            minCnt=500
            maxCntVec=maxCntVecRaw

            minCnt=10
            minCnt=1
            maxCntVec=maxCntVecCpm

            for (minCnt in c(0,1,2,5,10,20,50)) {
            #for (minCnt in c(1)) {
                cat("\n------------- SVA: minCnt ",minCnt," --------------\n",sep="")
                i=which(maxCntVec>=minCnt/10^6)
                i=which(maxCntVec>=minCnt)

                dgeF=dgeT[i,samId]
                #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
                #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
                ##dgeF$samples$group=as.factor(as.character(dgeF$samples$treatment2))
                if (varType=="categorical") {
                    dgeF$samples$group=as.factor(as.character(dgeF$samples[,varList]))
                } else {
                    dgeF$samples$group=dgeF$samples[,varList]
                }
                fName3=paste("_",ifelse(is.null(grpUniq),varList,paste(grpName[2],"Vs",capWords(grpName[1]),sep="")),covFlag,subsetFlag,fName2,"_voom_sva_minCnt",minCnt,sep="")
                if (covFlag=="") {
                    modelThis=paste("~",varList,sep="")
                } else {
                    #modelThis=paste("~",varList,"+as.factor(",varList2,")",sep="")
                    modelThis=paste("~",varList,paste("+as.factor(",varList2,")",collapse=""),sep="")
                }
                if (varType=="categorical") {
                    modelThis=sub("varList","as.factor(varList)",modelThis)
                }
                modelThis=as.formula(modelThis)
                mod = model.matrix(modelThis, data=dgeF$samples)
                if (covFlag=="") {
                    modelThis=as.formula(paste("~1",sep=""))
                } else {
                    #modelThis=as.formula(paste("~as.factor(",varList2,")",sep=""))
                    modelThis=as.formula(paste("~",paste("as.factor(",varList2,")",collapse="+"),sep=""))
                }
                #modelThis=as.formula(paste("~as.factor(",varList2,")",sep=""))
                mod0 = model.matrix(modelThis, data=dgeF$samples)
                #mod0 = model.matrix(~as.factor(experimentNo), data=dgeF$samples)
                #svObj=sva(lcpmT[i,samId],mod,mod0)
                #cat(colnames(mod),"\n",sep=", ")
                #cat(colnames(mod0),"\n",sep=", ")
                svObj=try(sva(lcpmT[i,samId],mod,mod0))
                if (class(svObj)=="try-error") next
                #svObj = svaseq(10^6*cpmT[i,samId],mod,mod0)
                nm=NULL
                if (is.matrix(svObj$sv)) {
                    nm=c(colnames(mod),paste("sv",1:ncol(svObj$sv),sep=""))
                } else if (is.numeric(svObj$sv)) {
                    if (length(svObj$sv)>1) {
                        nm=c(colnames(mod),"sv1")
                    } else {
                        nm=NULL
                    }
                }
                if (!is.null(nm)) {
                    #dat=svObj$sv; colnames(dat)=paste("sv",1:ncol(dat),sep="")
                    design=cbind(mod,svObj$sv)
                    colnames(design)=nm
                    #design=model.matrix(~group,data=dgeF$samples)
                    dat <- voom(dgeF,design,save.plot=F)
                    #save(dat,design,file=paste("voom",fName3,".RData",sep=""))
                    fit <- eBayes(lmFit(dat,design))
                    fit0=fit
                    rm(dat)
                    colId=2
                    top=cbind(geneId=rownames(fit$coef),logFC=fit$coef[,colId],PValue=fit$p.value[,colId])
                    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                } else {
                    cat("Could not run SVA !!!\n",sep="")
                }
            }
        }

        if (!is.null(varList2)) {
            minCnt=1
            if (compFlag=="_gammaVsSham" & subsetFlag=="_A1014A1402") {
                j=which(dgeT$samples$treatment2%in%c("gamma-rays","0Sham") & dgeT$samples$experimentNo%in%c("A1014"))
                dgeF=dgeT[i,j]
                #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
                #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
                dgeF$samples$group=as.factor(as.character(dgeF$samples$treatment2))
                fName3=paste("_gammaVsSham_A1014",fName2,"_voom_minCnt",minCnt,sep="")
                #design = model.matrix(~as.factor(treatment2), data=dgeF$samples)
                mod = model.matrix(~as.factor(treatment2), data=dgeF$samples)
                mod0 = model.matrix(~1, data=dgeF$samples)
                svObj=sva(lcpmT[i,j],mod,mod0)
                #svObj = svaseq(10^6*cpmT[i,j],mod,mod0)
                dat=svObj$sv; colnames(dat)=paste("sv",1:ncol(dat),sep="")
                design = cbind(mod,dat)
                dat <- voom(dgeF,design,save.plot=F)
                #save(dat,design,file=paste("voom",fName3,".RData",sep=""))
                fit <- eBayes(lmFit(dat,design))
                fit1=fit
                rm(dat)
                colId=2
                top=cbind(geneId=rownames(fit$coef),logFC=fit$coef[,colId],PValue=fit$p.value[,colId])
                write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)


                i=which(maxCntVec>=minCnt)
                j=which(dgeT$samples$treatment2%in%c("gamma-rays","0Sham") & dgeT$samples$experimentNo%in%c("A1402"))
                dgeF=dgeT[i,j]
                #dgeF$samples$treatment2[which(dgeF$samples$treatment2=="Sham")]="0Sham"
                #dgeF$samples$treatExpt=paste(dgeF$samples$treatment2,dgeF$samples$experimentNo)
                dgeF$samples$group=as.factor(as.character(dgeF$samples$treatment2))
                fName3=paste("_gammaVsSham_A1402",fName2,"_voom_minCnt",minCnt,sep="")
                #design = model.matrix(~as.factor(treatment2), data=dgeF$samples)
                mod = model.matrix(~as.factor(treatment2), data=dgeF$samples)
                mod0 = model.matrix(~1, data=dgeF$samples)
                svObj=sva(lcpmT[i,j],mod,mod0)
                #svObj = svaseq(10^6*cpmT[i,j],mod,mod0)
                dat=svObj$sv; colnames(dat)=paste("sv",1:ncol(dat),sep="")
                design = cbind(mod,dat)
                dat <- voom(dgeF,design,save.plot=F)
                #save(dat,design,file=paste("voom",fName3,".RData",sep=""))
                fit <- eBayes(lmFit(dat,design))
                fit2=fit
                rm(dat)
                colId=2
                top=cbind(geneId=rownames(fit$coef),logFC=fit$coef[,colId],PValue=fit$p.value[,colId])
                write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

                fName3=paste("_",ifelse(is.null(grpUniq),varList,paste(grpName[2],"Vs",capWords(grpName[1]),sep="")),subsetFlag,fName2,sep="")

                colId=2
                png(paste("coef",fName3,".png",sep=""))
                par(mfrow=c(2,2))
                lim=range(c(fit0$coef[,colId],fit1$coef[,colId],fit2$coef[,colId]),na.rm=T)
                x1=fit1$coef[,colId]; x2=fit2$coef[,colId]
                plot(x1,x2,xlim=lim,ylim=lim,main=paste("Gamma vs. sham\nPearson corr ",round(cor(x1,x2,use="complete.obs",method="pearson"),2),sep=""),xlab="A1014: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                x1=fit0$coef[,colId]; x2=fit1$coef[,colId]
                plot(x1,x2,xlim=lim,ylim=lim,main=paste("Gamma vs. sham\nPearson corr ",round(cor(x1,x2,use="complete.obs",method="pearson"),2),sep=""),xlab="A1402A1014 SV: Coef",ylab="A1014: Coef")
                abline(c(0,1))
                x1=fit0$coef[,colId]; x2=fit2$coef[,colId]
                plot(x1,x2,main=paste("Gamma vs. sham\nPearson corr ",round(cor(x1,x2,use="complete.obs",method="pearson"),2),sep=""),xlab="A1014A1402 SV: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                dev.off()


                colId=2
                png(paste("pv",fName3,".png",sep=""))
                par(mfrow=c(2,2))
                lim=range(c(fit0$p.value[,colId],fit1$p.value[,colId],fit2$p.value[,colId]),na.rm=T)
                plot(fit1$p.value[,colId],fit2$p.value[,colId],xlim=lim,ylim=lim,main="Gamma vs. sham",xlab="A1014: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                plot(fit0$p.value[,colId],fit2$p.value[,colId],xlim=lim,ylim=lim,main="Gamma vs. sham",xlab="A1402A1014 SV: Coef",ylab="A1014: Coef")
                abline(c(0,1))
                plot(fit0$p.value[,colId],fit1$p.value[,colId],main="Gamma vs. sham",xlab="A1014A1402 SV: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                dev.off()

                colId=2
                png(paste("histogram_pv",fName3,".png",sep=""))
                par(mfrow=c(2,2))
                lim=range(c(fit0$p.value[,colId],fit1$p.value[,colId],fit2$p.value[,colId]),na.rm=T)
                hist(fit0$p.value[,colId],main="A1402A1014 SV: Gamma vs. sham",xlab="P-value")
                hist(fit1$p.value[,colId],main="A1014: Gamma vs. sham",xlab="P-value")
                hist(fit2$p.value[,colId],main="A1402:Gamma vs. sham",xlab="P-value")
                dev.off()


                colId=2
                png(paste("pv",fName3,".png",sep=""))
                par(mfrow=c(2,2))
                lim=c(-3,3)
                plot(alcpmT[i1],log10(stat_1$pv)-log10(statW$pv[i1]),ylim=lim,main=paste(subsetName,compName,"\nN = ",nrow(stat_1),sep=""),xlab="Avg logCPM",ylab="log10(PV-limma)-log10(PV-wilcox)",cex=.2,pch=20,cex.main=2.5,cex.axis=2,cex.lab=1.5)
                abline(h=0,col="red")
                lim=range(c(fit0$p.value[,colId],fit1$p.value[,colId],fit2$p.value[,colId]),na.rm=T)
                plot(fit1$p.value[,colId],fit2$p.value[,colId],xlim=lim,ylim=lim,main="Gamma vs. sham",xlab="A1014: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                plot(fit0$p.value[,colId],fit2$p.value[,colId],xlim=lim,ylim=lim,main="Gamma vs. sham",xlab="A1402A1014 SV: Coef",ylab="A1014: Coef")
                abline(c(0,1))
                plot(fit0$p.value[,colId],fit1$p.value[,colId],main="Gamma vs. sham",xlab="A1014A1402 SV: Coef",ylab="A1402: Coef")
                abline(c(0,1))
                dev.off()


                pThres=0.05
                table(fit1$p.value[,colId]<pThres,fit2$p.value[,colId]<pThres)
                table(fit0$p.value[,colId]<pThres,fit1$p.value[,colId]<pThres)
                table(fit0$p.value[,colId]<pThres,fit2$p.value[,colId]<pThres)
            }
        }
    }
}


## -----------------------------------------
## NOT USED

datadir=""
tbl1=read.table(paste(datadir,"stat_slope_hzeTe_rnaSeq_p53_voom_minCnt1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl2=read.table(paste(datadir,"stat_slope_hzeTe_rnaSeq_p53_voom_sva_minCnt1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl3=read.table(paste(datadir,"stat_slope_cape_rnaSeq_p53_voom_minCnt1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

datadir="results/rnaSeq/slope/"
x1=dir(datadir)
datadir="results/rnaSeq/deGene/slope/"
x2=dir(datadir)
datadir="results/rnaSeq/treatAndExptNoAdj/slope/"
x3=dir(datadir)

