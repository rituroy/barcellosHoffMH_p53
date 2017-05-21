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


## ----------------------
datadir="docs/coral/"

fName="150903_D00372_0318_AC7FNUANXX_Combo_HSQ_98.txt"; offset=4; nrow=189
tmp=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1,skip=offset)
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=nrow-offset-1,skip=offset+1)
paste(unlist(tmp),collapse=",")
names(tbl)=c("lane","sampleId","sampleRef","index","description","control","project","yield","percPF","numReads","percRawClustPerLane","percPerfectIndexReads","percOneMismatchReads","percQ30Bases","meanQualityScore")
phen1=tbl

fName="150903_D00372_0318_AC7FNUANXX_Combo_HSQ_98.txt"; offset=192; nrow=257
tmp=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1,skip=offset)
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=nrow-offset-1,skip=offset+1)
paste(unlist(tmp),collapse=",")
names(tbl)=c("id","recipe","operator","directory")
phen2=tbl

fName="RNA sequencing_Coral Omene_Barcellos-Hoff.txt"
tmp=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=1)
paste(unlist(tmp),collapse=",")
names(tbl)=c("tubeLabel","sampleId","username","dateTime","nucleicAcid","unit","A260Abs","A280Abs","260_280","260_230","volume")
phen3=tbl

phen=phen1
phen=phen2
for (k in 1:ncol(phen)) {
    print(range(nchar(phen[,k])))
}
phen=phen2
k2=c()
for (k in 1:ncol(phen)) {
    if (any(!is.na(phen[,k]))) k2=c(k2,k)
}
phen=phen[,k2]
phen2=phen

fileList=dir(datadir,pattern="fastq.gz")
fileListC=fileList

datadir=""
fileList3=read.table(paste(datadir,"tmp3.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
fileList3=fileList3[,1]
fileList3=fileList3[grep("fastq",fileList1)]
table(fileList3%in%fileList)
table(fileList%in%fileList3)

phen=phen1
id1=paste(phen$sampleId,"_",phen$index,"_L",formatC(phen$lane,width=3,flag="0"),sep="")
id2=sub("_R1_001.fastq.gz","",fileList,fixed=T)
table(id2%in%id1)
x1=unique(phen$sampleId[id1%in%id2])
length(x1) ## No. of samples
x2=unique(phen$sampleId[!id1%in%id2])
x1[x1%in%x2] ## Should be NONE
x2
phen=cbind(cbind(id=id1,phen1,stringsAsFactors=F)[match(id2,id1),],fileList,stringsAsFactors=F)
phenC=phen

table(phen2$id%in%phen1$sampleId)


## ----------------------
datadir="docs/haoxu/"

fName="Haoxu-1020Aged-Data-RNASeq-08115150803_D00372_0306_AC7N2JANXX_Combo_HSQ_86.txt"; offset=3; nrow=150
tmp=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1,skip=offset)
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=nrow-offset-1,skip=offset+1)
paste(unlist(tmp),collapse=",")
names(tbl)=c("lane","sampleId","sampleRef","index","description","control","project","yield","percPF","numReads","percRawClustPerLane","percPerfectIndexReads","percOneMismatchReads","percQ30Bases","meanQualityScore")
phen1=tbl

fName="Haoxu-1020Aged-Data-RNASeq-08115150803_D00372_0306_AC7N2JANXX_Combo_HSQ_86.txt"; offset=154; nrow=228
tmp=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1,skip=offset)
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=nrow-offset-1,skip=offset+1)
paste(unlist(tmp),collapse=",") ## ignore
names(tbl)=c("id","recipe","operator","directory")
phen2=tbl

phen=phen1
phen=phen2
for (k in 1:ncol(phen)) {
    print(range(nchar(phen[,k])))
}

phen=phen1
for (k in 1:ncol(phen)) {
    if (is.character(phen[,k])) {
        phen[,k]=sub(" +$","",phen[,k])
    }
}
phen1=phen

phen=phen2
k2=c()
for (k in 1:ncol(phen)) {
    if (any(!is.na(phen[,k]))) {
        k2=c(k2,k)
        if (is.character(phen[,k])) {
            phen[,k]=sub(" +$","",phen[,k])
        }
    }
}
phen=phen[,k2]
phen2=phen

fileList=dir(datadir,pattern="fastq.gz")
fileListH=fileList

datadir="tmp2/"
fileList1=read.table(paste(datadir,"tmp1.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
fileList2=read.table(paste(datadir,"tmp2.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
fileList3=read.table(paste(datadir,"tmp4.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
fileList1=fileList1[,1]
fileList2=fileList2[,1]
fileList1=fileList1[grep("fastq",fileList1)]
fileList2=fileList2[grep("fastq",fileList2)]
fileList1=paste(fileList1,".gz",sep="")
table(fileList1%in%fileList2)
table(fileList2%in%fileList1)
fileList1[!fileList1%in%fileList2]

table(fileList1%in%unique(c(fileList2,fileList3)))
table(unique(c(fileList2,fileList3))%in%fileList1)

fileList=fileList2
fileList=fileList1

phen=phen1
id1=paste(phen$sampleId,"_",phen$index,"_L",formatC(phen$lane,width=3,flag="0"),sep="")
id2=sub("_R1_001.fastq.gz","",fileList,fixed=T)
table(id2%in%id1)
x1=unique(phen$sampleId[id1%in%id2])
length(x1) ## No. of samples
x2=unique(phen$sampleId[!id1%in%id2])
x1[x1%in%x2] ## Should be NONE
x2
phen=cbind(cbind(id=id1,phen1,stringsAsFactors=F)[match(id2,id1),],fileList,stringsAsFactors=F)
phenH=phen

table(phen2$id%in%phen1$sampleId)

id3=sapply(fileList,function(x) {paste(strsplit(x,"_")[[1]][1:2],collapse="_")},USE.NAMES=F)
table(duplicated(id3))

#write.table(cbind(geneId=rownames(cnt_1),cnt_1), file="count_norm_Homo_sapiens.txt",col.names=T,row.names=F, sep="\t",quote=F)


paste(sub(".gz","",fileList1[!fileList1%in%fileList2],fixed=T),collapse=" ")

