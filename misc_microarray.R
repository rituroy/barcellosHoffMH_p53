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

#################################################
## Create phen data
#################################################

dataset="A1020"

celFileName <- dir(path=paste(getwd(),"data/cel/",dataset,"/",sep="/"))
id <- sub(".CEL","",celFileName,fixed=T)
id <- gsub(" +","_",sub(".CEL","",celFileName,fixed=T))

datadir <- paste("docs/",dataset,"/",sep="")
phen1 <- read.table(file = paste(datadir,"Annotation-1.txt", sep=""), header = T, sep = "\t", quote="", comment.char="",as.is=T)
phen2 <- read.table(file = paste(datadir,"Annotation.txt", sep=""), header = T, sep = "\t", quote="", comment.char="",as.is=T)

names(phen1)[match(c("X","File.name","Treatment","Group"),names(phen1))]=c("order","fileName","treatment","group")
names(phen2)[match(c("order","Sample.ID","Sample.Name..unique.sample.name.","Experiment..","Treatment"),names(phen2))]=c("order","sampleId","sampleName","experimentNo","treatment")

phen1$treatment[which(phen1$treatment=="Sham")]="sham"
phen2$treatment[which(phen2$treatment=="Sham")]="sham"

table(phen2$experimentNo)
j=phen2$sampleId%in%phen2$sampleId[duplicated(phen2$sampleId)]
phen2[j,][order(phen2$sampleId[j],phen2$experimentNo[j]),]
table(phen1$treatment,phen1$group)
table(phen1$treatment,phen1$group, phen2$experimentNo)

phen <- data.frame(id=phen2$sampleName,phen1,phen2[,which(!names(phen2)%in%names(phen1))],celFileName)
write.table(phen,file=paste("phen_",dataset,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

#################################################
## Output normalized data
#################################################

library(affy)
load("eset.Rdata")

tbl=data.frame(affyId=geneNames(eset),exprs(eset))
datadir <- paste(getwd(),"/",sep="")
write.table(tbl,file=paste(datadir, "normalizedData.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)

#################################################
