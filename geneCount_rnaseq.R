argv=commandArgs(TRUE)
samId=argv[1]

library("R.utils")
use("aroma.seq (>= 0.6.2)")

print(packageVersion("aroma.seq"))
print(packageVersion("R.utils"))
print(packageVersion("R.filesets"))

print(capabilitiesOf(aroma.seq))
#isCapableOf(aroma.seq, "tophat2")

dataset <- "tmp"
organism <- "Mus_musculus"
assembly <- "GRCm38"
release <- 79

dataset <- "Haoxu1020Aged081115"
organism <- "Mus_musculus"
assembly <- "GRCm38"
release <- 79

dataset <- "CoralCapeA1402"
organism <- "Mus_musculus"
assembly <- "GRCm38"
release <- 79

cat("------------ geneCount.R ------------\n")
cat("------------ ",dataset,", ",organism,", ",samId," ------------\n",sep="")

## Setup annotation data
pathA <- file.path("annotationData", "organisms", organism)
filename <- sprintf("%s.%s.dna.chr=1-MT.fa", organism, assembly)
fa <- FastaReferenceFile(filename, path=pathA)
print(fa)

filename <- sprintf("%s.%s.%s.gtf", organism, assembly, release)
gtf <- GtfDataFile(filename, path=pathA)
print(gtf)


## Setup raw data set
pathR <- file.path("fastqData", dataset, organism)
fqs <- IlluminaFastqDataSet$byPath(pathR, paired=FALSE, recursive=TRUE, struct="<rootpath>/<dataset>/<organism>/<sample>/")
print(fqs)

#nbrOfSamples <- length(unique(getNames(fqs)))
#printf("Number of samples: %d\n", nbrOfSamples)

fqs <- fqs[getNames(fqs) == samId]
print(fqs)

#save(fqs,file=paste("bams_",organism,"_",samId,".RData",sep=""))

bams <- doTopHat2(fqs, reference=fa, transcripts=gtf, groupBy="name", tags=c("*"), verbose=-100)
print(bams)

counts <- doHTSeqCount(bams, transcripts = gtf, verbose = -100)
print(counts)
