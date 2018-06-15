source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","IRanges"))
biocLite("methylKit",dependencies=TRUE)

setwd("/media/donglei/New Volume/Users/Donglei/Documents")
library(methylKit)

indata_4OHT_1=processBismarkAln( location = "./HMLE-4OHT-24h-1_bismark.deduplicated.bam.sort.bam",
                                 sample.id="test1", assembly="grch38", 
                                 read.context="CpG", save.folder=getwd())

indata_4OHT_2=processBismarkAln( location = "./HMLE-4OHT-24h-2_bismark.deduplicated.bam.sort.bam",
                                 sample.id="test2", assembly="grch38", 
                                 read.context="CpG", save.folder=getwd())

indata_EtOH_1=processBismarkAln( location = "./HMLE-EtOH-24h-1_bismark.deduplicated.bam.sort.bam",
                                 sample.id="ctrl1", assembly="grch38", 
                                 read.context="CpG", save.folder=getwd())

indata_EtOH_2=processBismarkAln( location = "./HMLE-EtOH-24h-2_bismark.deduplicated.bam.sort.bam",
                                 sample.id="ctrl2", assembly="grch38", 
                                 read.context="CpG", save.folder=getwd())

start.time <- Sys.time()

file.list=list("./test1_CpG.txt","./test2_CpG.txt","./ctrl1_CpG.txt", "./ctrl2_CpG.txt")


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
               sample.id=list("test1","test2","ctrl1","ctrl2"),
               assembly="GRCh38",
               treatment=c(1,1,0,0),
               context="CpG"
)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Descriptive statistics on samples

getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj[[3]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[3]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj[[4]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[4]],plot=TRUE,both.strands=FALSE)

# Filter samples based on coverage

# discards bases that have coverage:

# 1. below 10X (a high enough read coverage will increase the power of the statistical tests)
# 2. more than 99.9th percentile of coverage in each sample (samples might be suffering from PCR bias)

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

# merge all samples to one object for base-pair locations that are covered in all sample
meth=unite(myobj, destrand=TRUE)
head(meth)

# check the correlation between samples
getCorrelation(meth,plot=TRUE)

# cluster the samples based on the similarity of their methylation profiles (dendrogram)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# PCA plots: PC1 vs PC2
PCASamples(meth)


# summarize methylation information over tiling windows rather than doing base-pair resolution analysis.
# tiles the genome with windows 1000bp length and 1000bp step-size 

tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
head(tiles[[1]],3)

# Finding differentially methylated bases or regions (adjust for overdispersion)

# by base

system.time(myDiff<-calculateDiffMeth(meth, overdispersion="MN",test="Chisq",mc.cores=8))
# get hyper methylated bases
system.time(myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper"))
# get hypo methylated bases
system.time(myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo"))
# get all differentially methylated bases
system.time(myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01))

# by region

system.time(myDiff.region<-calculateDiffMeth(tiles, overdispersion="MN",test="Chisq",mc.cores=8))
# get hyper methylated bases
system.time(myDiff25p.hyper.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01,type="hyper"))
# get hypo methylated bases
system.time(myDiff25p.hypo.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01,type="hypo"))
# get all differentially methylated bases
system.time(myDiff25p.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01))

# Annotating differentially methylated bases or regions

library(genomation)

# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                            package = "methylKit"))

# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data

annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)