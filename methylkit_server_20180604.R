# source("https://bioconductor.org/biocLite.R")
# biocLite(c("GenomicRanges","IRanges"))
# biocLite("methylKit",dependencies=TRUE)
# 
# setwd("/media/donglei/New Volume/Users/Donglei/Documents")

library(methylKit)


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


# Filter samples based on coverage

# discards bases that have coverage:

# 1. below 10X (a high enough read coverage will increase the power of the statistical tests)
# 2. more than 99.9th percentile of coverage in each sample (samples might be suffering from PCR bias)

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

# merge all samples to one object for base-pair locations that are covered in all sample
meth=unite(myobj, destrand=TRUE)


# summarize methylation information over tiling windows rather than doing base-pair resolution analysis.
# tiles the genome with windows 1000bp length and 1000bp step-size 

tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)


# Finding differentially methylated bases or regions (adjust for overdispersion)

# by base

system.time(myDiff<-calculateDiffMeth(meth, overdispersion="MN",test="Chisq",mc.cores=24))
# get hyper methylated bases
system.time(myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper"))
# get hypo methylated bases
system.time(myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo"))
# get all differentially methylated bases
system.time(myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01))

# by region

system.time(myDiff.region<-calculateDiffMeth(tiles, overdispersion="MN",test="Chisq",mc.cores=24))
# get hyper methylated bases
system.time(myDiff25p.hyper.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01,type="hyper"))
# get hypo methylated bases
system.time(myDiff25p.hypo.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01,type="hypo"))
# get all differentially methylated bases
system.time(myDiff25p.region=getMethylDiff(myDiff.region,difference=25,qvalue=0.01))
