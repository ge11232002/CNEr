library(CNEr)
library(rtracklayer)
library(GenomicRanges)
# Prepare the Axt alignments for the regions of barhl2 and sox14
## hg38 vs danRer10
axtFnHg38DanRer10 <- "/Users/gtan/Repos/github/CNEr/inst/extdata/hg38.danRer10.net.axt"
axtFnDanRer10Hg38 <- "/Users/gtan/Repos/github/CNEr/inst/extdata/danRer10.hg38.net.axt"
axtHg38DanRer10 <- readAxt(axtFnHg38DanRer10)
axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38)
qSize <- fetchChromSizes("hg38")
qSize <- seqlengths(qSize["chr6"])
## subAxt
danRer10.hg38.net.axt <- subAxt(axtDanRer10Hg38, chr="chr6", 
                                start=24000000, end=27000000, select="target")
writeAxt(danRer10.hg38.net.axt, "~/danRer10.hg38.net.axt")
hg38.danRer10.net.axt <- subAxt(axtHg38DanRer10, chr="chr6", 
                                start=24000000, end=27000000, select="query", qSize)
writeAxt(hg38.danRer10.net.axt, "~/hg38.danRer10.net.axt")

# Prepare the filter bed files for the regions of barhl2 and sox14
hg38Filter <- readBed("/Users/gtan/Repos/github/CNEr/inst/extdata/filter_regions.hg38.bed")
danRer10Filter <- readBed("/Users/gtan/Repos/github/CNEr/inst/extdata/filter_regions.danRer10.bed")
hits <- findOverlaps(hg38Filter, reduce(c(targetRanges(hg38.danRer10.net.axt),
                                          queryRanges(danRer10.hg38.net.axt))),
                     ignore.strand=TRUE)
hg38Filter <- hg38Filter[unique(queryHits(hits))]
danRer10Filter <- danRer10Filter[seqnames(danRer10Filter) == "chr6" &
                                 start(danRer10Filter) >= 24000000 & 
                                 end(danRer10Filter) <= 27000000]
export.bed(hg38Filter, con="~/filter_regions.hg38.bed")
export.bed(danRer10Filter, con="~/filter_regions.danRer10.bed")

# Prepare the CNE data 
axtFnHg38DanRer10 <- file.path(system.file("extdata", package="CNEr"), 
                               "hg38.danRer10.net.axt")
axtHg38DanRer10 <- readAxt(axtFnHg38DanRer10)
axtFnDanRer10Hg38 <- file.path(system.file("extdata", package="CNEr"), 
                               "danRer10.hg38.net.axt")
axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38)
bedHg38Fn <- file.path(system.file("extdata", package="CNEr"), 
                       "filter_regions.hg38.bed")
bedHg38 <- readBed(bedHg38Fn)
bedDanRer10Fn <- file.path(system.file("extdata", package="CNEr"), 
                           "filter_regions.danRer10.bed")
bedDanRer10 <- readBed(bedDanRer10Fn)
qSizesHg38 <- fetchChromSizes("hg38")
qSizesDanRer10 <- fetchChromSizes("danRer10")
CNEHg38DanRer10 <- ceScan(x=axtHg38DanRer10, tFilter=bedHg38,
                          qFilter=bedDanRer10, tSizes=qSizesHg38,
                          qSizes=qSizesDanRer10,
                          window=c(50,50,50), identity=c(45, 48, 49))
save(CNEHg38DanRer10, file="~/CNEHg38DanRer10.rda")
CNEDanRer10Hg38 <- ceScan(x=axtDanRer10Hg38, tFilter=bedDanRer10,
                          qFilter=bedHg38, tSizes=qSizesDanRer10,
                          qSizes=qSizesHg38,
                          window=c(50,50,50), identity=c(45, 48, 49))
save(CNEDanRer10Hg38, file="~/CNEDanRer10Hg38.rda")

## Prepare danRer10CNE.sqlite
library(CNEr)
cneFns <- file.path("/mnt/biggles/data/CNEs/blatFiltered_19-07-2013/",
                    c("cne2wBf_AstMex102_danRer10_48_50",
                      "cne2wBf_cteIde1_danRer10_75_75",
                      "cne2wBf_danRer10_hg38_21_30",
                      "cne2wBf_danRer10_hg38_45_50",
                      "cne2wBf_danRer10_hg38_49_50"))
dbName <- "danRer10CNE.sqlite"
readAncoraIntoSQLite(cneFns, dbName=dbName, overwrite=FALSE)