library(CNEr)
library(rtracklayer)

# Prepare the files under data/
axtFn <- file.path(system.file("extdata", package="CNEr"), 
                 "hg19.danRer7.net.axt")
axtHg19DanRer7 <- readAxt(axtFn)
save(axtHg19DanRer7,
     file="/Users/gtan/Repos/github/CNEr/data/axtHg19DanRer7.rda")

axtFn <- file.path(system.file("extdata", package="CNEr"), 
                   "danRer7.hg19.net.axt")
axtDanRer7Hg19 <- readAxt(axtFn)
save(axtDanRer7Hg19, 
     file="/Users/gtan/Repos/github/CNEr/data/axtDanRer7Hg19.rda")

# Prepare the Axt alignments for the regions of barhl2 and sox14
## hg38 vs danRer10
axtFnHg38DanRer10 <- "/Users/gtan/Repos/github/CNEr/inst/extdata/hg38.danRer10.net.axt.gz"
axtFnDanRer10Hg38 <- "/Users/gtan/Repos/github/CNEr/inst/extdata/danRer10.hg38.net.axt.gz"
axtHg38DanRer10 <- readAxt(axtFnHg38DanRer10)
axtDanRer10Hg38 <- readAxt(axtFnDanRer10Hg38)
qSize <- fetchChromSizes("hg38")
qSize <- seqlengths(qSize["chr6"])
## subAxt
danRer10.hg38.net.axt <- subAxt(axtDanRer10Hg38, chr="chr6", start=24000000, end=27000000L, select="target")
writeAxt(danRer10.hg38.net.axt, "~/danRer10.hg38.net.axt")
hg38.danRer10.net.axt <- subAxt(axtHg38DanRer10, chr="chr6", start=24000000, end=27000000L, select="query", qSize)
writeAxt(hg38.danRer10.net.axt, "~/hg38.danRer10.net.axt")
