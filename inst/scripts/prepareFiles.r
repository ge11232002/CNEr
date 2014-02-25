library(CNEr)
library(rtracklayer)
axtFilesHg19DanRer7 = list.files(path="/Users/gtan/CSC/CNEr/axtNet",
                                                                  pattern=".*hg19\\.danRer7\\.*", full.names=TRUE)
axtHg19DanRer7 = readAxt(axtFilesHg19DanRer7)
axtFilesDanRer7Hg19 = list.files(path="/Users/gtan/CSC/CNEr/axtNet",
                                                                  pattern=".*danRer7\\.hg19\\.*", full.names=TRUE)
axtDanRer7Hg19 = readAxt(axtFilesDanRer7Hg19)

qSize = fetchChromSizes("hg19")
qSize = seqlengths(qSize["chr11"])
type="any"
select="target"

axt1 = subAxt(axtHg19DanRer7, chr="chr11", start=31000000L, end=32500000L, select="target", type="any")
hg19.danRer7.net.axt = axt1
writeAxt(hg19.danRer7.net.axt, "~/hg19.danRer7.net.axt")
axt2 = subAxt(axtDanRer7Hg19, chr="chr11", start=31000000L, end=32500000L, select="query", type="any", qSize)
danRer7.hg19.net.axt = axt2
writeAxt(danRer7.hg19.net.axt, "~/danRer7.hg19.net.axt")


