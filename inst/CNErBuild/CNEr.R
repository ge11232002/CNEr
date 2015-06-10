############################ build package  ######################################

R_dev CMD build CNEr
R_dev CMD build --no-build-vignettes --no-manual CNEr
R_dev CMD INSTALL CNEr_1.5.3.tar.gz
R_dev CMD check CNEr_1.5.3.tar.gz

library(devtools)
reload(inst("CNEr"))

### readAxt
library(CNEr)
axtFileshg19danRer7 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", pattern=".*hg19\\.danRer7\\.*", full.names=TRUE)
axtshg19danRer7 = readAxt(axtFileshg19danRer7)
axtFilesdanRer7hg19 = list.files(path="/Users/gtan/CSC/CNEr/axtNet", pattern=".*danRer7\\.hg19\\.*", full.names=TRUE)
axtsdanRer7hg19 = readAxt(axtFilesdanRer7hg19)

## subAxt
  dyn.load("/Users/gtan/Repositories/Bitbucket/CNEr/src/alignment.so")
  .Call("subAlignment", c(92822L, 95115L), c(92873L, 95180L),
  c("CAAAACCAGATGCTGTGAGAATACTTTATTAGT----CAAAACCGCATA---CTATAAA", "TGCCAACCTTGGCGCCGATCTGATTCCCGCACTGCCCGATCTGCGTGAGCACGATCTCCCTCATGG"),
  c(1812440L, 47866092L), c(1812495L, 47866157L),
  c("CAAAAC---ATATCATAACTGTACCTTGTTTGTTCCACAAGATTGCATCTTTCCTTAAA",
  "TTCCTACCTTGGCACCAATCTGGTTGCCGCACTGTCCAGCCTGTAAATGCACGATCTCCCTCATTG"),
  c(92863L, 95115L), c(92873L, 95120L),
  c(59L, 66L))


subAxt(x, "chr10", 92866L, 92873L, select="target")
mm10ChromSizes <- fetchChromSizes("mm10")
subAxt(x, chr="chr2", start=1812441, end=1812494, select="query", qSize=seqlengths(mm10ChromSizes["chr2"]))
subAxt(x, chr="chr14", start=77037081, end=77039623, select="query", qSize=seqlengths(mm10ChromSizes["chr14"]))

subAxt(axtHg19DanRer7, chr="chr11", start=31500113, end=31500120, select="target")

foo2 = .subAxtMultiple(axtHg19DanRer7, chr="chr11", start=31500113, end=31500120, select="target")

foo3 = .subAxtMultiple(axtHg19DanRer7, chr="chr11", start=c(31082021, 32461267), end=c(31082862,32461581), select="target")


## readBedToGRanges
bedhg19 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.hg19.bed")
beddanRer7 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.danRer7.bed")
library(rtracklayer)
qSizesdanRer7 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/danRer7.2bit"))
qSizeshg19 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/hg19.2bit"))

## ceScan
# debug
axts = axtshg19danRer7
tFilter= bedhg19
qFilter= beddanRer7
qSizes= qSizesdanRer7
winSize=50
minScore=45
resFiles = tempfile(pattern = paste(minScore, winSize, "ceScan", sep="-"), tmpdir = tempdir(), fileext = "")
foo = .Call("myCeScan",  as.character(seqnames(tFilter)), 
      start(tFilter), end(tFilter), 
      as.character(seqnames(qFilter)), start(qFilter), end(qFilter),
              as.character(seqnames(qSizes)), as.integer(seqlengths(qSizes)),
              as.character(seqnames(targetRanges(axts))),
              start(targetRanges(axts)), end(targetRanges(axts)),
              as.character(strand(targetRanges(axts))),
              as.character(targetSeqs(axts)),
              as.character(seqnames(queryRanges(axts))),
              start(queryRanges(axts)), end(queryRanges(axts)),
              as.character(strand(queryRanges(axts))),
              as.character(querySeqs(axts)),
              score(axts), symCount(axts), winSize, minScore,
              as.character(resFiles))
foo = .Call("myCeScan",  NULL, 
      NULL, NULL, 
      as.character(seqnames(qFilter)), start(qFilter), end(qFilter),
              as.character(seqnames(qSizes)), as.integer(seqlengths(qSizes)),
              as.character(seqnames(targetRanges(axts))),
              start(targetRanges(axts)), end(targetRanges(axts)),
              as.character(strand(targetRanges(axts))),
              as.character(targetSeqs(axts)),
              as.character(seqnames(queryRanges(axts))),
              start(queryRanges(axts)), end(queryRanges(axts)),
              as.character(strand(queryRanges(axts))),
              as.character(querySeqs(axts)),
              score(axts), symCount(axts), winSize, minScore,
              as.character(resFiles))
              
# end of debug
CNEhg19_danRer7 = ceScan(axtshg19danRer7, qFilter=beddanRer7, qSizes=qSizesdanRer7, thresholds=c("45,50", "48,50", "49,50"))
CNEhg19_danRer7_2 = ceScan(axtshg19danRer7, bedhg19, beddanRer7, qSizesdanRer7, thresholds=c("45,50", "48,50", "49,50"))

CNEdanRer7_hg19 = ceScan(axtsdanRer7hg19, beddanRer7, bedhg19, qSizeshg19, thresholds=c("45,50", "48,50", "49,50"))

## ceScan File
bedhg19File = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
beddanRer7File = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
CNEhg19_danRer7 = ceScanFile(axtFileshg19danRer7, bedhg19File, beddanRer7File , qSizesdanRer7, thresholds=c("45,50", "48,50", "49,50"))
CNEdanRer7_hg19 = ceScanFile(axtFilesdanRer7hg19, beddanRer7File, bedhg19File, qSizeshg19, thresholds=c("45,50", "48,50", "49,50"))


## ceMerge
data(CNEDanRer7Hg19)
data(CNEHg19DanRer7)
cneMerge(CNEDanRer7Hg19[[1]], CNEHg19DanRer7[[1]])

## blatCNE
assemblyhg19Twobit = "/export/data/goldenpath/hg19/assembly.2bit"
assemblydanRer7Twobit = "/export/data/goldenpath/danRer7/assembly.2bit"
cneBlateddanRer7_hg19 = list()
for(i in 1:length(cneMergeddanRer7hg19)){
  cneBlateddanRer7_hg19[[names(cneMergeddanRer7_hg19)[i]]] = blatCNE(cneMergeddanRer7_hg19[[i]], sub("\\d+_", "", names(cneMergeddanRer7_hg19)[i]), 8, 4, assemblydanRer7Twobit, assemblyhg19Twobit)
}

## saveCNE
tableName = "cne_twoWay_danRer7_hg19_len50_id900_v1"
dbName = "/mnt/biggley/home/gtan/work/debug/CNEr-29-07-2013/cne.sqlite"
for(i in 1:length(cneBlateddanRer7_hg19)){
  tableName = paste0("cne_twoWay_danRer7_hg19_len50_id", as.integer(sub("_\\d+$", "", names(cneBlateddanRer7_hg19)[i]))*20, "_v1")
  saveCNEToSQLite(cneBlateddanRer7_hg19[[i]], dbName, tableName, overwrite=TRUE)
}

# readCNERangesFromSQLite
chr = "chr6"
CNEstart = 19900000
CNEend =  28000000
minLength = 50
dbName = "/Users/gtan/CSC/CNEr/CNESQL/cne.sqlite"
tableName = "cne2wBf_danRer7_hg19_27_30"
fetchedCNERanges = readCNERangesFromSQLite(dbName, tableName, chr, CNEstart, CNEend, whichAssembly="1", minLength)

##CNEAnnotate
windowSize= 300
cne2wBf_danRer7_hg19_21_30 = CNEAnnotate(dbName, "cne2wBf_danRer7_hg19_21_30", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_hg19_40_50 = CNEAnnotate(dbName, "cne2wBf_danRer7_hg19_40_50", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_hg19_49_50 = CNEAnnotate(dbName, "cne2wBf_danRer7_hg19_49_50", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)

cne2wBf_danRer7_tetNig2_21_30 = CNEAnnotate(dbName, "cne2wBf_danRer7_tetNig2_21_30", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_tetNig2_40_50 = CNEAnnotate(dbName, "cne2wBf_danRer7_tetNig2_40_50", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_tetNig2_49_50 = CNEAnnotate(dbName, "cne2wBf_danRer7_tetNig2_49_50", whichAssembly="1", chr, CNEstart, CNEend, windowSize, minLength)
### new Gviz way
genome = "danRer7"
strand = "+"
dataMatrix= cne2wBf_danRer7_hg19_21_30
dTrack1 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=1, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="hg19 21/30")
dataMatrix = cne2wBf_danRer7_hg19_40_50
dTrack2 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=1, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="hg19 45/50")
dataMatrix = cne2wBf_danRer7_hg19_49_50
dTrack3 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=1, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="hg19 49/50")

dataMatrix= cne2wBf_danRer7_tetNig2_21_30
dTrack4 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=2, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="tetNig2 21/30")
dataMatrix = cne2wBf_danRer7_tetNig2_40_50
dTrack5 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=2, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="tetNig2 45/50")
dataMatrix = cne2wBf_danRer7_tetNig2_49_50
dTrack6 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horiz", horizon.scale=2, fill.horizon=c("#B41414", "#E03231", "#F7A99C", "yellow", "orange", "red"), name="tetNig2 49/50")

axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "danRer7", chromosome = chr)
refGeneAnnotation = queryAnnotationSQLite(dbname="/Users/gtan/CSC/CNEr/annotationSQL/geneAnnotation.sqlite", tablename="danRer7_refGene", chr=chr, start= CNEstart, end= CNEend)
refgrtrack <- GeneRegionTrack(refGeneAnnotation, genome = "danRer7", chromosome = chr, name = "refGene")
ensGeneAnnotation = queryAnnotationSQLite(dbname="/Users/gtan/CSC/CNEr/annotationSQL/geneAnnotation.sqlite", tablename="danRer7_ensGene", chr=chr, start= CNEstart, end= CNEend)
ensgrtrack = GeneRegionTrack(ensGeneAnnotation, genome = "danRer7", chromosome = chr, name = "ensGene")

cpgIslands = UcscTrack(genome = "danRer7", chromosome = chr, track = "cpgIslandExt", from=CNEstart, to=CNEend, trackType = "AnnotationTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box", fill = "#006400", name = "CpG Islands")
plotTracks(list(axisTrack,ideoTrack,  refgrtrack, cpgIslands, dTrack1, dTrack2, dTrack3, dTrack4, dTrack5, dTrack6), collapseTranscripts = TRUE, shape = "arrow", from= CNEstart, to= CNEend, showId = TRUE, extend.left = 20000)

### old Gviz way


### prepare gene annotation sqlite3
makeGeneDbFromUCSC(genome="danRer7", tablename="refGene", dbnameSQLite="/Users/gtan/CSC/CNEr/annotationSQL/geneAnnotation.sqlite")
makeGeneDbFromUCSC(genome="danRer7", tablename="ensGene", dbnameSQLite="/Users/gtan/CSC/CNEr/annotationSQL/geneAnnotation.sqlite")


#############################################   readBinary##############
axtFiles ="/mnt/biggley/data/pairwiseAlignments/ucsc/axtNet/hg19.danRer7.net.axt" 
fn = file(axtFiles, "rb")
foo = readBin(fn, raw(), file.info(axtFiles)$size)
rawToChar(foo[1:16]) 
index = grepRaw("\n", foo, fixed=TRUE)

targetRanges="GRanges", targetSeqs="DNAStringSet",queryRanges="GRanges", querySeqs="DNAStringSet", score="integer", symCount="integer"
                        )
                        
                        

