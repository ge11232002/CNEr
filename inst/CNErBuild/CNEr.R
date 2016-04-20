## readBedToGRanges
bedhg19 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.hg19.bed")
beddanRer7 = readBed("/Users/gtan/CSC/CNEr/filters/filter_regions.danRer7.bed")
library(rtracklayer)
qSizesdanRer7 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/danRer7.2bit"))
qSizeshg19 = seqinfo(TwoBitFile("/Users/gtan/CSC/CNEr/2bit/hg19.2bit"))

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


                        

