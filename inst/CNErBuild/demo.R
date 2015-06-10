### Human VS Dog
selfDir = "~/Repos/CSC/CNEr/R"
selfScripts = list.files(path=selfDir, pattern='.*\\.r$', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}

axtFiles = list.files(path="/export/downloads/ucsc/axtNet/hg19", pattern=".*hg19\\.canFam3\\.*", full.names=TRUE)
axtFiles = axtFiles[1:5]
axts_hg19_canFam3 = readAxt(axtFiles)
axtFiles = list.files(path="/export/downloads/ucsc/axtNet/canFam3", pattern=".*canFam3\\.hg19\\.*", full.names=TRUE)
axts_canFam3_hg19 = readAxt(axtFiles)

bed_hg19 = readBedToGRanges("/export/data/CNEs/hg19/filters/filter_regions.hg19.bed")
bed_canFam3 = readBedToGRanges("/export/data/CNEs/canFam3/filters/filter_regions.canFam3.bed")
qSizes_canFam3 = seqinfo(TwoBitFile("/export/data/goldenpath/canFam3/assembly.2bit"))
qSizes_hg19 = seqinfo(TwoBitFile("/export/data/goldenpath/hg19/assembly.2bit"))

CNE_hg19_canFam3 = ceScan(axts_hg19_canFam3, bed_hg19, bed_canFam3, qSizes_canFam3, thresholds=c("29,30", "30,30", "35,50", "40,50", "45,50", "48,50", "49,50"))
CNE_canFam3_hg19 = ceScan(axts_canFam3_hg19, bed_canFam3, bed_hg19, qSizes_hg19, thresholds=c("29,30", "30,30", "35,50", "40,50", "45,50", "48,50", "49,50"))

cneMerged_canFam3_hg19 = mapply(ceMerge, CNE_canFam3_hg19, CNE_hg19_canFam3, SIMPLIFY=FALSE)

assembly_hg19_Twobit = "/export/data/goldenpath/hg19/assembly.2bit"
assembly_canFam3_Twobit = "/export/data/goldenpath/canFam3/assembly.2bit"

cneBlated_canFam3_hg19 = list()
for(i in 1:length(cneMerged_canFam3_hg19)){
  cneBlated_canFam3_hg19[[names(cneMerged_canFam3_hg19)[i]]] = blatCNE(cneMerged_canFam3_hg19[[i]], sub("\\d+_", "", names(cneMerged_canFam3_hg19)[i]), 4, 4, assembly_canFam3_Twobit, assembly_hg19_Twobit)
}


# readCNERangesFromSQLite
chr = "chr16"
CNEstart = 45000000
CNEend = 60000000
minLength = 50
fetchedCNERanges = readCNERangesFromSQLite(dbName, tableName, chr, CNEstart, CNEend, whichAssembly="2", minLength)

##CNEAnnotate
windowSize= 300
dbName = "/Users/gtan/Dropbox/Project/CSC/CNEr/cne.sqlite"
hg19_canFam3_len50_id900_300 = CNEAnnotate(dbName, "cne_twoWay_canFam3_hg19_len50_id900_v1", whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)
hg19_canFam3_len50_id960_300 = CNEAnnotate(dbName, "cne_twoWay_canFam3_hg19_len50_id960_v1", whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)
hg19_canFam3_len50_id980_300 = CNEAnnotate(dbName, "cne_twoWay_canFam3_hg19_len50_id980_v1", whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)
listToPlot = list(hg19_canFam3_len50_id900_300 = hg19_canFam3_len50_id900_300, hg19_canFam3_len50_id960_300 = hg19_canFam3_len50_id960_300, hg19_canFam3_len50_id980_300 = hg19_canFam3_len50_id980_300)
p = plotCNE(listToPlot, horizonscale=4, nbands=5)
zoomLevel = c(45000000, 60000000)
p+xlim(range(zoomLevel))

## genomic features
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
aldoa.gr <- GRanges(chr, IRanges(zoomLevel[1], zoomLevel[2]))
library(ggbio)
p1 <- autoplot(txdb, which = aldoa.gr, stat = "reduce")
p2 = tracks(knownGene = p1, CNE = p) + xlim(zoomLevel[1], zoomLevel[2])


### Gviz way
genome = "hg19"
strand = "+"
chr = "chr16"
dataMatrix= hg19_canFam3_len50_id900_300
dTrack1 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horizon", fill.horizonScale=4, ylim=c(0,4), name="hg19_canFam3_len50_id900_300")
dataMatrix = hg19_canFam3_len50_id960_300
dTrack2 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horizon", fill.horizonScale=4, ylim=c(0,4), name="hg19_canFam3_len50_id960_300")
dataMatrix = hg19_canFam3_len50_id980_300
dTrack3 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome, type="horizon", fill.horizonScale=4, ylim=c(0,4), name="hg19_canFam3_len50_id980_300")

axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
geneAnnotation = queryAnnotationSQLite(dbname="/Users/gtan/Dropbox/Project/CSC/CNEr/geneAnnotation.sqlite", tablename="hg19_refGene", chr="chr16", start= 45000000, end= 60000000)
grtrack <- GeneRegionTrack(geneAnnotation, genome = "hg19", chromosome = chr, name = "refGene")
plotTracks(list(ideoTrack, axisTrack, grtrack, dTrack1, dTrack2, dTrack3), collapseTranscripts = TRUE, shape = "arrow", from= 45000000, to=60000000, showId = TRUE, extend.left = 20000)


