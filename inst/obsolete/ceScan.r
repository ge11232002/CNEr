
## The core program of ceScan



tFilterFile = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
qFilterFile = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
queryFile = "/export/data/goldenpath/hg19/assembly.2bit"

tFilter = readBed(tFilterFile)
qFilter = readBed(qFilterFile)

axtFile = "/export/downloads/ucsc/axtNet/danRer7/danRer7.oryLat2.net.axt.gz"


## we start from targetSeqs and targetRanges, querySeqs and queryRanges
tFilter = readBed(tFilterFile)
qFilter = readBed(qFilterFile)
chromSizes = seqinfo(TwoBitFile(queryFile))
chromSizes = seqlengths(chromSizes)
qFilterRev = makeReversedFilter(qFilter, chromSizes)
targetOverlaps = findOverlaps(targetRanges, tFilter)
qFilter = c(qFilter, qFilterRev)
queryOverlaps = findOverlaps(queryRanges, qFilter)

targetStarts = start(targetRanges)
targetEnds = end(targetRanges)
queryStarts = start(queryRanges)
queryEnds = end(queryRanges)
tFilterStarts = start(tFilter)
tFilterEnds = end(tFilter)
qFilterStarts = start(qFilter)
qFilterEnds = end(qFilter)
# process one alignment
i = 1
system.time(for(i in 1:1000){
tFilterNow = subjectHits(targetOverlaps)[queryHits(targetOverlaps) == i]
if(length(tFilterNow) == 0){
  tFilterStart = -1
  tFilterEnd = -1
}else{
  tFilterStart = tFilterStarts[tFilterNow]
  tFilterEnd = tFilterEnds[tFilterNow]
}
qFilterNow = subjectHits(queryOverlaps)[queryHits(queryOverlaps) == i]
if(length(qFilterNow) == 0){
  qFilterStart = -1
  qFilterEnd = -1
}else{
  qFilterStart = qFilterStarts[qFilterNow]
  qFilterEnd = qFilterEnds[qFilterNow]
}
})
dyn.load("ceScan.so")
.C("ceScan", 
      as.character(targetSeqs[i]), 
      as.character(querySeqs[i]),
      targetStart=targetStarts[i],
      targetEnd=targetEnds[i],
      queryStart=queryStarts[i],
      queryEnd=queryEnds[i],
      tFilterStart=tFilterStart,
      tFilterEnd=tFilterEnd,
      qFilterStart=qFilterStart,
      qFilterEnd=qFilterEnd,
      output=as.character("output")
      )


####################
dyn.load("src/")
