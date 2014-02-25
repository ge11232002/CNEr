

tFilterFile = "/export/data/CNEs/danRer7/filters/filter_regions.danRer7.bed"
qFilterFile = "/export/data/CNEs/hg19/filters/filter_regions.hg19.bed"
queryFile = "/export/data/goldenpath/hg19/assembly.2bit"

tFilter = readBed(tFilterFile)
qFilter = readBed(qFilterFile)

axtFile = "/export/downloads/ucsc/axtNet/danRer7/danRer7.oryLat2.net.axt.gz"

tFilter = readBed(tFilterFile)
qFilter = readBed(qFilterFile)
chromSizes = seqinfo(TwoBitFile(queryFile))
chromSizes = seqlengths(chromSizes)
qFilterRev = makeReversedFilter(qFilter, chromSizes)
#tFilter = intersect(targetRanges, tFilter)
targetOverlaps = findOverlaps(targetRanges, tFilter)
targetOverlaps = subjectHits(targetOverlaps)

targetOverlaps = split(subjectHits(targetOverlaps), queryHits(targetOverlaps))
qFilter = c(qFilter, qFilterRev)
#qFilter = intersect(queryRanges, qFilter)
queryOverlaps = findOverlaps(queryRanges, qFilter)
queryOverlaps = split(subjectHits(queryOverlaps), queryHits(queryOverlaps))
targetStarts = start(targetRanges)
targetEnds = end(targetRanges)
queryStarts = start(queryRanges)
queryEnds = end(queryRanges)
tFilterStarts = start(tFilter)
tFilterEnds = end(tFilter)
qFilterStarts = start(qFilter)
qFilterEnds = end(qFilter)

## 
comparisonAlignment = compDNAStringSet(targetSeqs, querySeqs)
comparisonAlignment = lapply(comparisonAlignment, as.integer)
liftTarget = seqToAlignment(targetSeqs)
liftQuery = seqToAlignment(querySeqs)

## mask the comparisonAlignment
system.time(
for(i in 1:length(comparisonAlignment[1:1000])){
  tmasksStarts = pmax(tFilterStarts[targetOverlaps[[as.character(i)]]] - targetStarts[i] + 1, 1)
  tmasksEnds = pmin(tFilterEnds[targetOverlaps[[as.character(i)]]], targetEnds[i]) - targetStarts[i] + 1
  tmasksStarts = liftTarget[[i]][tmasksStarts]
  tmasksEnds = liftTarget[[i]][tmasksEnds]
  qmasksStarts = pmax(qFilterStarts[queryOverlaps[[as.character(i)]]] - queryStarts[i] + 1, 1)
  qmasksEnds = pmin(qFilterEnds[queryOverlaps[[as.character(i)]]], queryEnds[i]) - queryStarts[i] + 1
  qmasksStarts = liftQuery[[i]][qmasksStarts]
  qmasksEnds = liftQuery[[i]][qmasksEnds]
  comparisonAlignment[[i]][unlist(mapply(seq, c(tmasksStarts, qmasksStarts), c(tmasksEnds, qmasksEnds)))] = NA
})
       


