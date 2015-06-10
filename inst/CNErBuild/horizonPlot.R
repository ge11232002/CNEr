> xyplot(y~x, panel = function(x,y,...){for (i in 0:round(max(y)/2,0)) panel.xyarea(x,y=ifelse(y>0,y,NA)-(scale * i),col="green",border="green", alpha=0.4)})

fooColors=c("red","blue","yellow")
scale=2
xyplot(y~x, xlab=NULL,ylab=NULL, ylim=c(0, scale),origin=0,
panel = function(x,y,...){for (i in 0:round(max(y)/scale,0)) panel.xyarea(x,y=ifelse(y>0,y,NA)-(scale * i),col=fooColors[i+1],border="green")})
xyplot(y~x, xlab=NULL,ylab=NULL, ylim=c(0, scale),origin=0,
panel = function(x,y,...){for (i in 1:(round(max(y)/scale,0)+1)) panel.xyarea(x,y=ifelse(y>0,y,NA)-(scale * (i-1)),col=fooColors[i],border="green")})


rm Gviz_1.4.4.tar.gz
R CMD build  Gviz/
R CMD INSTALL Gviz_1.4.4.tar.gz
remove.packages("Gviz")
library(Gviz)
data(twoGroups)
dTrack <- DataTrack(twoGroups, name = "uniform")
plotTracks(dTrack)

library(lattice)
library(GenomicRanges)
library(latticeExtra)
set.seed(10)
foo = GRanges(seqnames=seqnames(twoGroups), ranges=ranges(twoGroups), y=runif(25,1,6))
dTrack <- DataTrack(foo, name = "uniform", type="horizon", fill.horizonScale=2, ylim=c(0,2))
plotTracks(dTrack)


ans = queryAnnotationSQLite(dbname="/Users/gtan/Dropbox/Project/CSC/CNEr/geneAnnotation.sqlite", tablename="hg19_refGene", chr="chr3", start=158500001, end=160200000)
grtrack <- GeneRegionTrack(ans, genome = "hg19", chromosome = "chr3", name = "foo")
axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr3")
plotTracks(list(ideoTrack,axisTrack,grtrack), collapseTranscripts = TRUE, shape = "arrow", extend.left = 20000, showId = TRUE)






