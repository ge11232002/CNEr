### -----------------------------------------------------------------
### chromosome to colour mapping
### Not exported!
chr2colour <- function(chromosomes){
  ### diverging colour scheme from http://colorbrewer2.org
  divergingColour <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
                       "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
                       "#5e4fa2")
  uniqChromosomes <- as.character(unique(chromosomes))
  chrColours <- rep(divergingColour,
                    ceiling(length(uniqChromosomes) / length(divergingColour)))
  mapping <- chrColours[1:length(uniqChromosomes)]
  names(mapping) <- uniqChromosomes
  return(mapping[as.character(chromosomes)])
}

### -----------------------------------------------------------------
### readAncoraCNEs: this is for internal use in Lenhard group.
### It reads the CNE file (filename format: cne2wBf_hg38_mm10_50_50)
### and return a GRanges.
### Exported!
readAncora <- function(fn, assembly){
  assembly1 <- strsplit(basename(fn), split="_")[[1]][2]
  assembly2 <- strsplit(basename(fn), split="_")[[1]][3]
  cne <- read_tsv(fn, col_names=FALSE)
  
  if(assembly == assembly1){
    ans <- GRanges(seqnames=cne$X1,
                   ranges=IRanges(start=cne$X2+1,
                                  end=cne$X3),
                   strand="*",
                   name=paste0(cne$X4, ":", (cne$X5+1), "-", cne$X6),
                   itemRgb=chr2colour(cne$X4)
                   )
  }else if(assembly == assembly2){
    ans <- GRanges(seqnames=cne$X4,
                   ranges=IRanges(start=cne$X5+1,
                                  end=cne$X6),
                   strand="*",
                   name=paste0(cne$X1, ":", (cne$X2+1), "-", cne$X3),
                   itemRgb=chr2colour(cne$X1)
                   )
  }else{
    stop("Wrongly specified assembly!")
  }
  return(ans)
}
