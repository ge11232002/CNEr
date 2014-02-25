### -----------------------------------------------------------------
### Fetch the CNE coordinates from SQL and compute the densities
### Exported!
setMethod("CNEDensity", 
          signature(tableName="character", assembly1="character",
                    assembly2="missing", threshold="missing"),
          function(dbName, tableName, assembly1, assembly2, threshold,
                   chr, start, end, windowSize, minLength=NULL){
            if(!grepl("^.+_.+_\\d+_\\d+$", tableName))
              stop("The tableName should be in the format danRer7_hg19_49_50.")
            assemblyNames <- strsplit(tableName, "_")[[1]][1:2]
            stopifnot(assembly1 %in% assemblyNames)
            if(which(assembly1 == assemblyNames) == 1){
              whichAssembly = "L"
            }else{
              whichAssembly = "R"
            }
            .CNEDensityInternal(dbName=dbName, tableName=tableName,
                                whichAssembly=whichAssembly,
                                chr=chr, start=start, end=end,
                                windowSize=windowSize, minLength=minLength)
          }
          )
setMethod("CNEDensity", signature(tableName="missing", assembly1="character",
                                  assembly2="character", threshold="character"),
          function(dbName, tableName, assembly1, assembly2, threshold,
                   chr, start, end, windowSize, minLength=NULL){
            stopifnot(length(threshold) == 1L)
            builtTableName = paste(paste(sort(c(assembly1, assembly2)), 
                                         sep="_", collapse="_"), threshold, 
                                   sep="_", collapse="_")
            if(which(assembly1 == sort(c(assembly1, assembly2))) == 1){
              whichAssembly = "L"
            }else{
              whichAssembly = "R"
            }
            .CNEDensityInternal(dbName=dbName, tableName=builtTableName,
                                whichAssembly=whichAssembly,
                                chr=chr, start=start, end=end,
                                windowSize=windowSize, minLength=minLength)
          }
          )



.CNEDensityInternal <- function(dbName, tableName, whichAssembly=c("L","R"), 
                       chr, start, end, windowSize, 
                       minLength=NULL){
  nrGraphs <- 1
  CNEstart <- start
  CNEend <- end
  # This is the pipeline of doing the density plot
  # The windowSize is in kb.
  whichAssembly <- match.arg(whichAssembly)
  if(!is(windowSize, "integer"))
    stop("windowSize must be an integer!")
  windowSize <- windowSize * 1000
  CNElength <- CNEend - CNEstart + 1
  pixel_width <- 2048
  if(CNElength <= pixel_width) {
    step_size <- 1
  }else{
    step_size <- as.integer(CNElength/pixel_width)
    if(step_size > windowSize/10)
      step_size <- windowSize/10
    while(windowSize %% step_size){
      step_size <- step_size - 1
    }
  }
  # make things easier
  if(windowSize %% 2 == 0)
    windowSize <- windowSize - 1L
  context_start <- as.integer(max(CNEstart - (windowSize-1L)/2, 1))
  context_end <- as.integer(CNEend + (windowSize-1)/2)
  #win_nr_steps = windowSize / step_size
  #context_start = CNEstart - as.integer(((win_nr_steps-1)*step_size)/2+0.5)
  #if(context_start < 1)
  #  context_start = 1
  #context_end = CNEend + 
  #  as.integer(((win_nr_steps-1)*step_size)/2+step_size+0.5)
  ranges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                    context_start, context_end, 
                                    whichAssembly, minLength)
  # Implement get_cne_ranges_in_region_partitioned_by_other_chr later!!!
  ranges <- reduce(ranges)
  covAll <- coverage(ranges, width=context_end)
  runMeanAll <- runmean(covAll, k=windowSize, "constant")
  resStart <- max(CNEstart, (windowSize-1)/2+1)
  resEnd <- min(CNEend, context_end-(windowSize-1)/2)
  resCoords <- seq(resStart, resEnd, by=step_size)
  if(nrGraphs == 1){
    runMeanRes <- runMeanAll[resCoords]*100
    res <- cbind(resCoords, as.numeric(runMeanRes))
    colnames(res) <- c("coordinates", "y")
  }else{
    runMeanRes <- lapply(runMeanAll, "[", resCoords)
    runMeanRes <- lapply(runMeanRes, "*", 100)
    res <- list()
    for(i in 1:length(runMeanRes)){
      res[[names(runMeanRes)[i]]] <- cbind(resCoords, 
                                           as.numeric(runMeanRes[[i]]))
      colnames(res[[names(runMeanRes)[i]]]) <- c("coordinates", "y")
    }
  }
  return(res)
}


#calc_window_scores = function(CNEstart, CNEend, ranges, 
                        # win_nr_steps, step_size){
#  ## Here the starts and ends are 1-based.
#  CNElength = CNEend - CNEstart + 1
#  win_size = win_nr_steps * step_size
#  offsetBlk = as.integer(((win_nr_steps-1)*step_size)/2+0.5)
#  context_start = CNEstart - offsetBlk
#  if(context_start < 1)
#    context_start = 1
#  context_end = CNEend + offsetBlk
#  context_size = context_end - context_start + 1
#  #nr_blocks = as.integer(context_size/step_size) + 
   # ifelse(context_size%%step_size, 1, 0)
#  #blk_scores = numeric(ifelse(nr_blocks>win_nr_steps, 
  # nr_blocks, win_nr_steps+1))
#
#  covAll = coverage(ranges, width=context_end)
#   
#  #runMeanAll = runmean(covAll, k=windowSize, "constant")
#  #resStart = max(CNEstart, (windowSize-1)/2+1)
#  #resEnd = min(CNEend, computeEnd-(windowSize-1)/2)
#  #height = runMeanAll[resStart:resEnd]*100
#}

#listToPlot = list(a=res, b=res)

#plotCNE = function(listToPlot, horizonscale=2, nbands=3){
#  mergedDf = as.data.frame(do.call(rbind, listToPlot))
#  mergedDf$grouping = rep(names(listToPlot), sapply(listToPlot, nrow))
#  mergedDf = mergedDf[ ,c("coordinates", "grouping", "y")]
#  p = horizon.panel.ggplot(mergedDf, horizonscale=horizonscale, nbands=nbands)
#  #if(!is.null(file)){
#  #  postscript(file=file)
#  #  on.exit(dev.off())
#  #}
#  return(p)
#}

#horizon.panel.ggplot = function(mergedDf, horizonscale=2, 
#                                nbands=3, my.title="fun"){
#  #require(ggplot2)
#  #require(reshape2)
#  origin = 0
#  #require(RColorBrewer)
#  #col.brew = brewer.pal(name="RdBu",n=10)
#  #col.brew = c("#67001F", "#B2182B", "#D6604D", 
#                #"#F4A582", "#FDDBC7", "#D1E5F0", 
#                #"#92C5DE", "#4393C3", "#2166AC", "#053061")
#  col.brew = c("yellow", "orange", "red", "chartreuse", "blue")
#  colnames(mergedDf) = c("coordinates", "grouping", "y")
#  for(i in 1:nbands){
#    #do positive
#    mergedDf[ ,paste("ypos",i,sep="")] = 
#      ifelse(mergedDf$y > origin,
#             ifelse(abs(mergedDf$y) > horizonscale * i,
#                    horizonscale,
#                    ifelse(abs(mergedDf$y) - (horizonscale * (i - 1) - origin) 
#                           > origin, abs(mergedDf$y) - (horizonscale * (i - 1) 
#                                                        - origin), origin)),
#             origin)
#  }
#  mergedDf.melt = melt(mergedDf[,c(1,2,4:8)],id.vars=1:2)
#  colnames(mergedDf.melt) = c("coordinates","grouping","band","value")
#  p = ggplot(data=mergedDf.melt) +
#    geom_area(aes(x = coordinates, y = value, fill=band),
#                position="identity") +
#    scale_fill_manual(values=c("ypos1"=col.brew[1],
#                               "ypos2"=col.brew[2],
#                               "ypos3"=col.brew[3],
#                               "ypos4"=col.brew[4],
#                               "ypos5"=col.brew[5]))+
#    ylim(origin,horizonscale) +
#    facet_grid(grouping ~ .) +
#    theme_bw() +
#    theme(legend.position = "none",
#         strip.text.y = element_text(),
#         #axis.text.y = element_blank(), ## remove the y lables
#         axis.ticks = element_blank(),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         plot.title = element_text(size=16, face="bold", hjust=0))+
#    ggtitle(my.title)
#    return(p)
#}

#prepareCNETracks = function(dataMatrix, chr, strand, genome){
#  dTrack = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], 
#  data=dataMatrix[ ,2], chromosome=chr, strand=strand, genome=genome)

#}


