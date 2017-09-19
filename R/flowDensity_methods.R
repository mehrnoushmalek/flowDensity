setMethod(f="flowDensity",
          signature=c("flowFrame", "ANY", "logical", "missing"),
          definition=function(obj, channels, position, ...){
              .deGate2D(obj, channels=channels, position=position, ...)
          }
          )

setMethod(f="flowDensity",
	  signature=c("CellPopulation", "ANY", "logical", "missing"),
	  definition=function(obj, channels, position, ...){
              .deGate2D(obj@flow.frame, channels=channels, position=position, ...)
	  }
          )
setMethod(f="flowDensity",
          signature=c("GatingHierarchy", "ANY", "logical", "ANY"),
          definition=function(obj, channels, position, node,...){
            .deGate2D(obj, channels=channels, position=position,node=node, ...)
          }
)



setMethod(f="getflowFrame",
	  signature=c("CellPopulation"),
	  definition=function(obj){
              .getDataNoNA(obj)
	  }
          )

setMethod(f="plot", signature=c("flowFrame", "CellPopulation"),
          definition=function(x, y, ...){
              .deGatePlot(f=x, cell.population=y, ...)
          }
          )


deGate <- function(obj,channel, n.sd = 1.5, use.percentile = FALSE,  percentile =NA,use.upper=FALSE, upper = NA,verbose=TRUE,twin.factor=.98,
                   bimodal=F,after.peak=NA,alpha = 0.1, sd.threshold = FALSE, all.cuts = FALSE,
                   tinypeak.removal=1/25,node=NA, adjust.dens = 1,count.lim=20,magnitude=.3, ...){

    ##========================================================================================================================================
    ## 1D density gating method
    ## Args:
    ##   obj: a 'FlowFrame' object, 'CellPopulation' or 'GatingHierarchy'
    ##   channel: a channel's name or an integer to specify the channel
    ##   n.sd: an integer that is multiplied to the standard deviation to determine the place of threshold if 'sd.threshold' is 'TRUE'
    ##   use.percentile: if TRUE, forces to return the 'percentile'th threshold
    ##   percentile: a value in [0,1] that is used as the percentile. the default is NA. If set to a value(n) and use.percentile=F, it returns the n-th percentile, for 1-peak populations. 
    ##   use.upper: if TRUE, forces to return the inflection point based on the first (last) peak if upper=F (upper=T)
    ##   upper: if TRUE, finds the change in the slope at the tail of the density curve, if FALSE, finds it at the head.
    ##   verbose: When FALSE doesn't print any messages
    ##   twin.factor: reverse of tinypeak.removal for ignoring twinpeaks
    ##   bimodal: if TRUE, splits population closer to 50-50 when there are more than 2 peaks
    ##   after.peak: if TRUE (FALSE) and bimodal=FALSE, it chooses a cuttoff after(before) the maximum peak.
    ##   alpha: a value in [0,1) specifying the significance of change in the slope which would be detected
    ##   sd.threshold: if TRUE, it uses 'n.sd' times standard deviation for gating
    ##   all.cuts: if TRUE, it returns all the cutoff points whose length+1 can roughly estimate the number of cell subsets in that dimension
    ##   tiny.peak.removal: a values in [0,1] for ignoring tiny peaks in density. Default is 1/25
    ##   node: A character of the parent name to be excluded from gating hierarchy. Default is NA when x is flowFrame or CellPopulation object. Check getNodesof flowWorkspace
    ##   adjust.dens: The smoothness of density, it is same as adjust in density(.). The default value is 1 and should not be changed unless necessary
    ##   count.lim: minimum limit for event count in order to calculate the threshold. Default is 20
    ## Value:
    ##   cutoffs, i.e. thresholds on the 1D data
    ## Author:
    ##   M. Jafar Taghiyar & Mehrnoush Malek
    ##-----------------------------------------------------------------------------------------------------------------------------------------

  
  if (class(obj)=="GatingHierarchy")
      obj<-getData(obj,node)
    if (class(obj)=="CellPopulation")
      obj<-.getDataNoNA(obj)
    .densityGating(flow.frame = obj, channel, n.sd = n.sd, use.percentile = use.percentile, percentile = percentile, use.upper=use.upper,upper = upper,verbose=verbose,
                   twin.factor=twin.factor,bimodal=bimodal,after.peak=after.peak,alpha = alpha, sd.threshold = sd.threshold,all.cuts = all.cuts,
                   tinypeak.removal=tinypeak.removal, adjust.dens=adjust.dens,count.lim=count.lim,magnitude=magnitude, ...)
}

getPeaks <-  function(obj, channel,tinypeak.removal=1/25,verbose=F, adjust.dens=1,...){
 ##===================================================
  #Finding peaks for flowFrame objects or numeric vectors
 ##==================================================
 if(class(obj)=="numeric"|class(obj)=="vector"){
    x<-obj
    channel <-NA
  }else  if (class(obj)=="GatingHierarchy")
      {
      obj<-getData(obj,node)
      x <- exprs(obj)[, channel]
   }else if (class(obj)=="CellPopulation")
     { 
     obj<-.getDataNoNA(obj)
     x <- exprs(obj)[, channel]
     }else{
   x <- exprs(obj)[, channel]
  }
  n<- which(!is.na(x))
  if (length(n)< 3)
  {
    if(verbose)
      cat("Less than 3 cells, returning NA as a Peak.","\n")
    return(NA)
  }
  dens <- .densityForPlot(data = x, adjust.dens=adjust.dens,...)
  if(is.numeric(channel))
    channel <- colnames(obj)[channel]
  all.peaks <- .getPeaks(dens, peak.removal=tinypeak.removal)
  return(all.peaks)
}



plotDens <- function(obj, channels,node=NA ,col, main, xlab, ylab, pch = ".", ...){

    ##===================================================
    ## Plot flowCytometry data with density-based color
    ##---------------------------------------------------
  
      if (class(obj)=="GatingHierarchy")
      {
        if(!is.na(node))
        {
          flow.frame <-flowWorkspace::getData(obj,node)
        }else
        {
          warning("For gatingHierarchy objects, node is required, otherwise flowFrame at the root node will be used.")
          flow.frame <-flowWorkspace::getData(obj)
        }
      }
   if(class(obj)=="CellPopulation"){
        flow.frame<-obj@flow.frame
        
      }else{
    flow.frame<-obj
  }
    f.exprs <- exprs(flow.frame)
    na.ind <-which(is.na(f.exprs[,1]))
    if(length(na.ind)>0)
      f.exprs<-f.exprs[-na.ind,]
    f.data <- pData(parameters(flow.frame))
    f.col.names <- colnames(flow.frame)
    if (length(f.exprs[,channels[1]])<20)
        {
         col<-1
         pch<-20
        }
    if(missing(col)){
        colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
        col <- densCols(f.exprs[,channels], colramp = colPalette)
    }
    if(is.numeric(channels))
        channels <- f.col.names[channels]
    if(missing(xlab))
        xlab <- paste("<", channels[1], ">:", f.data$desc[which(f.col.names==channels[1])], sep = "")
    if(missing(ylab))
        ylab <- paste("<", channels[2], ">:", f.data$desc[which(f.col.names==channels[2])], sep = "")
    if(missing(main))
        main <- "All Events"
        plot(f.exprs[,channels], col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ...)
}


notSubFrame <- function(obj, channels, position = NA, gates, filter){

    ##===============================================================
    ## Excludes the subframe gated by 'gates' form given 'flow.frame'
    ##---------------------------------------------------------------
  if (class(obj)=="CellPopulation")
    flow.frame <- obj@flow.frame
  else if (class(obj)=="flowFrame")
    flow.frame <- obj
  else
    stop("Input should be either a cellPopulation or flowFrame.")
   if(missing(filter))
        .notSubFrame(flow.frame=flow.frame, channels=channels, position=position, gates=gates)
    else
          .notSubFrame(flow.frame=flow.frame, channels=channels, filter=filter)
}

nmRemove <- function(flow.frame, channels, neg=FALSE, verbose=FALSE,return.ind=FALSE){
    ##===============================
    ## Removes negatives and margines and replaces them with NA
    ##-------------------------------
  
    inds<-vector("list",length(channels))
    names(inds)<-channels
    new.f <- flow.frame
    for(chan in channels){
       j<- ifelse(is.character(chan),yes= which(colnames(flow.frame)==chan),no=chan)
        if(length(j)==0)
            next
        d <- exprs(new.f)[,j]
        m <- new.f@parameters@data$maxRange[j]
        if(neg)
            i <- which(d > 0 & d != m)
        else
            i <- which(d != m)
        exprs(new.f)[i, ] <- exprs(new.f)[i,]
        if (verbose)
          cat(length(which(!is.na(d)))-length(i), "Margin events removed from",colnames(new.f)[j],"\n")
        exprs(new.f)[-i,] <- NA
        inds[[chan]]<-which( is.na(exprs(new.f)[,1]))
    }
    if(return.ind)
       return(inds)
    exprs(new.f) <- exprs(new.f)[which(!is.na(exprs(new.f))[,1]),]
    return(new.f)
}


