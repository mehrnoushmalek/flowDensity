setMethod(f="flowDensity",
          signature=c("flowFrame", "ANY", "logical", "missing"),
          definition=function(obj, channels, position, ...){
              .deGate2D(f=obj, channels=channels, position=position, ...)
          }
          )

setMethod(f="flowDensity",
	  signature=c("CellPopulation", "ANY", "logical", "missing"),
	  definition=function(obj, channels, position, ...){
              .deGate2D(f=obj@flow.frame, channels=channels, position=position, ...)
	  }
          )

setMethod(f="flowDensity",
	  signature=c("CellPopulation", "missing", "missing", "logical"),
	  definition=function(obj, singlet.gate){
              if(singlet.gate)
                  .sngltGate(obj)
              else
                  return(obj)
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

deGate <- function(flow.frame, channel, n.sd=1.5, use.percentile = FALSE, percentile = .95, use.upper=FALSE, upper = NA,talk=TRUE,
                   alpha = 0.1, sd.threshold = FALSE, graphs = FALSE, all.cuts = FALSE,tinypeak.removal=1/25, adjust.dens=1){

    ##========================================================================================================================================
    ## 1D density gating method
    ## Args:
    ##   flow.frame: a 'FlowFrame' object
    ##   channel: a channel's name or an integer to specify the channel
    ##   n.sd: an integer that is multiplied to the standard deviation to determine the place of threshold if 'sd.threshold' is 'TRUE'
    ##   use.percentile: if TRUE, forces to return the 'percentile'th threshold
    ##   percentile: a value in [0,1] that is used as the percentile. the default is 0.95.
    ##   graphs: if TRUE, it plots the threshold on the density curve
    ##   all.cuts: if TRUE, it returns all the cutoff points whose length+1 can roughly estimate the number of cell subsets in that dimension
    ##   upper: if TRUE, finds the change in the slope at the tail of the density curve, if FALSE, finds it at the head.
    ##   use.upper: if TRUE, forces to return the inflection point based on the first (last) peak if upper=F (upper=T)
    ##   tiny.peak.removal: a values in [0,1] for ignoring tiny peaks in density. Default is 1/25
    ##   adjust.dens: The smoothness of density, it is same as adjust in density(.). The default value is 1 and should not be changed unless necessary
    ## Value:
    ##   cutoffs, i.e. thresholds on the 1D data
    ## Author:
    ##   M. Jafar Taghiyar & Mehrnoush Malek
    ##-----------------------------------------------------------------------------------------------------------------------------------------
    .densityGating(flow.frame, channel, n.sd = n.sd, use.percentile = use.percentile, percentile = percentile, use.upper=use.upper,upper = upper,talk=talk,
                   alpha = alpha, sd.threshold = sd.threshold, graphs = graphs, all.cuts = all.cuts,tinypeak.removal=tinypeak.removal, adjust.dens=adjust.dens)
}


plotDens <- function(flow.frame, channels, col, main, xlab, ylab, pch = ".", s = FALSE, outdir, file.name, ...){

    ##===================================================
    ## Plot flowCytometry data with density-based color
    ##---------------------------------------------------
    f.exprs <- exprs(flow.frame)
    f.data <- pData(parameters(flow.frame))
    f.col.names <- colnames(flow.frame)

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
    if(s){
        if(!missing(file.name))
            fn <- paste(file.name, channels[2], channels[1], sep="_")
        else
            fn <- paste(channels[2], channels[1], sep="_")
        if(missing(outdir))
            outdir <- getwd()
        png(filename=paste(outdir,"/" ,fn, ".png",sep=""), width=800, height=800);
        plot(f.exprs[,channels], col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ...)
        dev.off()
    }else{
        plot(f.exprs[,channels], col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ...)
    }
}


notSubFrame <- function(flow.frame, channels, position = NA, gates, filter){

    ##===============================================================
    ## Excludes the subframe gated by 'gates' form given 'flow.frame'
    ##---------------------------------------------------------------
    if(missing(filter))
        .notSubFrame(flow.frame=flow.frame, channels=channels, position=position, gates=gates)
    else
        .notSubFrame(flow.frame=flow.frame, channels=channels, filter=filter)
}

nmRemove <- function(flow.frame, channels, neg=FALSE, talk=FALSE){
    ##===============================
    ## Removes negatives and margines and replaces them with NA
    ##-------------------------------
    new.f <- flow.frame
    for(chan in channels){
        if(is.character(chan))
            chan <- which(colnames(flow.frame)==chan)
        if(length(chan)==0)
            next
        d <- exprs(new.f)[,chan]
        m <- new.f@parameters@data$maxRange[chan]
        if(neg)
            i <- which(d > 0 & d != m)
        else
            i <- which(d != m)
        exprs(new.f)[i, ] <- exprs(new.f)[i,]
        if (talk)
          cat(length(which(!is.na(d)))-length(i), "Margin events removed from",colnames(new.f)[chan],"\n")
        exprs(new.f)[-i,] <- NA
    }
    return(new.f)
}


