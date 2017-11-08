##########################################################################################
#These wrapper functions are used to integerate flowDensity into openCyto framework
#####################################################################################

#' wrapper function used as a dispatcher to either 1d or 2d gating function
#' 
#' It can be registered as a plugin to openCyto framework
#' 
#' @param fr \code{flowFrame}
#' @param pp_res \code{list} the preprocesssing results returned from preprocessing module. (Not used)
#' @param xChannel, yChannel \code{character} specifies channels/markers to be used 
#'                                             for 1d gating, leave xChannel as \code{NA}
#' @param ... other arguments passed to flowDensity 
#' @return a \code{rectangleGate}(1d) or \code{polygonGate} (2d)
.flowDensity <- function(fr, pp_res, xChannel = NA, yChannel = NA, ...){
  
  chnls <- yChannel
  if(!is.na(xChannel)){
    chnls <- c(xChannel,yChannel)
  }
  
  if(length(chnls)==2)
    .flowDensity.2d(fr, channels = chnls, ...)
  else
    .flowDensity.1d(fr, channel = chnls, ...)
  
}

#' a wrapper for flowDensity::deGate
#' 
#' It mainly convert the cut point returned by deGate to a flowCore gate
#' 
#' @param fr \code{flowFrame}
#' @param channel \code{character} specifies channel/marker to be used for 1d gating,
#' @param positive \code{logical} indicate which side of cut point to keep. TRUE for right (+) , FALSE for left (-)
#' @param ... other arguments passed to \link{deGate}
#' @return a \code{rectangleGate}
.flowDensity.1d <- function(fr, channel, filterId = "", positive, ...){
  cutpoint <- flowDensity::deGate(fr, channel = channel, ...)

  #construct 1d filter 
  if (positive) {
    gate_coordinates <- list(c(cutpoint, Inf))
  } else {
    gate_coordinates <- list(c(-Inf, cutpoint))
  }
  
  names(gate_coordinates) <- channel
  
  rectangleGate(gate_coordinates, filterId = filterId)
  
}

#' lightweight version of flowDensity:::.deGate2D
#' 
#' We simply copy \link{.deGate2D} and modify the results so that 
#' it only returns polygon gate without carrying the flow data
#' and convert the matrix to a polygon filter
#' 
#' @return a \code{polygonGate}

.flowDensity.2d <- function (f, channels, position, n.sd=c(1.5,1.5), use.percentile=c(F,F), percentile=c(NA,NA), use.upper=c(F,F),
   upper=c(NA,NA), verbose=c(TRUE,TRUE), twin.factor=c(.98,.98), bimodal=c(F,F),filter=NA,tinypeak.removal=c(1/25,1/25),
   after.peak=c(NA,NA),alpha=c(0.1,0.1), sd.threshold=c(F,F), use.control=c(F,F), control=c(NA,NA),
   gates=c(NA,NA),all.cuts=c(F,F),node=NA,remove.margins=F,count.lim=5, ...)
{
  if (class(f)=="GatingHierarchy")
  {
    if(node!='root')
    {
      f <-flowWorkspace::getData(f,node)
    }else
    {
      warning("For gatingHierarchy objects, node is required, otherwise flowFrame at the root node will be used.")
      f <-flowWorkspace::getData(f,node)
    }
  }
  col.nm <- c(grep(colnames(f), pattern = "FSC"), grep(colnames(f), 
          pattern = "SSC"))
  if (length(col.nm) == 0) 
    warning("No forward/side scatter channels found, margin events not removed.")
  else f <- flowDensity:::nmRemove(f, col.nm)
  i <- which(!is.na(exprs(f)[, channels[1]]))
  if (length(i) == 0) 
    stop("invalid flowFrame input: all NA values")
  args <- names(list(...))
  eligible.args <- c('use.percentile', 'use.upper','upper', 'percentile', 'sd.threshold','filter',
                     'n.sd','twin.factor','bimodal','after.peak','use.control', 'control', 'alpha',
                     'scale', 'ellip.gate','tinypeak.removal','remove.margins','count.lim','verbose')
  if (length(setdiff(args, eligible.args) != 0)) 
    warning("unused argument(s): ", setdiff(args, eligible.args))
  col.nm <- c(grep(colnames(f), pattern="FSC"), grep(colnames(f), pattern="SSC"))
  if(remove.margins)
  {
    if(length(col.nm)==0)
      warning('No forward/side scatter channels found in control data for first channel, margin events not removed.')
    else
      f <- flowDensity:::nmRemove(f, col.nm)
  }
  if (is.na(gates[1]) & is.na(gates[2]) & is.na(position[1]) & 
      is.na(position[2])) 
    stop("Improper 'position' value, one position must be not 'NA' or 'gates' should be provided")
  if (is.matrix(filter))
  {
    inds <-.subFrame(f=f, channels, position=NA, gates, filter,include.equal=F)
    if(is.numeric(channels))
      channels <- c(colnames(f)[channels[1]], colnames(f)[channels[2]])
    cell.population  <- new("CellPopulation", flow.frame=f, channels=channels, position=c(NA,NA))
    exprs(cell.population@flow.frame)[-inds, ] <- NA
    cell.population@cell.count <- length(inds)
    cell.population@position <-position

    cell.population@proportion <- length(inds)/length(i)*100
    cell.population@filter <- filter
    cell.population@index<- inds
    g1<-ifelse(is.na(position[1]),yes = NA,no=ifelse(test = position[1],yes =min(filter[,1]),no = max(filter[,1]) ))
    g2<-ifelse(is.na(position[2]),yes = NA,no=ifelse(test = position[2],yes =min(filter[,2]),no = max(filter[,2]) ))
    cell.population@gates <- c(g1,g2)
    return(cell.population)
    
  }else{
  if (is.na(gates[1])) {
    if (use.control[1] & !is.na(position[1])) {
      if (is.na(control[1])) 
        stop("Missing 'control' data for first channel, set 'use.control' to FALSE if control should not be used")
      f.control1 <- NA
      control1.class <- class(control[1][[1]])
      if (control1.class == "flowFrame") 
        f.control1 <- control[1][[1]]
      else if (control1.class == "CellPopulation") 
        f.control1 <- control[1][[1]]@flow.frame
      else stop("The 'control' input for the first channel must be a flowFrame or CellPopulation object")
      col.nm.control1 <- c(grep(colnames(f.control1), pattern = "FSC"), 
          grep(colnames(f.control1), pattern = "SSC"))
      if (length(col.nm.control1) == 0) 
        warning("No forward/side scatter channels found in control data for first channel, margin events not removed.")
      else f.control1 <- flowDensity:::nmRemove(f.control1, col.nm.control1)
      gates[1] <- flowDensity:::.densityGating(flow.frame=f.control1, channel=channels[1], n.sd=n.sd[1], use.percentile=use.percentile[1], percentile=percentile[1], use.upper=use.upper[1], upper=upper[1],
                                 verbose=verbose[1],twin.factor=twin.factor[1],bimodal=bimodal[1],after.peak = after.peak[1],alpha=alpha[1],
                                 sd.threshold=sd.threshold[1],all.cuts=all.cuts[1],tinypeak.removal=tinypeak.removal[1],count.lim=count.lim)
    }
    else if(!is.na(position[1]))
      gates[1] <- flowDensity:::.densityGating(flow.frame=f, channel=channels[1], n.sd=n.sd[1], use.percentile=use.percentile[1], percentile=percentile[1], use.upper=use.upper[1], upper=upper[1],
                                 verbose=verbose[1],twin.factor=twin.factor[1],bimodal=bimodal[1],after.peak = after.peak[1],alpha=alpha[1],
                                 sd.threshold=sd.threshold[1],all.cuts=all.cuts[1],tinypeak.removal=tinypeak.removal[1],count.lim=count.lim)

    else gates[1] <- 0
  }
  if (is.na(gates[2])) {
    if (use.control[2] & !is.na(position[2])) {
      if (is.na(control[2])) 
        stop("Missing 'control' data for second channel, set 'use.control' to FALSE if control should not be used")
      f.control2 <- NA
      control2.class <- class(control[2][[1]])
      if (control2.class == "flowFrame") 
        f.control2 <- control[2][[1]]
      else if (control2.class == "CellPopulation") 
        f.control2 <- control[2][[1]]@flow.frame
      else stop("The 'control' input for the second channel must be a flowFrame or CellPopulation object")
      col.nm.control2 <- c(grep(colnames(f.control2), pattern = "FSC"), 
          grep(colnames(f.control2), pattern = "SSC"))
      if (length(col.nm.control2) == 0) 
        warning("No forward/side scatter channels found in control data for second channel, margin events not removed.")
      else f.control2 <- flowDensity:::nmRemove(f.control2, col.nm.control2)
    gates[2] <- flowDensity:::.densityGating(flow.frame=f.control2, channel=channels[2], n.sd=n.sd[2], use.percentile=use.percentile[2], percentile=percentile[2], use.upper=use.upper[2], upper=upper[2],
                               verbose=verbose[2],twin.factor=twin.factor[2],bimodal=bimodal[2],after.peak = after.peak[2],alpha=alpha[2],
                               sd.threshold=sd.threshold[2],all.cuts=all.cuts[2],tinypeak.removal=tinypeak.removal[2],count.lim=count.lim)
  }
    else if(!is.na(position[2]))
      gates[2] <- flowDensity:::.densityGating(flow.frame=f, channel=channels[2], n.sd=n.sd[2], use.percentile=use.percentile[2], percentile=percentile[2], use.upper=use.upper[2], upper=upper[2],
                                 verbose=verbose[2],twin.factor=twin.factor[2],bimodal=bimodal[2],after.peak = after.peak[2],alpha=alpha[2],
                                 sd.threshold=sd.threshold[2],all.cuts=all.cuts[2],tinypeak.removal=tinypeak.removal[2],count.lim=count.lim)
    else gates[2] <- 0
  }
  }
  new.f <- .ellipseGating(f, channels=channels, position=position, gates=gates, ...)
  cell.index <- which(!is.na(exprs(new.f)[,channels[1]]))
  cell.count <- length(cell.index)
  proportion <- cell.count/length(i) * 100
  if (cell.count != 0) {
    X <- exprs(new.f)[index, channels]
    filter <- chull(X)
    filter <- X[c(filter, filter[1]), ]

  }else {
    filter <- as.matrix(NA)
  }

  colnames(filter) <- channels
  polygonGate(.gate = filter)
  
}
