.densityGating <- function(f, channel, n.sd = 1.5, use.percentile = FALSE,  percentile = 0.95,use.upper=FALSE, upper = NA, avg = FALSE,talk=TRUE,
                           alpha = 0.1, sd.threshold = FALSE, graphs = FALSE, all.cuts = FALSE, debris.gate = FALSE,tinypeak.removal=tinypeak.removal, adjust.dens = 1){
  
  ##========================================================================================================================================
  ## 1D density gating method
  ## Args:
  ##   flow.frame: a 'FlowFrame' object
  ##   channel: a channel's name or an integer to specify the channel
  ##   n.sd: an integer that is multiplied to the standard deviation to determine the place of threshold if 'sd.threshold' is 'TRUE'
  ##   use.percentile: if TRUE, forces to return the 'percentile'th threshold
  ##   use.upper: if True, forces to return the inflection point based on the first (last) peak if upper=F (upper=T)
  ##   percentile: a value in [0,1] that is used as the percentile. The default value is 0.95.
  ##   upper: if 'TRUE', it finds the change in the slope after the peak with index 'peak.ind'
  ##   avg: if TRUE, the mean of all identified cutoff points are used
  ##   alpha: a value in [0,1) specifying the significance of change in the slope which would be detected
  ##   sd.threshold: if TRUE, it uses 'n.sd' times standard deviation for gating
  ##   graphs: if TRUE, it plots the density as well as the threshold on the same plot
  ##   all.cuts: if TRUE, it returns all the cutoff points whose length+1 can roughly estimate the number of cell subsets in that dimension
  ##   debris.gate: if TRUE, it would try to remove the debris and gate the lymphocyte cells
  ##   tiny.peak.removal: a values in [0,1] for ignoring tiny peaks in density
  ##   adjust.dens: The smoothness of density, it is same as adjust in density(.).The default value is 1 and should not be changed unless necessary
  ## Value:
  ##   cutoffs, i.e. thresholds on the 1D data
  ## Authors:
  ##   M. Jafar Taghiyar & Mehrnoush Malek
  ##-----------------------------------------------------------------------------------------------------------------------------------------
  x <- exprs(f)[, channel]
  if (length(x)<2)
  {
    print("Less than 2 cells, returning NA as a threshold.")
    return(NA)
  }
  dens <- .densityForPlot(x, adjust.dens)
  stdev <- sd(x)
  med <- median(dens$x)
  if(is.numeric(channel))
    channel <- colnames(f)[channel]
  cutoffs <- c()
  peaks <- .getPeaks(dens, tinypeak.removal=tinypeak.removal)
  peak.ind <- peaks[[2]]
  peaks <- peaks[[1]]
  len <- length(peaks)
  for(i in 1:(len-1))
    cutoffs[i] <- .getIntersect(dens, peaks[i], peaks[i+1])
  if(all.cuts & length(peaks)>1)
    return(cutoffs)
  if(avg)
    return(mean(cutoffs))
  if(use.percentile)
    cutoffs <- quantile(x, percentile, na.rm=T)
  else if(debris.gate) # Removes debris and gate Lymphocytes
    cutoffs <- max(quantile(x, 0.1, na.rm=T), min(cutoffs))
  else if(use.upper){
    peak.ind <- ifelse(upper, peak.ind[len], peak.ind[1])
    track.slope.cutoff <- .trackSlope(dens=dens, peak.ind=peak.ind, upper=upper, alpha=alpha)
    if(is.infinite(track.slope.cutoff))
      cutoffs <- ifelse(upper, max(dens$x), min(dens$x))
    else
      cutoffs <- track.slope.cutoff
  }
  else if(len==1){
    if(talk)
      print("Only one peak is found, set standard deviation, percentile, or upper arguments accordingly")
    if(!is.na(percentile)){
      cutoffs <- quantile(x, percentile, na.rm=T)
    } else if (!is.na(upper))
    {
      track.slope.cutoff <- .trackSlope(dens=dens, peak.ind=peak.ind, upper=upper, alpha=alpha)
      if(is.infinite(track.slope.cutoff))
        cutoffs <- ifelse(upper, max(dens$x), min(dens$x))
      else
        cutoffs <- track.slope.cutoff
    }else if (sd.threshold)
    {
      upper <- as.logical(ifelse(peaks>med, FALSE, TRUE))
      cutoffs <- ifelse(upper, peaks+n.sd*stdev, peaks-n.sd*stdev)
    }
    else{
      flex.point <- .getFlex(dens=dens, peak.ind=peak.ind)
      if(!is.na(flex.point))
      { 
        if(talk)
          print("Cutoff is based on inflection point.")
        cutoffs <- flex.point
      }else{
        if(is.na(upper))
          upper <- as.logical(ifelse(peaks>med, FALSE, TRUE))
        track.slope.cutoff <- .trackSlope(dens=dens, peak.ind=peak.ind, upper=upper, alpha=alpha)
        stdev.cutoff <- ifelse(upper, peaks+n.sd*stdev, peaks-n.sd*stdev)
        cutoffs <- ifelse(is.infinite(track.slope.cutoff)|sd.threshold, stdev.cutoff, track.slope.cutoff) # use the 's.threshold' if the gate from 'trackSlope()' is too loose
        if(talk)
          print(" cuttoff is based on tracking the slope and cpmparing the position of the peak and mean of the population.")
      }
    }
  }
  else{
    distance <- .getDist(peaks)
    index <- .getScoreIndex(dens, cutoffs, distance)
    cutoffs <- cutoffs[index]
  }
  
  if(graphs){
    plot(dens, type="l", main=paste(channel, pData(parameters(f))$desc[which(colnames(f)==channel)], sep=": "))
    abline(v=cutoffs, lty="dashed", lwd=2, col=2)
  }
  return(cutoffs)
}

.deGate2D <- function(f, channels, position, n.sd=c(1.5,1.5), use.percentile=c(F,F), use.upper=c(F,F), percentile=c(NA,NA),
                      upper=c(NA,NA), avg=c(F,F), alpha=c(0.1,0.1), sd.threshold=c(F,F), use.control=c(F,F), control=c(NA,NA),
                      debris.gate=c(F,F), gates=c(NA,NA),all.cuts=c(F,F),remove.neg=F, ...){
  
  ##=========================================================================================================================
  ## 2D density gating method
  ##  Args:
  ##   f: a 'FlowFrame' object
  ##   channels: a vector of two channel names or their equivalent integer number
  ##   position: a vector of two logical values specifying the position of the cell subset of interest on the 2D plot
  ##   use.percentile, use.upper, upper, percentile, sd.threshold, n.sd, alpha: refer to the arguments of the '.densityGating(.)' function
  ##   ...: provides access to the arguments of the '.ellipseGating(.)' function
  ## Value:
  ##   a 'CellPopulation' object that includes the result of 2D gating
  ## Author:
  ##   M. Jafar Taghiyar
  ##--------------------------------------------------------------------------------------------------------------------------
  col.nm <- c(grep(colnames(f), pattern="FSC"), grep(colnames(f), pattern="SSC"))
  if(length(col.nm)==0)
    warning('No forward/side scatter channels found, margin events not removed.')
  else
    f <- nmRemove(f, col.nm,neg = remove.neg)
  i <- which(!is.na(exprs(f)[,channels[1]]))
  if(length(i)==0)
    stop('invalid flowFrame input: This flowFrame has 0 cells')
  if(length(i)<3)
  {
    if(is.numeric(channels))
      channels <- c(colnames(f)[channels[1]], colnames(f)[channels[2]])
    cell.population <- new("CellPopulation",
                           flow.frame=f,
                           proportion=100,
                           cell.count=length(i),
                           channels=channels,
                           position=position,
                           gates=c(0,0),
                           filter=as.matrix(NA),
                           index=i)
    return(cell.population)
  }
  args <- names(list(...))
  eligible.args <- c('use.percentile', 'use.upper','upper', 'avg', 'percentile', 'sd.threshold', 'n.sd',
                     'use.control', 'control', 'alpha', 'debris.gate', 'scale', 'ellip.gate', 'graphs')
  if(length(setdiff(args, eligible.args)!=0))
    warning('unused argument(s): ', setdiff(args, eligible.args))
  
  if(is.na(gates[1])&is.na(gates[2])&is.na(position[1])&is.na(position[2]))
    stop("Improper 'position' value, one position must be not 'NA' or 'gates' should be provided" )
  if(is.na(gates[1])){
    if(use.control[1]&!is.na(position[1])){
      if(is.na(control[1]))
        stop("Missing 'control' data for first channel, set 'use.control' to FALSE if control should not be used")
      
      f.control1 <- NA
      control1.class <- class(control[1][[1]])
      if(control1.class=="flowFrame")
        f.control1 <- control[1][[1]]
      else if(control1.class=="CellPopulation")
        f.control1 <- control[1][[1]]@flow.frame
      else
        stop("The 'control' input for the first channel must be a flowFrame or CellPopulation object")
      
      col.nm.control1 <- c(grep(colnames(f.control1), pattern="FSC"), grep(colnames(f.control1), pattern="SSC"))
      if(length(col.nm.control1)==0)
        warning('No forward/side scatter channels found in control data for first channel, margin events not removed.')
      else
        f.control1 <- nmRemove(f.control1, col.nm.control1)
      
      gates[1] <- .densityGating(f=f.control1, channel=channels[1], use.percentile=use.percentile[1], use.upper=use.upper[1], percentile=percentile[1], upper=upper[1], avg=avg[1],
                                 all.cuts=all.cuts[1],sd.threshold=sd.threshold[1], n.sd=n.sd[1], alpha=alpha[1], debris.gate=debris.gate[1],tinypeak.removal=1/25)
    }
    else if(!is.na(position[1]))
      gates[1] <- .densityGating(f=f, channel=channels[1], use.percentile=use.percentile[1], percentile=percentile[1],use.upper=use.upper[1], upper=upper[1], avg=avg[1],
                                 all.cuts=all.cuts[1],sd.threshold=sd.threshold[1], n.sd=n.sd[1], alpha=alpha[1], debris.gate=debris.gate[1],tinypeak.removal=1/25)
    else
      gates[1] <- 0
  }
  if(is.na(gates[2])){
    if(use.control[2]&!is.na(position[2])){
      if(is.na(control[2]))
        stop("Missing 'control' data for second channel, set 'use.control' to FALSE if control should not be used")
      
      f.control2 <- NA
      control2.class <- class(control[2][[1]])
      if(control2.class=="flowFrame")
        f.control2 <- control[2][[1]]
      else if(control2.class=="CellPopulation")
        f.control2 <- control[2][[1]]@flow.frame
      else
        stop("The 'control' input for the second channel must be a flowFrame or CellPopulation object")
      
      col.nm.control2 <- c(grep(colnames(f.control2), pattern="FSC"), grep(colnames(f.control2), pattern="SSC"))
      if(length(col.nm.control2)==0)
        warning('No forward/side scatter channels found in control data for second channel, margin events not removed.')
      else
        f.control2 <- nmRemove(f.control2, col.nm.control2)
      
      gates[2] <- .densityGating(f=f.control2, channel=channels[2], use.percentile=use.percentile[2], use.upper=use.upper[2],percentile=percentile[2], upper=upper[2], avg=avg[2],
                                 all.cuts=all.cuts[2],sd.threshold=sd.threshold[2], n.sd=n.sd[1], alpha=alpha[1], debris.gate=debris.gate[2],tinypeak.removal=1/25)
    }
    else if(!is.na(position[2]))
      gates[2] <- .densityGating(f=f, channel=channels[2], use.percentile=use.percentile[2], percentile=percentile[2],use.upper=use.upper[2], upper=upper[2], avg=avg[2],
                                 all.cuts=all.cuts[2],sd.threshold=sd.threshold[2], n.sd=n.sd[1], alpha=alpha[1], debris.gate=debris.gate[2],tinypeak.removal=1/25)
    else
      gates[2] <- 0
  }
  
  new.f <- .ellipseGating(flow.frame=f, channels=channels, position=position, gates=gates, ...)
  index <- which(!is.na(exprs(new.f)[,channels[1]]))
  cell.count <- length(index)
  proportion <- cell.count/length(i) * 100
  
  if(cell.count>1){
    X <- exprs(new.f)[index,channels]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
  } else 
    filter <- as.matrix(NA)
  
  if(is.numeric(channels))
    channels <- c(colnames(f)[channels[1]], colnames(f)[channels[2]])
  
  cell.population <- new("CellPopulation",
                         flow.frame=new.f,
                         proportion=proportion,
                         cell.count=cell.count,
                         channels=channels,
                         position=position,
                         gates=gates,
                         filter=filter,
                         index=index
  )
  return(cell.population)
}

.deGatePlot <- function(f, cell.population, adjust.dens = 1,...){
  
  ##========================================================================================
  ## Plots the output of 2D density-estimate gating method
  ##  Args:
  ##   f: a 'flowFrame' object
  ##   cell.population: output of 'flowDensity(.)',i.e. an object of class 'CellPopulation'
  ## Value:
  ##   a plot showing the density-estimate gating results
  ## Author:
  ##   M. Jafar Taghiyar
  ##----------------------------------------------------------------------------------------
  channels <- cell.population@channels
  pos <- cell.population@position
  gates <- cell.population@gates
  f.col.names <- colnames(f)
  f.data <- pData(parameters(f))
  exprs(f) <- exprs(f)[which(!is.na(exprs(f)[,1])),]
  
  #    require(gplots)
  #    dev.new(pointsize=18)
  .layoutCreate(1,1)
  par(mar=rep(0,4))
  par(yaxt="n", cex.axis=1.25, bty="n", ann=F)
  
  ## Display statistics of the cell population
  obj <- sprintf("Cell count: %d\nProportion: %.2f", cell.population@cell.count, cell.population@proportion)
  textplot(obj, cex=1.5)
  
  ## Plot the pdf of channels[1]
  data.chan1 <- exprs(f)[,channels[1]]
  dens.chan1 <- .densityForPlot(data.chan1, adjust.dens)
  plot(dens.chan1, xlim=range(data.chan1), type="l", frame.plot=F, axes=F, ann=F)
  
  ## Fill the pdf of channels[1]
  pts <- .densRange(dens.chan1$x, dens.chan1$y, gates[1], pos[1])
  polygon(pts$x, pts$y, col="#00ff0044", border= "#00ff00ff")
  
  ## Add the bottom line and the density-estimated gate
  abline(h=0)
  if(!is.na(pos[1]))
    abline(v=gates[1], lty="dashed", col=colors()[290], lwd=1.5)
  
  ## Plot the pdf of channels[2]
  data.chan2 <- exprs(f)[,channels[2]]
  dens.chan2 <- .densityForPlot(data.chan2, adjust.dens)
  plot(dens.chan2$y, dens.chan2$x, ylim=range(data.chan2), xlim=rev(range(dens.chan2$y)),
       type="l", col=1, frame.plot=F, axes=F, ann=F)
  
  ## Fill the pdf of channels[2]
  pts <- .densRange(dens.chan2$x, dens.chan2$y, gates[2], pos[2])
  polygon(pts$y, pts$x, col="#ff450055", border="#ff4500ff")
  
  ## Add the bottom line and the density-estimated gate
  abline(v=0)
  if(!is.na(pos[2]))
    abline(h=gates[2], lty="dashed", col=colors()[290], lwd=1.5)
  
  ## Plot the scatter plot of the flowframe
  plotDens(flow.frame=f, channels=channels, devn=F, ...)
  
  ## Add labels of the axes
  index1 <- which(f.col.names==channels[1])
  index2 <- which(f.col.names==channels[2])
  text.chan1 <- f.data$desc[index1]
  text.chan2 <- f.data$desc[index2]
  if(is.na(text.chan1))
    text.chan1 <- channels[1]
  if(is.na(text.chan2))
    text.chan2 <- channels[2]
  mtext(text=text.chan1, line=3, side=1, font=2, cex=2)
  mtext(text=text.chan2, line=5, side=2, font=2, cex=2)
  
  ## Add the gates
  if(!is.na(pos[1]))
    abline(v=gates[1], lty="dashed", col=colors()[290], lwd=1.5)
  if(!is.na(pos[2]))
    abline(h=gates[2], lty="dashed", col=colors()[290], lwd=1.5)
  
}

.layoutCreate <- function(m, n, l.show = FALSE){
  
  ##=========================================================
  ## Create suitable layout for the '.deGatePlot(.)' function
  ##---------------------------------------------------------
  mat <- c()
  tmp <- c()
  x <- matrix(c(1,3,3,0,2,4,4,0,2,4,4,0), nrow=4)
  for(i in 1:n){
    tmp <- cbind(tmp,x+(i-1)*4)
  }
  for(j in 1:m){
    mat <- rbind(mat,tmp+(j-1)*4*n)
  }
  for(i in seq(4,4*m, 4))
    mat[i,] <- 0;
  l <- layout(mat)
  if(l.show)
    layout.show(l)
}

.getPeaks <- function(dens,tinypeak.removal=tinypeak.removal){
  
  ##=====================================================================================================================
  ## Finds the peaks in the given density
  ## Args:
  ##   dens: density of the channel whose peaks are going to be found. It is a variable of class 'density'.
  ##   w: the length of the window where the function searches for the peaks. If all peaks required, use the default w=1.
  ## Value:
  ##   peaks in the density of the provided channel
  ##---------------------------------------------------------------------------------------------------------------------
  d <- dens$y
  w <- 1
  peaks <- c()
  peaks.ind <- c()
  for(i in 1:(length(d)-w)){
    if(d[i+w] > d[(i+w+1):(i+2*w)] && d[i+w] > d[i:(i+w-1)] && d[i+w] > tinypeak.removal*max(d)){ # also removes tiny artificial peaks less than ~%4 of the max peak
      peaks <- c(peaks, dens$x[i+w])
      peaks.ind <- c(peaks.ind, i+w)
    }
  }
  return(list(Peaks=peaks, P.ind=peaks.ind))
}

.getIntersect <- function(dens, p1, p2){
  
  ##=================================================================
  ## Returns the min intersection between two peaks
  ## Args:
  ##   dens: density of the channel whose peaks are going to be found
  ##   p1: firs peak
  ##   p2: second peak
  ## Value:
  ##   the min intersection between two peaks
  ##-----------------------------------------------------------------
  if(missing(p2) | is.na(p2))
    p2 <- min(dens$x)
  fun <- splinefun(dens)
  index <- intersect(which(dens$x < max(c(p1,p2))), which(dens$x > min(c(p1,p2))))
  min.intsct <- dens$x[index][which.min(fun(dens$x[index]))]
  #if(length(min.intsct==0) | is.infinite(min.intsct))
  #  min.intsct <- min(dens$x)
  return(min.intsct)
}

.getDist <- function(peaks){
  
  ##============================================
  ## Returns the distance between adjacent peaks
  ## Args:
  ##   peaks: peaks of the density
  ## Value:
  ##   distance between adjacent peaks
  ##--------------------------------------------
  len <- length(peaks)-1
  d <- vector('numeric', len)
  for(i in 1:len)
    d[i] <- abs(peaks[i] - peaks[i+1])
  return(d)
}

.getScoreIndex <- function(dens, cutoffs, d){
  
  ##==========================================================================================
  ## Returns the score value based upon which the .densityGating() function decides on the thresholds
  ## Args:
  ##   dens: density of the channel to be investigated
  ##   cutoffs: the min intersections returned by getIntersect() function
  ##   d: the distance between adjacent peaks returned by the getDist() function
  ## Value:
  ##   the index of the corresponding metric value
  ##------------------------------------------------------------------------------------------
  h <- c()
  fun <- splinefun(dens)
  h <- fun(cutoffs);
  max.metric <- 0;
  max.h <- max(h);
  max.d <- max(d);
  if(length(which(h==0))==0)
    return(which.max(d/h))
  else
    return(which.max(d[which(h==0)]))
  
}

.trackSlope <- function(dens, peak.ind, alpha, upper=T, w=10, return.slope=F){
  
  ##==========================================================================================================================
  ## Returns the point of the large change in the slope of the density, i.e., 'dens', where the threshold is likely to be there
  ## Args:
  ##   dens: density of the channel to investigate
  ##   peak.ind: index of the peak in the density distribution of the channel
  ##   alpha: a value in [0,1) specifying the significance of change in the slope which would be detected
  ##   upper: if 'TRUE', it finds the change in the slope after the peak with index 'peak.ind'
  ##   w: specifies the length of the window within which the change of slope is investigated
  ##   return.slope: if 'TRUE' returns all the slope values only
  ## Value:
  ##   the point where the slope changes dramatically
  ## Author:
  ##   M. Jafar Taghiyar
  ##---------------------------------------------------------------------------------------------------------------------------
  d.y <- dens$y
  d.x <- dens$x
  slope <- c()
  heights<- c()
  lo <- ifelse(upper, peak.ind+w, 1)
  up <- ifelse(upper, length(d.y)-w, peak.ind-w)
  
  for(i in seq(from=lo,to=up,by=w)){
    slope <- c(slope, abs((d.y[i+w]-d.y[i])/(d.x[i+w]-d.x[i])))
    heights <-c(heights,d.y[i])
  }
  ## try to estimate the optimum alpha using the distribution of slopes
  #   if(missing(alpha))
  #     alpha <- ifelse(sd(slope) < max(slope)/10, 0.01, 0.1)
  if(return.slope)
    return(slope)
  small.slopes <- which(slope < (max(slope)*alpha))
  small.slopes <- small.slopes[which(heights[small.slopes] < .3*d.y[peak.ind])]
  len <- length(small.slopes)
  if(len==0)
    return(ifelse(upper, Inf, -Inf))
  ind <- ifelse(upper, small.slopes[1], small.slopes[len]) #if upper, the first large change of slope is given else the last change of slope is returned
  
  return(ifelse(upper, d.x[peak.ind+w*(ind-1)], d.x[ind*w]))
}

.getFlex <- function(dens, peak.ind, w=3){
  
  ##==========================================================================================================================
  ## Returns the point of inflection or flex
  ## Args:
  ##   dens: density of the channel to investigate
  ##   peak.ind: index of the peak in the density distribution of the channel
  ##   w: specifies the length of the window within which the change of slope is investigated
  ## Value:
  ##   the point of inflection
  ## Author:
  ##   M. Jafar Taghiyar
  ##---------------------------------------------------------------------------------------------------------------------------
  upper.slope <- .trackSlope(dens=dens, peak.ind=peak.ind, alpha=0.1, upper=T, w=w, return.slope=T)
  lower.slope <- .trackSlope(dens=dens, peak.ind=peak.ind, alpha=0.1, upper=F, w=w, return.slope=T)
  for(i in 2: (length(upper.slope)-1)){
    if(upper.slope[i]<upper.slope[i+1] & upper.slope[i]<upper.slope[i-1]){
      upper.ind <- (i-1)*w
      break
    }else
      upper.ind <- Inf
  }
  for(i in seq(from=length(lower.slope)-1, to=2, by=-1)){
    if(lower.slope[i]<lower.slope[i+1] & lower.slope[i]<lower.slope[i-1]){
      lower.ind <- (i-1)*w
      break
    }else
      lower.ind <- Inf
  }
  if(is.infinite(upper.ind)&is.infinite(lower.ind))
    return(NA)
  if(upper.ind<lower.ind)
    return(dens$x[peak.ind+upper.ind])
  else
    return(dens$x[peak.ind-lower.ind])
}

.densRange <- function(x, y, gate, pos = FALSE){
  
  ##==================================================================
  ## Plots the output of density-estimate gating method
  ##  Args:
  ##   x: 'x' slot of the density returned by the 'density()' function
  ##   y: 'y' slot of the density returned by the 'density()' function
  ##   gate: a threshold given by the '.densityGating()' function
  ##   pos: refer to the 'pos' in 'deGatePlot()' function
  ## Value:
  ##   (x,y) coordinates
  ##------------------------------------------------------------------
  pts <- list()
  if(is.na(pos))
    return(list(x=c(x,tail(x,1),x[1]), y=c(y,min(y),min(y))))
  if(pos){
    x.pts <- c(x[which(x>=gate)], tail(x[which(x>=gate)],1), x[which(x>=gate)][1])
    y.pts <- c( y[which(x>=gate)], min(y[which(x>=gate)]), min(y[which(x>=gate)]))
  }else{
    x.pts <- c(x[which(x<gate)][1], x[which(x<gate)], gate)
    y.pts <- c(min(y[which(x>=gate)]),y[which(x<gate)], min(y[which(x<gate)]))
  }
  pts$x <- x.pts
  pts$y <- y.pts
  return(pts)
}

.densityForPlot <- function(data, adjust.dens=1){
  dens <- density(data[which(!is.na(data))], adjust=adjust.dens)
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

.subFrame <- function(f, channels, position, gates, filter){
  
  ##================================================================================================================================
  ## Returns the new flowframe out of flowFrame 'f' gated by the '.densityGating()' thresholds on 'channels'
  ## Args:
  ##   f: a 'FlowFrame' object
  ##   channels: a vector specifying which 2 channels used for 2D gating, eg. c(chan1, chan2)
  ##   pos: a vector of two logical values that determines positive/negative for each channel, i.e. the position of the ellipse gate
  ##   gates: a two-element vector specifying the thresholds for each channel, usually the output of the .densityGating(.) function
  ##   filter: the boundary of a 'CellPopulation' object. Its value is stored in the 'filter' slot of the object
  ## Value:
  ##   indices of the subset of input flowframe
  ##--------------------------------------------------------------------------------------------------------------------------------
  f.exprs <- exprs(f)
  if(!missing(filter)){
    POK <- list(x=filter[,1], y=filter[,2])
    ind <- which(is.na(f.exprs[,channels[1]]))
    if(length(ind) != 0)
      f.exprs <- f.exprs[-ind, ]
    in.p <-  inpoly(x=f.exprs[,channels[1]], y=f.exprs[,channels[2]], POK=POK)
    tmp <- vector(mode='numeric', length=length(exprs(f)[,channels[1]]))
    if(length(ind) != 0)
      tmp[-ind] <- in.p
    else
      tmp <- in.p
    f.sub.ind <- which(tmp == 1)
    return(f.sub.ind)
  }
  else if(is.na(position[1])&is.na(position[2]))
    stop("Improper 'position' value, one position must be not 'NA'")
  if(is.na(position[1])){
    if(position[2])
      f.sub.ind <- which(f.exprs[,channels[2]] > gates[2])
    else
      f.sub.ind <- which(f.exprs[,channels[2]] < gates[2])
  }
  else if(is.na(position[2])){
    if(position[1])
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1])
    else
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1])
  }
  else{
    if(position[1]&position[2])
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1] & f.exprs[,channels[2]] > gates[2])
    else if (position[1]&!position[2])
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1] & f.exprs[,channels[2]] < gates[2])
    else if (!position[1]&position[2])
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1] & f.exprs[,channels[2]] > gates[2])
    else
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1] & f.exprs[,channels[2]] < gates[2])
  }
  return(f.sub.ind)
}

.notSubFrame <- function(flow.frame, channels, position = NA, gates, filter){
  
  ##===============================================================
  ## Excludes the subframe gated by 'gates' form given 'flow.frame'
  ## Args:
  ##   please refer to the '.subFrame(.)' function
  ## Value:
  ##   a 'CellPopulation' object from which the undesirable cell population
  ##   has been removed
  ## Author:
  ## M. Jafar Taghiyar
  ##---------------------------------------------------------------
  
  if(!missing(filter))
    ind <- .subFrame(f=flow.frame, channels=channels, filter=filter)
  else
    ind <- .subFrame(f=flow.frame, channels=channels, position=position, gates=gates)
  
  if(is.numeric(channels))
    channels <- c(colnames(flow.frame)[channels[1]], colnames(flow.frame)[channels[2]])
  
  cell <- new("CellPopulation", flow.frame=flow.frame, channels=channels, position=position)
  exprs(cell@flow.frame)[ind, ] <- NA
  index <- which(!is.na(exprs(cell@flow.frame)[,channels[1]]))
  cell@cell.count <- length(index)
  cell@proportion <- length(index)/length(exprs(flow.frame)[,channels[1]]) * 100
  if(length(index)>1){
    X <- exprs(cell@flow.frame)[index,channels]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
  }else
    filter <- as.matrix(NA)
  cell@filter <- filter
  cell@index <- index
  if(!missing(gates))
    cell@gates <- gates
  return(cell)
}

.getDataNoNA <- function(obj){
  
  ##==========================================================
  ## Returns the flow date of a given 'CellPopulation' object
  ## Args:
  ##   cell.population: an object of class 'CellPopulation'
  ## Value:
  ##   a flowFrame
  ##----------------------------------------------------------
  f <- obj@flow.frame
  ind <- which(!is.na(exprs(f))[,1])
  if (length(ind)==1)
  {
    print("Cannot make frames with 1 cell, duplicating the cell to make a flowFrame")
    exprs(f) <- rbind(exprs(f)[ind,],exprs(f)[ind,])
  }else if(length(ind)==0){
    print("Cannot make frames with 0 cell, making a flowFrame of size 2, all values euqal to -10")
    temp <- rbind(rep(-10,ncol(f)),rep(-10,ncol(f)))
    colnames(temp) <- colnames(exprs(f))
    exprs(f) <- temp
  }else{
    
  exprs(f) <- exprs(f)[which(!is.na(exprs(f))[,1]),]
  }
  return(f)

}

.sngltGate <- function(cell.population){
  
  ##==========================================================
  ## Returns singlets cell subset of a given 'flowFrame' object
  ## Args:
  ##   flow.frame: an object of class 'flow.frame'
  ## Value:
  ##   singlet cells, an object of class 'CellPopulation'
  ##----------------------------------------------------------
  flow.frame <- cell.population@flow.frame
  col.names <- colnames(flow.frame)
  chans <- c('FSC-A', 'SSC-A')
  if('FSC-H' %in% col.names)
    chans <- c(chans[1], 'FSC-H')
  else if('FSC-W' %in% col.names)
    chans <- c(chans[1], 'FSC-W')
  else if('SSC-H' %in% col.names)
    chans <- c(chans[2], 'SSC-H')
  else if('SSC-W' %in% col.names)
    chans <- c(chans[2], 'SSC-W')
  else{
    warning('None of these channels exist: "FSC-H", "FSC-W", "SSC-H", "SSC-W" \nThe same frame was returned.' )
    return(cell.population)
  }
  snglt <- flowDensity(obj=flow.frame, channels=chans, position=c(F,NA),
                       use.percentile=c(T,F), percentile=c(0.99,NA), ellip.gate=T, scale=0.99)
  #  sf <- snglt@filter
  #  x.mean <- (min(sf[,1])+max(sf[,1]))/2
  #  ind1 <- which.min(abs(sf[,1]-x.mean))
  
  #  ## Make criterion
  #  calcSlope <- function(p1,p2) (p1[2]-p2[2])/(p1[1]-p2[1])
  #  p.max <- c(max(sf[,1]),max(sf[,2]))
  #  p.min <- c(min(sf[,1]),min(sf[,2]))
  #  m <- calcSlope(p.max,p.min)
  #  crt1 <- .luline(point=p.max, m=m, up=T, tp=sf[ind1,])
  
  #  ## Find the ind2 based on crt1
  #  i <- ind1
  #  a <- sf
  #  repeat{
  #    a <- a[-i,]
  #    if(length(a[,1])!=0){
  #      ind2 <- which.min(abs(a[,1]-x.mean))
  #      crt2 <- .luline(point=p.max, m=m, up=T, tp=a[ind2,])
  #    }
  #    if(as.logical(ifelse(crt1, !crt2, crt2)) | length(a[,1])==0)
  #      break()
  #    i <- ind2
  #  }
  
  #  ## Calculate gates
  #  p1 <- sf[ind1+1,]
  #  q1 <- sf[ind1-1,]
  #  m1 <- calcSlope(p1,q1)
  #  p2 <- a[ind2+1,]
  #  q2 <- a[ind2-1,]
  #  m2 <- calcSlope(p2,q2)
  
  #  ## Extract the indexes inside the gates
  #  data.points <- exprs(flow.frame)[,chans]
  #  lu.p1 <- .luline(data.points, p1, m1, !crt1)
  #  lu.p2 <- .luline(data.points, p2, m2, !crt2)x[which(x>=
  #  ind <- intersect(lu.p1, lu.p2)
  
  #  ## Generating the output CellPopulation
  #  snglt@index <- ind
  #  snglt@cell.count <- length(ind)
  #  i <- which(!is.na(exprs(flow.frame)[,chans[1]]))
  #  snglt@proportion <- snglt@cell.count/length(i) * 100
  #  snglt@flow.frame <- flow.frame
  #  exprs(snglt@flow.frame)[-ind,] <- NA
  #  if(length(ind)!=0){
  #    X <- exprs(snglt@flow.frame)[ind,chans]
  #    filter <- chull(X)
  #    filter <- X[c(filter,filter[1]),]
  #  }else
  #    filter <- as.matrix(NA)
  #  snglt@filter <- filter
  
  return(snglt)
}

.ellipseGating <- function(flow.frame, channels, position, gates, scale = 0.95, ellip.gate = FALSE, ...){
  
  ##====================================================================================================================================
  ## Fits an ellipse to the output of density-estimate gating 'flowDensity(.)', and extracts the inside events
  ##  Args:
  ##   flow.frame: a 'flowFrame' object
  ##   channels: a vector of two channel names or their equivalent integer number
  ##   position: a vector of two logical values specifying the position of the cell subset of interest on the 2D plot
  ##   gates: the 'gates' slot in the output of the 'flowDensity(.)' function, can be 'missing'
  ##   scale: an integer in [0,1) that scales the size of ellipse
  ##   ellip.gate: if TRUE, the ellipse is used for the purpose of gating, otherwise returns the rectangle gating results
  ##  Value:
  ##     events inside the ellipse/rectangle gate
  ## Author:
  ## M. Jafar Taghiyar
  ##-------------------------------------------------------------------------------------------------------------------------------------
  f.exprs <- exprs(flow.frame)
  if(missing(gates))
    gates <- flowDensity(flow.frame, channels, position, ...)@gates
  if(ellip.gate){
    new.f <- flow.frame
    ind <- which(is.na(exprs(flow.frame)[,channels[1]]))
    if(length(ind) != 0)
      exprs(new.f) <- exprs(flow.frame)[-ind, ]
    f.sub.ind <- .subFrame(new.f, channels, position, gates)
    eg <- dataEllipse(exprs(new.f)[f.sub.ind,channels[1]],
                      exprs(new.f)[f.sub.ind,channels[2]],
                      pch=".",
                      center.cex=0,
                      grid=F,
                      levels=c(0,scale),
                      add=F,
                      draw=F
    )
    
    ## Extract the coordinates of the elliptic gate
    p <- list()
    tmp <- names(eg)[2]
    p$x <- eg[[2]][,1]
    p$y <- eg[[2]][,2]
    
    ## Gate based on the ellipse not the rectangle produced by 'gates'
    in.p <- inpoly(x=exprs(new.f)[,channels[1]], y=exprs(new.f)[,channels[2]], POK=p)
    
    tmp.in.p <- vector(mode='numeric', length=length(exprs(flow.frame)[,channels[1]]))
    tmp.in.p[ind] <- NA
    if(length(ind) != 0)
      tmp.in.p[-ind] <- in.p
    else
      tmp.in.p <- in.p
    exprs(flow.frame)[which(tmp.in.p!=1), ] <- NA
  }else{
    f.sub.ind <- .subFrame(flow.frame, channels, position, gates)
    if(length(f.sub.ind)==0)
      exprs(flow.frame)[,] <- NA
    else
      exprs(flow.frame)[-f.sub.ind,] <- NA
  }
  return(flow.frame)
}

.luline <- function(data.points, point, m, up = TRUE, tp=NULL){
  ##================================================================================================================
  ## Returns the indexes of the points above and below a given line with slope m and the given 'point' coordination
  ## Args:
  ##   data.points: a matrix of point with 2 columns specifying the x and y coordinates
  ##   point: a point on the line with c(x,y) coordination
  ##   m: slope of the line
  ##   up: if TRUE the points above the line are returned, else below
  ##   tp: a test point with c(x,y) coordination to be evaluated if it is above or below the line
  ## Value:
  ##   the indexes of the points below and above the line
  ##-----------------------------------------------------------------------------------------------------------------
  yCalc <- function(x, m, d) m*x+d
  d <- point[2]-m*point[1]
  if(missing(data.points)){
    dp <- yCalc(tp[1],m,d)
    return(as.logical(ifelse(up, tp[2]-dp>0, tp[2]-dp<0)))
  }else{
    dp <- sapply(data.points[,1], FUN=function(x) yCalc(x=x,m=m,d=d))
    if(up)
      res <- which(data.points[,2]-dp>0)
    else
      res <- which(data.points[,2]-dp<0)
    return(res)
  }
}

