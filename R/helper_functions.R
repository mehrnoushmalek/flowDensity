.densityGating <- function(obj, channel, n.sd = 1.5, use.percentile = FALSE,  percentile = NA,use.upper=FALSE, upper = NA,verbose=TRUE,twin.factor=.98,
                           bimodal=F,after.peak=NA,alpha = 0.1, sd.threshold = FALSE, all.cuts = FALSE,
                           tinypeak.removal=tinypeak.removal, adjust.dens = 1,count.lim=20,magnitude = .3,slope.w=4,seq.w=4, spar=.4,...){

  ##========================================================================================================================================
  ## 1D density gating method
  ## Args:
  ##   obj: a 'FlowFrame' object, vector or a density object
  ##   channel: a channel's name or an integer to specify the channel
  ##   n.sd: an integer that is multiplied to the standard deviation to determine the place of threshold if 'sd.threshold' is 'TRUE'
  ##   use.percentile: if TRUE, forces to return the 'percentile'th threshold
  ##   use.upper: if True, forces to return the inflection point based on the first (last) peak if upper=F (upper=T)
  ##   percentile: a value in [0,1] that is used as the percentile. The default value is 0.95.
  ##   upper: if 'TRUE', it finds the change in the slope after the peak with index 'peak.ind'
  ##   verbose: When FALSE doesn't print any messages
  ##   twin.factor: reverse of tinypeak.removal for ignoring twinpeaks
  ##   bimodal: if TRUE, splits population closer to 50-50 when there are more than 2 peaks
  ##   after.peak: if TRUE (FALSE) and bimodal=FALSE, it chooses a cuttoff after(before) the maximum peak.
  ##   alpha: a value in [0,1) specifying the significance of change in the slope which would be detected
  ##   sd.threshold: if TRUE, it uses 'n.sd' times standard deviation for gating
  ##   all.cuts: if TRUE, it returns all the cutoff points whose length+1 can roughly estimate the number of cell subsets in that dimension
  ##   tinypeak.removal: a values in [0,1] for ignoring tiny peaks in density
  ##   adjust.dens: The smoothness of density, it is same as adjust in density(.).The default value is 1 and should not be changed unless necessary
  ##  seq.w: window for generating the sequence of points in density to use in .trackSlope
  ## spar, is the smoothing factor in spline function, default is .4
  ## Value:
  ##   cutoffs, i.e. thresholds on the 1D data
  ## Authors:
  ##   M. Jafar Taghiyar & Mehrnoush Malek
  ##-----------------------------------------------------------------------------------------------------------------------------------------
  if(class(obj)=="numeric"|class(obj)=="vector"){
    x<-obj
    n<- which(!is.na(x))
    channel <-NA
  dens <- .densityForPlot(data = x, adjust.dens,spar=spar,...)
  stdev <- sd(x,na.rm = T)
  }else if(class(obj)=="density"){
   warning("Percentiles and SD would be meaningless when dealing with density objects. Make sure to not use them with density objects.")
    dens <- obj
    channel <-NA
    stdev <- sd(dens$y,na.rm = T)
    x<-obj
   n<- dens$y
  }else{
   x <- exprs(obj)[, channel]
   n<- which(!is.na(x))
  dens <- .densityForPlot(data = x, adjust.dens,spar=spar,...)
  stdev <- sd(x,na.rm = T)
  }

  if (length(n)< count.lim)
  { 
    if(verbose)
       cat("Less than", count.lim, "cells, returning -Inf as a threshold.","\n")
    return(-Inf)
  }

  
  med <- median(dens$x,na.rm = T)
  if(is.numeric(channel))
    channel <- colnames(obj)[channel]
  cutoffs <- c()
  all.peaks <- .getPeaks(dens, peak.removal=tinypeak.removal)
  peak.ind <- all.peaks[[2]]
  peaks <- all.peaks[[1]]
  len <- length(peaks)
  if (len==0)
  {
    if(verbose)
      print("Failed to find a peak, check the density.")
    return(NA)
  }
  for(i in 1:(len-1))
    cutoffs[i] <- .getIntersect(dens, p1 = peaks[i],p2 =  peaks[i+1])
 
  if(use.percentile)
   return(quantile(x, probs = percentile, na.rm=T))
  if(use.upper){
    peak.ind <- ifelse(upper, peak.ind[len], peak.ind[1])
    track.slope.cutoff <- .trackSlope(dens, peak.ind=peak.ind, upper=upper, alpha=alpha,magnitude = magnitude,w=slope.w,seq.w=seq.w)
    if(is.infinite(track.slope.cutoff))
      return(ifelse(upper, max(dens$x), min(dens$x)))
    else
     return(track.slope.cutoff)
  }

 if(len==1){
    if(verbose)
      print("Only one peak is found, set standard deviation, percentile, or upper arguments accordingly")
    if(!is.na(percentile)){
      cutoffs <- quantile(x, probs = percentile, na.rm=T)
    } else if (!is.na(upper))
    {
      track.slope.cutoff <- .trackSlope(dens, peak.ind=peak.ind, upper=upper, alpha=alpha, magnitude = magnitude,w=slope.w,seq.w=seq.w)
      if(is.infinite(track.slope.cutoff))
        cutoffs <- ifelse(upper, max(dens$x), min(dens$x))
      else
        cutoffs <- track.slope.cutoff
    }else if (sd.threshold)
    {
      upper <- as.logical(ifelse(peaks>med, FALSE, TRUE))
      cutoffs <- ifelse(upper, peaks+n.sd*stdev, peaks-n.sd*stdev)
    }else{
      flex.point <- .getFlex(dens, peak.ind=peak.ind,magnitude = magnitude)
      if(!is.na(flex.point) & !is.null(flex.point))
      { 
        if(verbose)
          print("Cutoff is based on inflection point.")
        cutoffs <- flex.point
      }else{
        if(is.na(upper))
          upper <- as.logical(ifelse(peaks>med, FALSE, TRUE))
        track.slope.cutoff <- .trackSlope(dens, peak.ind=peak.ind, upper=upper, alpha=alpha, magnitude = magnitude,w=slope.w,seq.w=seq.w)

        stdev.cutoff <- ifelse(upper, peaks+n.sd*stdev, peaks-n.sd*stdev)
        cutoffs <- ifelse(is.infinite(track.slope.cutoff)|sd.threshold, stdev.cutoff, track.slope.cutoff) # use the 's.threshold' if the gate from 'trackSlope()' is too loose
        if(verbose)
          print(" cuttoff is based on tracking the slope and cpmparing the position of the peak and mean of the population.")
      }
    }
   
   return(cutoffs)
 }else{
   if(all.cuts)
     return(cutoffs)

     cutoffs <- .getScoreIndex(dat=x, peaks = all.peaks,cutoffs = cutoffs,percentile = percentile,adjust.dens = adjust.dens,seq.w=seq.w,w=slope.w,sd.threshold,n.sd=n.sd,
                               upper=upper,alpha=alpha,twin.factor=twin.factor,bimodal=bimodal,after=after.peak,magnitude=magnitude, ...)

   }
   return(cutoffs)
 }





.deGate2D <- function(f, channels, position, n.sd=c(1.5,1.5), use.percentile=c(F,F), percentile=c(NA,NA), use.upper=c(F,F),
                      upper=c(NA,NA), verbose=c(TRUE,TRUE), twin.factor=c(.98,.98), bimodal=c(F,F),filter=NA,tinypeak.removal=c(1/25,1/25),
                      after.peak=c(NA,NA),alpha=c(0.1,0.1), sd.threshold=c(F,F), use.control=c(F,F), control=c(NA,NA),
                      gates=c(NA,NA),all.cuts=c(F,F),remove.margins=F,count.lim=3, spar=spar,seq.w=4,w=4,...){

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
 
  i <- which(!is.na(exprs(f)[,channels[1]]))
  if(length(i)==0)
    warning('invalid flowFrame input: This flowFrame has 0 cells, creating dummy gate.')
  if(length(i)<count.lim)
  {
    if(is.numeric(channels))
      channels <- c(colnames(f)[channels[1]], colnames(f)[channels[2]])
    filter <- cbind(c(-Inf, -Inf,-Inf,-Inf), c(-Inf, -Inf,-Inf,-Inf))
    colnames(filter) <- channels
    cell.population <- new("CellPopulation",
                           flow.frame=f,
                           proportion=100,
                           cell.count=length(i),
                           channels=channels,
                           position=position,
                           gates=c(0,0),
                           filter=filter,
                           index=i)
    return(cell.population)
  }
  args <- names(list(...))
  eligible.args <- c('use.percentile', 'use.upper','upper', 'percentile', 'sd.threshold','filter', 'count.lim','gates',
                     'n.sd','twin.factor','bimodal','after.peak','use.control', 'control', 'alpha', 'scale', 'ellip.gate','tinypeak.removal','verbose')
  if(length(setdiff(args, eligible.args)!=0))
    warning('unused argument(s): ', setdiff(args, eligible.args))
  col.nm <- c(grep(colnames(f), pattern="FSC"), grep(colnames(f), pattern="SSC"))
  if(remove.margins)
  {
    if(length(col.nm)==0)
      warning('No forward/side scatter channels found in control data for first channel, margin events not removed.')
    else
      f <- nmRemove(f, col.nm)
  }
  if(is.na(gates[1])&is.na(gates[2])&is.na(position[1])&is.na(position[2]))
    stop("Improper 'position' value, one position must be not 'NA' or 'gates' should be provided" )
  if (is.matrix(filter) | is.data.frame(filter) )
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
  if(!is.na(filter))
    warning("Filter provided, is not a matrix or data.frame. You can fix it otherwise, flowDensity will estimates the threshold.")
  if(is.na(gates[1])){
    if(use.control[1]&!is.na(position[1])){
      if(is.na(control[1]))
        stop("Missing 'control' data for first channel, set 'use.control' to FALSE if control should not be used")
      
      f.control1 <- NA
      control1.class <- class(control[1][[1]])
      if(control1.class=="flowFrame"){
        f.control1 <- control[1][[1]]
      }else if(control1.class=="CellPopulation"){
        f.control1 <- control[1][[1]]@flow.frame
      }else{
        stop("The 'control' input for the first channel must be a flowFrame or CellPopulation object")
      }
      col.nm.control1 <- c(grep(colnames(f.control1), pattern="FSC"), grep(colnames(f.control1), pattern="SSC"))
      if(remove.margins)
      {
        if(length(col.nm.control1)==0)
          warning('No forward/side scatter channels found in control data for first channel, margin events not removed.')
        else
          f.control1 <- nmRemove(f.control1, col.nm.control1)
      }
      gates[1] <- .densityGating(obj=f.control1, channel=channels[1], n.sd=n.sd[1], use.percentile=use.percentile[1], percentile=percentile[1], use.upper=use.upper[1], upper=upper[1],
                                 verbose=verbose[1],twin.factor=twin.factor[1],bimodal=bimodal[1],after.peak = after.peak[1],alpha=alpha[1],
                                 sd.threshold=sd.threshold[1],all.cuts=all.cuts[1],tinypeak.removal=tinypeak.removal[1],count.lim=count.lim,
                                 spar=spar,w=w,seq.w=seq.w)
    }
    else if(!is.na(position[1]))
      gates[1] <- .densityGating(obj=f, channel=channels[1], n.sd=n.sd[1], use.percentile=use.percentile[1], percentile=percentile[1], use.upper=use.upper[1], upper=upper[1],
                                 verbose=verbose[1],twin.factor=twin.factor[1],bimodal=bimodal[1],after.peak = after.peak[1],alpha=alpha[1],
                                 sd.threshold=sd.threshold[1],all.cuts=all.cuts[1],tinypeak.removal=tinypeak.removal[1],count.lim=count.lim,
                                 spar=spar,w=w,seq.w=seq.w)
    else
      gates[1] <- min(exprs(f)[,channels[1]],na.rm = T)
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
      if(remove.margins)
      {
        if(length(col.nm.control2)==0)
          warning('No forward/side scatter channels found in control data for first channel, margin events not removed.')
        else

        f.control2 <- nmRemove(f.control2, col.nm.control2)
      }
      gates[2] <- .densityGating(obj=f.control2, channel=channels[2], n.sd=n.sd[2], use.percentile=use.percentile[2], percentile=percentile[2], use.upper=use.upper[2], upper=upper[2],
                                 verbose=verbose[2],twin.factor=twin.factor[2],bimodal=bimodal[2],after.peak = after.peak[2],alpha=alpha[2],
                                 sd.threshold=sd.threshold[2],all.cuts=all.cuts[2],tinypeak.removal=tinypeak.removal[2],count.lim=count.lim,
                                 spar=spar,w=w,seq.w=seq.w)
    }
    else if(!is.na(position[2]))
      gates[2] <- .densityGating(obj=f, channel=channels[2], n.sd=n.sd[2], use.percentile=use.percentile[2], percentile=percentile[2], use.upper=use.upper[2], upper=upper[2],
                                 verbose=verbose[2],twin.factor=twin.factor[2],bimodal=bimodal[2],after.peak = after.peak[2],alpha=alpha[2],
                                 sd.threshold=sd.threshold[2],all.cuts=all.cuts[2],tinypeak.removal=tinypeak.removal[2],count.lim=count.lim,
                                 spar=spar,w=w,seq.w=seq.w)
    else
      gates[2] <- min(exprs(f)[,channels[2]],na.rm = T)
  }
  }
  new.f <- .ellipseGating(f, channels=channels, position=position, gates=gates, ...)
  cell.index <- which(!is.na(exprs(new.f)[,channels[1]]))
  cell.count <- length(cell.index)
  proportion <- cell.count/length(i) * 100
  
  if(cell.count>1){
    X <- exprs(new.f)[cell.index,channels]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
  } else{
     if(is.numeric(channels))
       channels <- c(colnames(f)[channels[1]], colnames(f)[channels[2]])
     filter <- cbind(c(-Inf, -Inf,-Inf,-Inf), c(-Inf, -Inf,-Inf,-Inf))
     colnames(filter) <- channels

   }

  
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
                         index=cell.index )

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
  dens.chan1 <- .densityForPlot(data.chan1, adjust.dens, ...)
  graphics::plot(dens.chan1, xlim=range(data.chan1), type="l", frame.plot=F, axes=F, ann=F)
  
  ## Fill the pdf of channels[1]
  pts <- .densRange(dens.chan1$x, dens.chan1$y, gates[1], pos[1])
  polygon(pts$x, pts$y, col="#00ff0044", border= "#00ff00ff")
  
  ## Add the bottom line and the density-estimated gate
  abline(h=0)
  if(!is.na(pos[1]))
    abline(v=gates[1], lty="dashed", col=colors()[290], lwd=1.5)
  
  ## Plot the pdf of channels[2]
  data.chan2 <- exprs(f)[,channels[2]]
  dens.chan2 <- .densityForPlot(data.chan2, adjust.dens, ...)
  graphics::plot(dens.chan2$y, dens.chan2$x, ylim=range(data.chan2), xlim=rev(range(dens.chan2$y)),
       type="l", col=1, frame.plot=F, axes=F, ann=F)
  
  ## Fill the pdf of channels[2]
  pts <- .densRange(dens.chan2$x, dens.chan2$y, gates[2], pos[2])
  polygon(pts$y, pts$x, col="#ff450055", border="#ff4500ff")
  
  ## Add the bottom line and the density-estimated gate
  abline(v=0)
  if(!is.na(pos[2]))
    abline(h=gates[2], lty="dashed", col=colors()[290], lwd=1.5)
  
  ## Plot the scatter plot of the flowframe
    plotDens(obj=f, channels=channels, ...)
  
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

.getPeaks <- function(d,peak.removal=1/25,twin.factor=1){
  
  ##=====================================================================================================================
  ## Finds the peaks in the given density
  ## Args:
  ##   dens: density of the channel whose peaks are going to be found. It is a variable of class 'density'.
  ##   w: the length of the window where the function searches for the peaks. If all peaks required, use the default w=1.
  ## Value:
  ##   peaks in the density of the provided channel
  ##---------------------------------------------------------------------------------------------------------------------
  if (all(is.na(d)))
     return(NA)
  d.y <- d$y
  w <- 1
  peaks <- c()
  peaks.ind <- c()
  for(i in 1:(length(d.y)-2*w)){
    if(d.y[i+w] > d.y[(i+w+1):(i+2*w)] && d.y[i+w] > d.y[i:(i+w-1)] && d.y[i+w] > peak.removal*max(d.y)){ # also removes tiny artificial peaks less than ~%4 of the max peak
      peaks <- c(peaks, d$x[i+w])
      peaks.ind <- c(peaks.ind, i+w)
      
    }
  }
  peaks.height <- d.y[peaks.ind]
if (all(is.na(peaks)))
{
warning("No peaks could be found, returning the maximum value of density.")
  peaks <-d$x[which.max(d$y)]
  peaks.height<-d$y[which.max(d$y)]
  peaks.ind <-which.max(d$y)
}
max.peak <- which.max(peaks.height)
if ( length(peaks)>1 & twin.factor<1)
{
  twins <-which(peaks.height>peaks.height[max.peak]*twin.factor)
 twins<-  twins[-which(twins==max.peak)]
  if (length(twins)>0)
{
   if( .getIntersect (d,peaks[tail(twins,1)],peaks[max.peak])>peaks.height[max.peak]*twin.factor*.98)
  {
   peaks<-peaks[-twins]
    peaks.height<-peaks.height[-twins]
    peaks.ind <- peaks.ind[-twins]
}
}
}
  return(list(Peaks=peaks, P.ind=peaks.ind,P.h=peaks.height))
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

return.bimodal<-function(x,cutoffs)
{
  scores<-c()
  for (i in 1:length(cutoffs))
  {
    p.1<-length(which(x>cutoffs[i]))/length(x)
    p.2<-length(which(x<cutoffs[i]))/length(x)
    scores <- c(scores,abs(p.1-p.2))
  }
  return(cutoffs[which.min(scores)])
}

.getScoreIndex <- function(dat, peaks,cutoffs,percentile,adjust.dens,upper,alpha,twin.factor,magnitude = magnitude,bimodal,spar=spar
                           ,after,w=w,seq.w=seq.w, sd.threshold,n.sd,...){
  ##==========================================================================================
  ## Returns the cutoff based upon which the .densityGating() function decides on the thresholds
  ## Args:
  ##   dens: density of the channel to be investigated
  ##   cutoffs: the min intersections returned by getIntersect() function
  ##   twin.factor: A value from 0-1 for ignoring twin peaks
  ##   percentile.cut: returns percentile when 1 peak is detected after twin-peaks removal
  ##   upper: returns upper when 1 peak is detected after twin-peaks removal
  ##   bimodal: chooses a cutoff which splits population closer to 50-50
  ##   after: default (NA), when is TRUE, the cutoff closes to the right of the max-peak will be chosen
  ## Value:
  ##   the selected cutoff based on defined arguments
  ##------------------------------------------------------------------------------------------

  if (class(dat)!="density"){
     dens <- .densityForPlot(dat, adjust.dens=1,spar=spar, ...)
  }else{
   warning("Percentiles and SD would be meaningless when dealing with density objects. Make sure to not use them with density objects.")
   dens <- dat
  dat <- dat$y
}
  mins <-c()
  h <- c()
  fun <- splinefun(dens)
  h.c <- fun(cutoffs);
  h.p <- fun(peaks$Peaks);
  d<- .getDist(peaks$Peaks)
  max.peak<-max(peaks$P.h)

  for(i in 1:(length(peaks$Peaks)-1)){
    valley <- cutoffs[i]

    if ((dens$y[which.min(abs(dens$x-valley))]>twin.factor*dens$y[peaks$P.ind[i]])& (dens$y[which.min(abs(dens$x-valley))]>twin.factor*dens$y[peaks$P.ind[i+1]]))
    {

      valley <- NA

    }else{
      mins<-c(mins,valley)
    }

  }
    if (suppressWarnings(all(is.na(mins))))
    {
      if(!is.na(percentile)){
        return(quantile(dat, percentile, na.rm=T))
      } else if (!is.na(upper))
      {
        return(.trackSlope(dens=dens, peak.ind=ifelse(upper,yes = tail(peaks$P.ind,1),no =peaks$P.ind[1]), 
                           upper=upper, alpha=alpha, magnitude = magnitude, w=w,seq.w=seq.w))

      }else if (!is.na(sd.threshold))
    {
      upper <- as.logical(ifelse(tail(peaks$Peaks,1)>mean(dens$x,na.rm=T), FALSE, TRUE))
     return(ifelse(upper, yes=tail(peaks$Peaks,1)+n.sd*sd(dens$x,na.rm=T), no=tail(peaks$Peaks,1)-n.sd*sd(dens$x,na.rm=T)))
    }else{
    stop("After using twin.factor, no threshold can be claculated, please set upper, percentile or sd.threshold")
}
      
    }else if(length(mins)>1)
    {

    if (bimodal)
    {
      return(return.bimodal(dat,mins))
    }else if(!is.na(after))
    {
      cutoff <-ifelse(after,yes = mins[mins>peaks$Peaks[which.max(peaks$P.h)]][1],
         no=tail(mins[which(mins<peaks$Peaks[which.max(peaks$P.h)])],1))
      if(is.na(cutoff))
        cutoff <-ifelse(after,yes =tail(mins,1),no=mins[1])
      return(cutoff)
    }else{
      ##Old method of flowDensity
      if(length(which(h.c==0))==0)
        return(cutoffs[which.max(d/h.c)])
      else
        return(cutoffs[which.max(d[which(h.c==0)])])
    }
      
    }else{
      return(mins)
    }
 
  
}

.trackSlope <- function(dens, peak.ind, alpha, upper=T, w,seq.w=4, return.slope=F,start.lo=1,rev=F,magnitude=.3){
  
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
  start.p <- ifelse(rev,yes = peak.ind-w,no = start.lo+w)
  end.p <- ifelse(rev,yes = start.lo+w,no = peak.ind-w)

  lo <- ifelse(upper, yes = peak.ind+w, no = start.p)
  up <- ifelse(upper, yes = length(d.y)-w, no = end.p)
  w <-ifelse(rev,yes = -w,no = w) 
  if(lo>up &w>0)
  {
    w<-2
    lo <- ifelse(upper, yes = peak.ind+w, no = start.p)
    up <- ifelse(upper, yes = length(d.y)-w, no = end.p)
  }
  if (mode(tryCatch(seq(from=lo,to=up,by=w), error = function(x) return("error")))=="character")
  {
    warning("not enough point to calculate upper, returning NA.")
    return(NA)
  }
  for(i in seq(from=lo,to=up,by=seq.w)){

    slope <- c(slope, abs((d.y[i+w]-d.y[i-w])/(d.x[i+w]-d.x[i-w])))
    heights <-c(heights,d.y[i])
  }
  ## try to estimate the optimum alpha using the distribution of slopes
  #   if(missing(alpha))
  #     alpha <- ifelse(sd(slope) < max(slope)/10, 0.01, 0.1)
  if(return.slope)
    return(slope)
  small.slopes <- which(slope < (max(slope)*alpha))
  small.slopes <- small.slopes[which(heights[small.slopes] < magnitude*d.y[peak.ind])]
  len <- length(small.slopes)
  if(len==0)
    return(ifelse(upper, Inf, -Inf))
  ind <- ifelse(upper, yes=small.slopes[1], no=small.slopes[len]) #if upper, the first large change of slope is given else the last change of slope is returned
  
  return(ifelse(upper,yes= d.x[peak.ind+w*(ind)], no=d.x[(lo+(ind-1)*w)]))
}

.getFlex <- function(dens, peak.ind, w=2,magnitude = magnitude){

  
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
  upper.slope <- .trackSlope(dens=dens, peak.ind=peak.ind, alpha=0.1, upper=T, w=w,seq.w = seq.w, return.slope=T,magnitude = magnitude)
  lower.slope <- .trackSlope(dens=dens, peak.ind=peak.ind, alpha=0.1, upper=F, w=w,seq.w = seq.w, return.slope=T,magnitude = magnitude)
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
{    
    warning("Inflection points cannot be calculated, returning the first peak + sd.")
    return(dens$x[which.max(dens$y)]+sd(dens$x,na.rm=T))
}
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

.densityForPlot <- function(data, adjust.dens=1,spar=.4,...){
  if (length(data)<2)
  {
      warning ("Less than 2 cells, returning NA.")
      return(NA)
 }
  dens <- density(data[which(!is.na(data))], adjust=adjust.dens)
  dens <- smooth.spline(dens$x, dens$y, spar=spar, ...)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

.subFrame <- function(f, channels, position, gates, filter,include.equal=T){
  
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
    #in.p <-  inpoly(x=f.exprs[,channels[1]], y=f.exprs[,channels[2]], POK=POK)
     in.p <- polyclip::pointinpolygon(list(x=f.exprs[,channels[1]], y=f.exprs[,channels[2]]),list(x=filter[,1],y=filter[,2]))
      tmp.in.p <- vector(mode='numeric', length=length(exprs(f)[,channels[1]]))
    if(length(ind) != 0){
      tmp.in.p[ind] <- NA
      tmp.in.p[-ind] <- in.p
    }else{
      tmp.in.p <- in.p
    }
     f.sub.ind <- which(tmp.in.p!=0)
    return(f.sub.ind)
  }else if(is.na(position[1])&is.na(position[2])){
    stop("Improper 'position' value, one position must be not 'NA'")
}else{
  if(is.na(position[1])){
    if(position[2])
    {
      f.sub.ind <- which(f.exprs[,channels[2]] > gates[2])
    }else{
      f.sub.ind <- which(f.exprs[,channels[2]] < gates[2])
    }
    if(include.equal)
      f.sub.ind <-c(f.sub.ind,which(f.exprs[,channels[2]]==gates[2]))
  }else if(is.na(position[2])){
    if(position[1]){
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1])
    }else{
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1])
    }
    if(include.equal)
      f.sub.ind <-c(f.sub.ind,which(f.exprs[,channels[1]]==gates[1]))
  }else{
    if(position[1]&position[2]){
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1] & f.exprs[,channels[2]] > gates[2])
    }else if (position[1]&!position[2]){
      f.sub.ind <- which(f.exprs[,channels[1]] > gates[1] & f.exprs[,channels[2]] < gates[2])
    }else if (!position[1]&position[2]){
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1] & f.exprs[,channels[2]] > gates[2])
    }else{
      f.sub.ind <- which(f.exprs[,channels[1]] < gates[1] & f.exprs[,channels[2]] < gates[2])
    }
    if(include.equal)
     {
       if(!is.na(gates[2]) & is.na(gates[1]))
         f.sub.ind <-c(f.sub.ind, which(f.exprs[,channels[2]]==gates[2]))
        if(!is.na(gates[1]) & is.na(gates[2]))
          
         f.sub.ind <-c(f.sub.ind, which(f.exprs[,channels[1]]==gates[1]))
        if(!is.na(gates[1]) & !is.na(gates[2]))
          f.sub.ind <-c(f.sub.ind, which(f.exprs[,channels[2]]==gates[2] &f.exprs[,channels[1]]==gates[1] ))
 
    }
  }
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
  
  if(!missing(filter)){
    ind <- .subFrame(f=flow.frame, channels=channels, filter=filter)
    sub.filter <- filter
  }else{
    ind <- .subFrame(f=flow.frame, channels=channels, position=position, gates=gates,include.equal=T)
    X1 <-exprs(flow.frame)[ind,channels]
    filt <- chull(X1) 
    sub.filter <- X1[c(filt,filt[1]),]
  }
  
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
    poly1 <- list(list(x=filter[,1],y=filter[,2])  )

    poly2 <- list(list(x=sub.filter[,1],y=sub.filter[,2])  )
     intersection <- polyclip::polyclip(poly1,poly2,op="intersection")
    if(length(intersection[[1]][1])>0){
       temp <-tryCatch(polyclip::polyclip(poly1,poly2,op="minus"), error=function(ex) return(1))
       if (mode(temp)=="numeric")
       {
         temp <-polyclip::polyclip(poly1,polyclip:polylineoffset(poly2,1)[[1]])
       }
       filter<-cbind(c(temp[[1]]$x,temp[[1]]$x[1]),c(temp[[1]]$y,temp[[1]]$y[1]))
    }
  }else

  filter <- cbind(c(-Inf, -Inf,-Inf,-Inf), c(-Inf, -Inf,-Inf,-Inf))
  colnames(filter) <- channels
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



.ellipseGating <- function(fcs, channels, position, gates, scale = 0.95, ellip.gate = FALSE, ...){
  
  ##====================================================================================================================================
  ## Fits an ellipse to the output of density-estimate gating 'flowDensity(.)', and extracts the inside events
  ##  Args:
  ##   fcs: a 'flowFrame' object
  ##   channels: a vector of two channel names or their equivalent integer number
  ##   position: a vector of two logical values specifying the position of the cell subset of interest on the 2D plot
  ##   gates: the 'gates' slot in the output of the 'flowDensity(.)' function, can be 'missing'
  ##   scale: an integer in [0,1) that scales the size of ellipse
  ##   ellip.gate: if TRUE, the ellipse is used for the purpose of gating, otherwise returns the rectangle gating results
  ##  Value:
  ##     events inside the ellipse/rectangle gate
  ## Author:
  ## M. Malek, J. Taghiyar
  ##-------------------------------------------------------------------------------------------------------------------------------------
  if(missing(gates))
    gates <- flowDensity(fcs, channels, position, ...)@gates
  if(ellip.gate){
    new.f <- fcs
    ind <- which(is.na(exprs(fcs)[,channels[1]]))
    if(length(ind) != 0)
      exprs(new.f) <- exprs(fcs)[-ind, ]
    f.sub.ind <- .subFrame(new.f, channels, position, gates)
    if (length(f.sub.ind)<5)
    {
        return(fcs[f.sub.ind,])
    }else{
      eg <- dataEllipse(exprs(new.f)[f.sub.ind,channels[1]],
                      exprs(new.f)[f.sub.ind,channels[2]],
                      pch=".",
                      center.cex=0,
                      grid=F,
                      levels=c(0,scale),
                      add=F,
                      draw=F)
    
    ## Extract the coordinates of the elliptic gate
   # p <- list()
    #tmp <- names(eg)[2]
    #p$x <- eg[[2]][,1]
    #p$y <- eg[[2]][,2]
    
    ## Gate based on the ellipse not the rectangle produced by 'gates'
    #in.p <- inpoly(x=exprs(new.f)[,channels[1]], y=exprs(new.f)[,channels[2]], POK=p)
       in.p <- polyclip::pointinpolygon(list(x=f.exprs[,channels[1]], y=f.exprs[,channels[2]]),
                                        list(x=eg[[2]][,1],y=eg[[2]][,2]))
       tmp.in.p <- vector(mode='numeric', length=length(exprs(fcs)[,channels[1]]))
      if(length(ind) != 0){
        tmp.in.p[ind] <- NA
        tmp.in.p[-ind] <- in.p
      }else{
        tmp.in.p <- in.p
      }
      exprs(fcs)[which(tmp.in.p==0), ] <- NA
   }
   }else{
    f.sub.ind <- .subFrame(fcs, channels, position, gates)
    if(length(f.sub.ind)==0)
      exprs(fcs)[,] <- NA
    else
      exprs(fcs)[-f.sub.ind,] <- NA
   }
   return(fcs)

}


