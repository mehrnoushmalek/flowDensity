\name{deGate}
\docType{methods}
\alias{deGate}
\title{
  1D density gating method
}
\description{
    Find the best threshold for a single channel in flow cytometry
    data based on its density distribution.
}
\usage{
deGate(obj,channel, n.sd = 1.5, use.percentile = FALSE,  percentile =NA,use.upper=FALSE, upper = NA,verbose=TRUE,twin.factor=.98,
                   bimodal=F,after.peak=NA,alpha = 0.1, sd.threshold = FALSE, all.cuts = FALSE,
                   tinypeak.removal=1/25, adjust.dens = 1,count.lim=20,magnitude=.3,slope.w=4,seq.w = 4, spar = 0.4, ...)
}
\arguments{
  \item{obj}{
  	 obj: a 'FlowFrame' object, 'CellPopulation' or 'GatingHierarchy'
  }
  \item{channel}{
    a channel's name or its corresponding index in the 'flow.frame'.
  }
  \item{n.sd}{
    an integer coefficient for the standard deviation
    to determine the threshold based on the standard deviation if 'sd.threshold' is TRUE.
  }
  \item{use.percentile}{
    if TRUE, forces to return the 'percentile'th threshold.
  }
  \item{percentile}{

     A value in [0,1] that is used as the percentile. The default is NA. If set to a value(n) and use.percentile=F, it returns the n-th percentile, for 1-peak populations. 
  }
    \item{use.upper}{
    Logical. If TRUE, forces to return the inflection point based on the first (last) peak if upper=F (upper=T). Default value is set to 'FALSE'
  }
  
   \item{upper}{
    if TRUE, finds the change in the slope at the tail of the density
    curve, if FALSE, finds it at the head. Default value is set to 'NA'.
  }
  
  \item{verbose}{
    Logical. If TRUE, Prints a message if only one peak is found, or when inflection point is used to set the gates.
  }
   \item{twin.factor}{
    a value in [0,1] that is used to exclude twinpeaks
  }
  
     \item{bimodal}{
     Logical. If TRUE, it returns a cutoff that splits population closer to 50-50, when there are more than two peaks.
  }
      \item{after.peak}{
     Logical. If TRUE, it returns a cutoff that is after the maximum peaks, when there are more than two peaks.
  }
  
  \item{alpha}{
    a value in [0,1) specifying the significance of change in the slope
    being detected. This is by default 0.1, and typically need not be changed.
  }
  \item{sd.threshold}{
    if TRUE, uses 'n.sd' times standard deviation as the threshold. Default value is set to 'FALSE'.
  }
 
  \item{all.cuts}{
    if TRUE, returns all the identified cutoff points, i.e. potential thresholds
    for that channel. Default value is set to 'FALSE'.
}
  \item{tinypeak.removal}{
    A number in [0,1] to exclude/include tiny peaks in density distribution.
  }
  
  \item{adjust.dens}{
    The smoothness of density in [0,Inf] to be used in density(.). The default value is 1 and should not be changed unless necessary

  }
  
  \item{count.lim}{
   minimum limit for events count in order to calculate the threshold. Default is 20, returning NA as threshold.
  }

 \item{magnitude}{
   A value between 0 and 1, for tracking a slope and reporting changes that are smaller than magnitude*peak_height
  }

\item{slope.w}{
 window.width for tracking slope. Default is 4, calculating a slope based on 4 points before and after the current point.
}

\item{seq.w}{
 value used for making the sequence of density points, used in trackSlope.
}

\item{spar}{
 value used in smooth.spline function, used in generating the density, default is 0.4.
}

\item{...}{
   Extra arguments to be passed to smoothSpline function.
  }

}

\value{
  an integer value (vector) of cutoff(s), i.e. threshold(s), on the specified channel
}

\details{

deGate works for GatingHierarchy, flowFrame, CellPopulation object or a numeric vector of data. In case the input is a numeric vector, channel doesn't need to provided, but the rest of arguments can be used to tune the outcome.
}

\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
#Find the threshold for CD20
cd19.gate <- deGate(f,channel="PerCP-Cy5-5-A")
# Gate out the CD20- populations using the notSubFrame
plotDens(f,c("APC-H7-A","PerCP-Cy5-5-A"))
abline(h=cd19.gate,lty=3,col=2)
}


\author{
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
}

\seealso{
\code{\link{getflowFrame}} \code{\link{notSubFrame}} \code{\link{flowDensity}}
}
