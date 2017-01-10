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
deGate(flow.frame, channel, n.sd=1.5, use.percentile = FALSE, 
		percentile = 0.95,use.upper=FALSE, upper = NA, talk=TRUE, alpha = 0.1, 
		sd.threshold = FALSE, graphs = FALSE, all.cuts = FALSE, tinypeak.removal=1/25, adjust.dens=1)
}
\arguments{
  \item{flow.frame}{
  	a 'FlowFrame' object.
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
    a value in [0,1] that is used as the percentile. The default value is 0.95.
  }
   \item{upper}{
    if TRUE, finds the change in the slope at the tail of the density
    curve, if FALSE, finds it at the head. Default value is set to 'NA'.
  }
  
  \item{use.upper}{
    Logical. If TRUE, forces to return the inflection point based on the first (last) peak if upper=F (upper=T). Default value is set to 'FALSE'
  }
  
  \item{talk}{
    Logical. If TRUE, Prints a message if only one peak is found, or when inflection point is used to set the gates.
  }
  
  \item{alpha}{
    a value in [0,1) specifying the significance of change in the slope
    being detected. This is by default 0.1, and typically need not be changed.
  }
  \item{sd.threshold}{
    if TRUE, uses 'n.sd' times standard deviation as the threshold. Default value is set to 'FALSE'.
  }
  \item{graphs}{
    if TRUE, generates density distribution plot plus its corresponding threshold.
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
}
\value{
  an integer value (vector) of cutoff(s), i.e. threshold(s), on the specified channel
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
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
}

