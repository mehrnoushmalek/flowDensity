\name{getPeaks}
\docType{methods}
\alias{getPeaks}
\title{
  Finding Peaks 
}
\description{
    Find all peaks in density along with their indices
}
\usage{
getPeaks(frame, channel,tinypeak.removal=1/25,...)
}
\arguments{
  \item{frame}{
  	a 'FlowFrame' object.
  }
  \item{channel}{
    a channel's name or its corresponding index in the 'flow.frame'.
  }
  
  \item{tinypeak.removal}{
    A number in [0,1] to exclude/include tiny peaks in density distribution.
  }

 \item{...}{
   Arguments passed to smoothSpline function, spar is 0.4.
  }
  
 
}
\value{
  a list, including peaks and their corresponding indices
}
\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
#Find the threshold for CD20
peaks <- getPeaks(f,channel="PerCP-Cy5-5-A",tinypeak.removal=1/30)
peaks
}


\author{
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
}

\seealso{
\code{\link{deGate}} \code{\link{notSubFrame}} \code{\link{flowDensity}}
}

