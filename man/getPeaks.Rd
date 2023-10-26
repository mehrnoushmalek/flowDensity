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
getPeaks(obj, channel,tinypeak.removal=1/25, adjust.dens=1,verbose=F,twin.factor=1,spar = 0.4,...)
}
\arguments{
  \item{obj}{
  	a 'FlowFrame', 'GatingHierarchy', 'CellPopulation' a density object or a numeric vector of  density.
  }
  \item{channel}{
    a channel's name or its corresponding index. If the input is numeric vector, channel is NA.
  }
  
  \item{tinypeak.removal}{
    A number in [0,1] to exclude/include tiny peaks in density distribution. Default is 1/25.
  }
 \item{adjust.dens}{
    The smoothness of density in [0,Inf] to be used in density(.). The default value is 1 and should not be changed unless necessary
  }

\item{verbose}{
 If TRUE, printing warnings.
}
\item{twin.factor}{
 If smaller than 1, peaks that are of greater than hieght as the maximum peak*twin.factor will be removed.
}
 \item{...}{
   Arguments passed to smoothSpline function, spar is 0.4.
  }
}
\value{
  a list, including peaks, their corresponding indices and height.
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

