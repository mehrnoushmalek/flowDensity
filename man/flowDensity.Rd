\name{flowDensity}
\docType{methods}
\alias{flowDensity}
\title{ Methods for Function \code{flowDensity} in Package \pkg{flowDensity} }
\usage{
	flowDensity(obj, channels, position, ...)
}
\description{
  flowDensity is an automated clustering algorithm which aims to emulate the current practice of manual sequential gating.
  It is designed to identify the predefined cell subsets based on the density distribution of the parent cell population by analyzing the peaks of the density curve.
  
  When the density distribution has only one peak:
  \itemize{
  \item the first argument that would be used is the percentile, the default value is 95th.
  \item If set to 'NA' then upper will be used when set to 'FALSE/TRUE'.
  \item If both upper and percentiles are set to 'NA' sd.threshold is used when set to 'TRUE'.
  \item If either all these arguments are 'NA', or one of them fails, the algorithm find the best cutoff based on the inflection points and by comparing the position of the peak and the mean of the population.
}
}

\arguments{
  \item{obj}{
  	a 'CellPopulation', 'flowFrame' or 'GatingHierarchy' object.
  }
  \item{channels}{
    a vector of two channel names or their corresponding indices.
  }
  \item{position}{
    a vector of two logical values specifying the position of the cell subset of interest on the 2D plot.
  }
 
  \item{...}{
  	This can be used to pass one of the following arguments:
  	\itemize{
  		\item 'use.percentile' if TRUE, returns the 'percentile'th threshold.

  		\item 'percentile' a value in [0,1] that is used as the percentile if 'use.percentile' is TRUE.

  		\item 'upper' if 'TRUE', it finds the change in the slope after the peak with index 'peak.ind'.
  		
  		\item 'use.upper' if 'TRUE', forces to return the inflection point based on the first (last) peak if upper=F (upper=T)
  		
  		\item 'twin.factor' a value in [0,1] that is used to exclude twinpeaks.
  		
  		\item 'bimodal' If TRUE, it returns a cutoff that splits population closer to 50-50, when there are more than two peaks.
  		
  		\item 'after.peak' If TRUE, it returns a cutoff that is after the maximum peaks, when there are more than two peaks.

  		\item 'sd.threshold' if TRUE, it uses 'n.sd' times standard deviation for gating.

  		\item 'n.sd' an integer that is multiplied to the standard deviation to determine the place 
  		of threshold if 'sd.threshold' is 'TRUE'.
  		
      \item 'tinypeak.removal' a vector of length 2, for sensitivity of peak finding for each channel. See deGate() for more information.
      
      \item 'filter' If provided it uses the given filter to gate the population.
      
  		\item 'use.control' if TRUE, it finds the threshold using a matched control population and uses it for gating.

  		\item 'control' a 'flowFrame' or 'CellPopulation' object used for calculating the gating threshold when 
      'use.control' is set to TRUE. If a control population is used, the other arguments ('upper', 'percentile', 
      etc.) are applied to the control data when finding the threshold (i.e. not to 'obj').

  		\item 'alpha' a value in [0,1) specifying the significance of change in the slope
    	which would be detected. This is by default 0.1, and typically need not be changed.

  		\item 'ellip.gate' if TRUE, it fits an ellipse on the data as a gate, otherwise the 
  		rectangle gating results are returned

      
  		\item 'scale' a value in [0,1) that scales the size of ellipse to fit if 'ellip.gate' is TRUE

  	}
  	}  
}
\value{
  A CellPopulation object
}

\author{
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
}
\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
lymph <- flowDensity(obj=f, channels=c('FSC-A', 'SSC-A'),
                     position=c(TRUE, FALSE), upper= c(NA, TRUE), debris.gate=c(TRUE, FALSE))
slotNames(lymph)
}

\seealso{
\code{\link{deGate}} \code{\link{getflowFrame}} \code{\link{notSubFrame}}
}

