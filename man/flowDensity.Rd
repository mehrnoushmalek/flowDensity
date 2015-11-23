\name{flowDensity}
\docType{methods}
\alias{flowDensity}
\alias{.flowDensity}
\title{
  Automated Sequential Gating Tool for Flow Cytometry
}
\description{
  flowDensity is an automated clustering algorithm which aims to emulate the current practice of manual sequential gating.
  It is designed to identify the predefined cell subsets based on the density distribution of the parent cell population by analyzing the peaks of the density curve.
}
\usage{
	flowDensity(obj, channels, position, singlet.gate, ...)
}
\arguments{
  \item{obj}{
  	a 'CellPopulation' or 'flowFrame' object.
  }
  \item{channels}{
    a vector of two channel names or their corresponding indices.
  }
  \item{position}{
    a vector of two logical values specifying the position of the cell subset of interest on the 2D plot.
  }
  \item{singlet.gate}{
  	if TRUE, singlet cell gate is derived and the corresponding cell population is returned. 
  }
  \item{...}{
  	This can be used to pass one of the following arguments:
  	\itemize{
  		\item 'use.percentile' if TRUE, returns the 'percentile'th threshold.

  		\item 'percentile' a value in [0,1] that is used as the percentile if 'use.percentile' is TRUE.

  		\item 'upper' if 'TRUE', it finds the change in the slope after the peak with index 'peak.ind'.

  		\item 'sd.threshold' if TRUE, it uses 'n.sd' times standard deviation for gating.

  		\item 'n.sd' an integer that is multiplied to the standard deviation to determine the place 
  		of threshold if 'sd.threshold' is 'TRUE'.

  		\item 'use.control' if TRUE, it finds the threshold using a matched control population and uses it for gating.

  		\item 'control' a 'flowFrame' or 'CellPopulation' object used for calculating the gating threshold when 
      'use.control' is set to TRUE. If a control population is used, the other arguments ('upper', 'percentile', 
      etc.) are applied to the control data when finding the threshold (i.e. not to 'obj').

  		\item 'alpha' a value in [0,1) specifying the significance of change in the slope
    	which would be detected. This is by default 0.1, and typically need not be changed.

  		\item 'debris.gate' if TRUE, it would try to remove the debris and gate the lymphocyte cells.
	
  		\item 'avg' if TRUE, it uses the mean of all identified cutoff points as a threshold for gating.
  		
  		\item 'ellip.gate' if TRUE, it fits an ellipse on the data as a gate, otherwise the 
  		rectangle gating results are returned
  		
      \item 'remove.neg' if TRUE, it removes negative values when removing margin events of scatter channels.
      
  		\item 'scale' a value in [0,1) that scales the size of ellipse to fit if 'ellip.gate' is TRUE
		
		\item 'graphs' if TRUE, the ellipse is added to the current plot of the output of the 'flowDensity',
		 otherwise a new plot is drawn and the ellipse is added on that.
  	}
  	}  
}
\value{
  A CellPopulation object
}

\author{
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
}
\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
lymph <- flowDensity(obj=f, channels=c('FSC-A', 'SSC-A'),
                     position=c(TRUE, FALSE), upper= c(NA, TRUE), debris.gate=c(TRUE, FALSE))
slotNames(lymph)
}

