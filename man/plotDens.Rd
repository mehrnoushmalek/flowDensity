\name{plotDens}
\docType{methods}
\alias{plotDens}
\alias{plot, flowFrame, CellPopulation-method}
\title{
  Plot flow cytometry data with density-based colors
}

\description{
    Generate a scatter dot plot with colors based on the distribution of
    the density of the provided channels.
}
\usage{
plotDens(obj, channels,node=NA ,  col, main, xlab, ylab, pch=".", ...)
}
\arguments{
  \item{obj}{
  	a 'FlowFrame','GatingHierarchy' or 'cellPopulation' object.
  }
  \item{channels}{
    a vector of two channel names or their corresponding indices in the 'flow.frame'.
  }
  \item{node}{
   Only used if gatingHierarchy is provided as an input.It subtractes the flowFrame from the node of gatingHierarchy. see getNodes() in flowWorkspace.
  }
  
  \item{col}{
 	A specification for the default plotting color: see '?par'. 	
  }
  \item{main}{
  	an overall title for the plot: see '?plot'
  }
  \item{xlab}{
  	a title for the x axis: see '?plot'
  }
  \item{ylab}{
  	a title for the y axis: see '?plot'  	
  }
  \item{pch}{
  	Either an integer specifying a symbol or a single character to be used as the default in plotting points: see '?par'.
  }
  \item{...}{
    can be used to provide desired arguments for the plot() function used to
    plot the output results.
  }
}


\value{
  a scatter dot plot  with density-based colors.
}

\author{
Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
  
}

\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
#Plot CD3 vs. CD19 to see the distribution of cell populations and their density
plotDens(f,c("V450-A","PerCP-Cy5-5-A"))
}
\keyword{ Automated gating }

\keyword{ FlowCytData }
