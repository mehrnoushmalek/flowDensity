\name{getflowFrame}
\docType{methods}
\alias{getflowFrame}
\alias{getflowFrame, CellPopulation-method}
\title{
	'CellPopulation' class accessor.
}
\description{
    an accessor for 'CellPopulation' class to get its 'FlowFrame' object. This will remove all the NA values in the frame.
}
\usage{
	getflowFrame(obj)
}
\arguments{
  \item{obj}{
  	a 'CellPopulation' object.
  }
}
\value{
  a 'FlowFrame' object.
}
\author{
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
}

\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
lymph <- flowDensity(obj=f, channels=c('FSC-A', 'SSC-A'),
                     position=c(TRUE, NA))
f.lymph <- getflowFrame(lymph)
}
