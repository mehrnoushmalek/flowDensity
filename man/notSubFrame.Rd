\name{notSubFrame}
\docType{methods}
\alias{notSubFrame}
\title{
  Removing a subset of a FlowFrame object
}

\description{
    Remove a subset of a FlowFrame object specified by gates from the
    flowDensity method. It comes in handy when one needs the complement
    of a cell population in the input flow cytometry data.
}
\usage{
notSubFrame(flow.frame, channels, position = NA, gates, filter)
}
\arguments{
  \item{flow.frame}{
  	a 'FlowFrame' object.
  }
  \item{channels}{
    a vector of two channel names or their corresponding indices in the 'flow.frame'.
  }
  \item{position}{
    a vector of two logical values specifying the position of the cell subset of interest on the 2D plot.
  }
  \item{gates}{
    the gates slot in the CellPoulation object which is output by
    flowDensity function. It can also be a vector of two integer values
    each of which specifies a threshold for the corresponding channel in
    'channels' argument.
  }
  \item{filter}{
  	boundary of the subset to be removed. This value is stored in the 'filter' slot of a 'CellPopulation' object.
  }
}
\value{
  a CellPopulation object.
}

\author{
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
}
\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
#Find the threshold for CD20
cd20.gate <- deGate(f,channel="APC-H7-A")
# Gate out the CD20- populations using the notSubFrame
CD20.pos <- notSubFrame(f,channels=c("APC-H7-A","PerCP-Cy5-5-A"),position=c(FALSE,NA),gates=c(cd20.gate,NA))
#Plot the CD20+ cells on same channels
plotDens(CD20.pos@flow.frame,c("APC-H7-A","PerCP-Cy5-5-A"))
}
