\name{CellPopulation-class}
\Rdversion{1.1}
\docType{class}
\alias{CellPopulation-class}
\alias{getflowFrame,CellPopulation-method}
\alias{plot,flowFrame,CellPopulation-method}

\title{Class \code{"CellPopulation"}}
\description{
	This class represents the output of 'flowDensity(.)' function from flowDensity package.
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("CellPopulation", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{flow.frame}:}{
    	Object of class \code{"flowFrame"} representing the flow cytometry data of the cell population
    	}
    \item{\code{proportion}:}{
    	Object of class \code{"numeric"} representing proportion of the cell population with respect to its parent cell population
    	}
    \item{\code{cell.count}:}{
    	Object of class \code{"numeric"} representing cell count of the cell population
    	}
    \item{\code{channels}:}{
    	Object of class \code{"character"} representing channel names corresponding to the 2 dimensions where the cell population is extracted
    	}
    \item{\code{position}:}{
    	Object of class \code{"logical"} representing position of the cell population in the 2-dimensional space
    	}
    \item{\code{gates}:}{
    	Object of class \code{"numeric"} representing thresholds on each channel used to gate the cell population
    	}
    \item{\code{filter}:}{
    	Object of class \code{"matrix"} representing boundary of the cell population using a convex polygon
    	}
    \item{\code{index}:}{
    	Object of class \code{"numeric"} representing indices of the data points in the cell population with respect to its parent cell population
    	 }
  }
}
\section{Methods}{
  \describe{
    \item{flowDensity}{\code{signature(obj = "CellPopulation", channels = "ANY", position = "logical", singlet.gate = "missing")}: ... }
    \item{flowDensity}{\code{signature(obj = "CellPopulation", channels = "missing", position = "missing", singlet.gate = "logical")}: ... }
    \item{getflowFrame}{\code{signature(obj = "CellPopulation")}: ... }
    \item{plot}{\code{signature(x = "flowFrame", y = "CellPopulation")}: ... }
	 }
}
\author{
	Jafar Taghiyar <email: <jtaghiyar@bccrc.ca>> 
}
\examples{
showClass("CellPopulation")
}
\keyword{classes}
