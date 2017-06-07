\name{flowDensity-methods}
\docType{methods}
\alias{flowDensity}
\alias{flowDensity,CellPopulation,ANY,logical,missing-method}
\alias{flowDensity,flowFrame,ANY,logical,missing-method}
\alias{flowDensity,GatingHierarchy,ANY,logical,ANY}
\alias{flowDensity-methods}
\title{ Methods for Function \code{flowDensity} in Package \pkg{flowDensity} }


\usage{
flowDensity(obj, channels, position,node, ...)
\S4method{flowDensity}{CellPopulation, ANY, logical, missing}(obj, channels, position, ...)
\S4method{flowDensity}{flowFrame,ANY, logical, missing)}(obj, channels, position, ...)
\S4method{flowDensity}{GatingHierarchy, ANY, logical, ANY}(obj, channels, position,node, ...)
}

\arguments{
\item{obj}{GatingHierarchy or \code{\linkS4class{CellPopulation}}object
}

 \item{channels}{
    a vector of two channel names or their corresponding indices.
  }
  \item{position}{
    a vector of two logical values specifying the position of the cell subset of interest on the 2D plot.
  }
  \item{node}{
   A character defining the parent population when extracting data from GatingHierarchy. 
  }
\item{...}{
 Check \code{flowDensity} description.
}
}
\value{
a CellPopulation object.
}
\description{
  Methods for function \code{flowDensity} in package \pkg{flowDensity} 
}


