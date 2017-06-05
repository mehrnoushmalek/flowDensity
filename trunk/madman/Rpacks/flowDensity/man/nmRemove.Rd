\name{nmRemove}
\docType{methods}
\alias{nmRemove}
\title{
	Preprocessing helper function for flow cytometry data
}
\description{
    Remove the margin events on the axes. Usually, these events are considered as debris or artifacts. This is specifically useful for 'FSC' and 'SSC' channels in a 'FlowFrame' object. However, any channel can be input as an argument. 
}
\usage{
	nmRemove(  flow.frame, channels, neg=FALSE, verbose=FALSE,return.ind=FALSE)
}
\arguments{
  \item{flow.frame}{
  	a 'FlowFrame' object.
  }
  \item{channels}{
  	a vector of channel names or their corresponding indices.
  }
  \item{neg}{
  	if TRUE, negative events are also removed
  }
\item{verbose}{
 if TRUE, it prints the margin event in each channel
}
\item{return.ind}{
 if TRUE, it return indices of margin events for each channel.
}
}
\value{
  a 'FlowFrame' object, or a 'list' of indices identifying margin events for each channel.
}


\author{
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
  Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
}

\examples{

data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_2', data_dir, full = TRUE))
#Removing margin events of FSC-A and SSC-A channels
no.margin <- nmRemove(f2, c("FSC-A","SSC-A"),talk=TRUE)
plotDens(f2, c("FSC-A","SSC-A"))
# Scatter plot of FSC-A vs. SSC-A after removing margins
plotDens(no.margin, c("FSC-A","SSC-A"))
}
