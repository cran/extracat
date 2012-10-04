\name{cfluctile}
\alias{cfluctile}
\title{
cluster fluctuation for 2-way tables after optile
 }
\description{
Plots a fluctuation diagram of a 2-way table and adds red rectangles for clusters. Uses Kendalls tau for the cluster identification.
}
\usage{
cfluctile(x, tau0 = NULL, method = "Kendall", col = "red", lwd =
                 2, lty = 1, gap.prop = 0.2, floor = 0, rev.y = FALSE,
                 add = FALSE, shape = "r", just = "c", dir = "b", ...)
}
\arguments{
  \item{x}{
A 2-way table
}
  \item{tau0}{
The minimum value of Kendalls tau or Cohens Kappa as a stop criterion for each split.
}
\item{method}{
Either \code{"Kendall"} for Kendalls Tau or \code{"Cohen"} for Cohens Kappa.
}
  \item{col}{
color of the rectangles which highlight the clusters
}
\item{lwd}{
	The line width.
}
\item{lty}{
The line type.
}

  \item{gap.prop}{
proportion of the gaps between rows/columns.
}
  \item{floor}{
floor censored zooming: all cases will be plotted but only those with a frequency of at least \code{floor} will be considered for the clustering.
}
  \item{rev.y}{
revert the y axis.
}
\item{add}{
	Whether to make a new plot or to add to an existing one.
}
\item{shape}{
	The shape of the objects. See \link[extracat]{fluctile}.
}
\item{just}{
	See \link[extracat]{fluctile}.
}
\item{dir}{
	See \link[extracat]{fluctile}.
}
  \item{\dots}{
dots}
}
\details{
This function calls \link[extracat]{fluctile} to create a 2-way fluctuation diagram and then adds cluster rectangles to it. The cluster rectangles are computed in the following way: 

First Kendalls tau is computed for the whole table if no tau0 value is defined.

Initially all columns and all rows are seen as one single cluster (1x1 table).

Then the columns as well as the rows are split into 2 groups each such that Kendalls tau is maximal for the resulting table (2x2).

This process is then repeated within the resulting diagonal blocks as long as the tau values are better than tau0.
}
\value{
\code{invisible(TRUE)}
}

\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{
This was part of the Google Summer of Code 2011.
}


\seealso{
\link[vcd]{mosaic}
}

\examples{
M <- arsim(1000,c(10,10),3)
M2 <- optile(M)
cfluctile(M2)

M2 <- optile(M, iter = 100)
cfluctile(M2)

M2 <- optile(M, iter = 100, presort = TRUE)
cfluctile(M2)

\dontrun{
M2 <- optile(M, fun = "hamm", iter = 20)
cfluctile(M2)
}

M2 <- optile(M, fun = "casort")
cfluctile(M2)

M2 <- optile(M, fun = "csvd")
cfluctile(M2)

M2 <- optile(M, fun = "preclass", foreign = ".Call")
cfluctile(M2)

M <- arsim(1000,c(10,10),3)
MX <- optile(M, iter = 100)
MX <- ceiling(M/3) + MX
cfluctile(MX)

}
