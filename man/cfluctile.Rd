\name{cfluctile}
\alias{cfluctile}
\title{
cluster fluctuation for 2-way tables after optile
 }
\description{
Plots a fluctuation diagram of a 2-way table and adds red rectangles for clusters. Uses Kendalls tau for the cluster identification.
}
\usage{
cfluctile(x, tau0 = NULL, col = "red", gap.prop = 0.2, floor = 0, rev.y = FALSE, ...)
}
\arguments{
  \item{x}{
a 2-way table
}
  \item{tau0}{
the minimum value of kendalls tau as a stop criterion for each split.
}
  \item{col}{
color of the rectangles which highlight the clusters
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
M <- as.table(bloma(4,3,1.5))
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

M <- as.table(bloma(4,3,1.5))
MX <- optile(M, iter = 100)
MX <- ceiling(M/3) + MX
cfluctile(MX)

}
