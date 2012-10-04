\name{tfluctile}
\alias{tfluctile}
\title{
Fluctuation diagram with additional dendrograms
}
\description{
Plots a fluctuation diagram via \link[extracat]{fluctile} and adds dendrograms for the rows and columns to it.
}
\usage{
tfluctile(x, tree = NULL, dims = c(1, 2), tw = 0.2, border = NULL, tile.col = hsv(0.1, 0.1, 0.1, alpha = 0.6), bg.col = NULL, vp = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The two-way table or matrix with the data.
}
  \item{tree}{
A list with tree objects. This may be NULL or will be disregarded if \code{x} has an attribute \code{attr(x, "tree")} which should also be a list. The latter way is the standard for objects returned by \code{optile.list} or \code{optile}.
}
  \item{dims}{
If \code{x} has more than two dimensions this vector of length 2 indicates which variables to plot.
}
  \item{tw}{
The proportion of the total space to the left and at the top which is used for the dendrogram.
}
  \item{border}{
How much space is left white around the dendrogram.
}
  \item{tile.col}{
The tile color.
}
  \item{bg.col}{
A background color for the cells.
}
  \item{vp}{
A viewport to which the plot should be added or \code{NULL}.
}
  \item{\dots}{
further args
}
}
\value{
\code{invisible(TRUE)}
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{fluctile}, \link[extracat]{cfluctile}
}
\examples{
library(amap)
hc1 <- hcluster(t(plants[,-1]), method="manhattan", link = "ward")
hc2 <- hcluster(t(plants[,-1]), method="manhattan", link = "complete")

hclist <- list(hc1, hc2)
tfluctile( tt<-optile(hclist, k= c(8,8) ) )

s1 <- subtree(hc1, k = 12)
s2 <- subtree(hc2, k = 10)

tfluctile( table(s1$data, s2$data), tree = list(s1,s2))

}

