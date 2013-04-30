\name{boxplot2g}
\alias{boxplot2g}
\title{
2D projection boxplots
}
\description{
Contours derived from the boxplot marks of data projections.
}
\usage{
boxplot2g(x, y = NULL, groups = NULL, smooth = loess,
smooth.args = list(span = 0.1), colv = NULL, alpha = 1, n = 360, ...)
}
\arguments{
  \item{x}{
Numeric x variable.
}
  \item{y}{
Numeric y variable.
}
  \item{groups}{
A grouping variable or \code{NULL}.
}
  \item{smooth}{
A smoothing function such as \code{loess} or \code{FALSE}.
}
  \item{smooth.args}{
Additional arguments for the smoother, such as span.
}
  \item{colv}{
A color vector with one color for each group.
}
  \item{alpha}{
The alpha-blending value for the points of the scatterplot.
}
  \item{n}{
The number of irection vectors. Defaults to 360 which means 1-degree steps.
}
  \item{\dots}{
Further args.
}
}
\details{
The data is first centered and standardized via multiplication with the inverse of the covariance matrix. Then projections of the data onto different direction vectors (0 to 360 degree) are computed. The boxplot marks (median, box, whisker) are computed for each such projection and finally retransformed
again using the covariance matrix. Optionally the projected marks are then smoothed, e.g. via \code{smooth=loess}.
}
\value{
The ggplot object.
}

\author{
Alexander Pilhoefer
}
\note{
Please note that I do not claim that the resulting contours have any statistically useful characteristics! Usage at one's own risk!
}


\seealso{
\link[graphics]{contour}
}
\examples{
\dontrun{
data(olives)
boxplot2g(olives$palmitoleic,olives$oleic,olives$Region)
boxplot2g(olives$palmitoleic,olives$oleic,olives$Area,
alpha=0.5,smooth.args=list(span=0.3), colv = rainbow(9))
}
}

