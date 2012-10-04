\name{BCImat}
\alias{BCImat}

\title{
BCI matrix
}
\description{
Computes pairwise BCI values via \link[extracat]{qBCI}.
}
\usage{
BCImat(x, k = 5, iter = 20, p = NULL)
}

\arguments{
  \item{x}{
A dataframe.
}
  \item{k}{
The minimum expected number for each cell after quantile binning. See \link[extracat]{qBCI}.
}
  \item{iter}{
An \link[extracat]{optile} parameter.
}
  \item{p}{
The quantile distance. See \link[extracat]{qBCI}.
}
}
\details{
Uses pairwise complete cases only!
}
\value{
A symmetric matrix with BCI values between 0 and 1.
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{qBCI}, See \link[extracat]{wdcor}.
}
\examples{
m1 <- BCImat(olives)
fluctile(1 - m1,shape="o")
}
