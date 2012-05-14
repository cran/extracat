\name{fluctile3d}
\alias{fluctile3d}

\title{
3d fluctuation diagram
}
\description{
Uses rgl to create a 3-diamensional fluctuation diagram.}
\usage{
fluctile3d(x, shape = "cube", col = "darkgrey", alpha = 0.8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 3-dimensional data table or array.
}
\item{shape}{
	Either \code{"cube"} for cubes or \code{"oct"} for octahedrons.
}
  \item{col}{
The color of the cubes.
}
  \item{alpha}{
The alpha-value of the cubes.
}
  \item{\dots}{
dots
}
}

\value{
TRUE
}

\author{
Alexander Pilhoefer
}



\seealso{
\link[extracat]{fluctile}
}
\examples{
\dontrun{
A <- arsim(1000,c(5,7,9),3)
fluctile3d(A)
}
}

