\name{fluctile3d}
\alias{fluctile3d}

\title{
3d fluctuation diagram
}
\description{
Uses rgl to create a 3-dimensional fluctuation diagram.}
\usage{
fluctile3d(x, shape = "cube", col = "darkgrey", col.array = NULL,
alpha = 0.8, add = FALSE, lab = TRUE,...)
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
  \item{col.array}{
An optional array with colors for the cubes.
}
  \item{alpha}{
The alpha-value of the cubes.
}
  \item{add}{
whether to create a new plot or to add the cases to an existing one.
}
  \item{lab}{
whether or not to draw labels
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

fluctile3d(A,shape="oct")

fluctile3d(A,shape="oct", col.array = f3dcol(A,dims=2,col.fun=rainbow_hcl))
}
}

