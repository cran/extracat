\name{f3dcol}
\alias{f3dcol}

\title{
Color array for fluctile3d
}
\description{
Returns an array with color values with the same dimensions as the table or array x.
Written for the \code{col.array} option in \link[extracat]{fluctile3d}.
}
\usage{
f3dcol(x, dims = c(1, 2), col.fun = rainbow_hcl, col.opt = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 3D data table or array.
}
  \item{dims}{
Which dimensions get differing color values. E.g. \code{dims = c(1,3)} means that the color values are the same for all cases with identical values in dimension 1 and 3 but possibly different values in dimension 2. \code{dims = c(1,2,3)} leads to a different color for each entry and if \code{dims = 1} only the first dimension is relevant for the color choice.  
}
  \item{col.fun}{
A function which produces a palette of \code{n} colors like \link[grDevices]{rainbow} or \link[colorspace]{sequential_hcl}, \link[colorspace]{rainbow_hcl}, \link[colorspace]{diverge_hcl}.
}
  \item{col.opt}{
Options passed to the color palette function. NOT YET IMPLEMENTED.
}
}

\value{
An array with hex colors.
}

\author{
Alexander Pilhoefer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{fluctile3d}
}
\examples{
\dontrun{
x <- arsim(2000,c(5,7,4),2,shuffle=FALSE)
fluctile3d(x,col.array = f3dcol(x,col.fun = rainbow))
fluctile3d(x,col.array = f3dcol(x,col.fun = function(z) sample(rainbow(z))))
fluctile3d(x,col.array = f3dcol(x,col.fun = diverge_hcl))
}
}

