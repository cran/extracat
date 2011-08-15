\name{spread}
\alias{spread}
\title{
  Matrix Expansion
}
\description{
Replaces each element of a matrix with a submatrix (containing this element) with a prespecified number of rows and columns.
}
\usage{spread(M, ncol = 1, nrow = 1)}
\arguments{
  \item{M}{The matrix to be expanded.}
  \item{ncol}{The number of column repetitions.}
  \item{nrow}{The number of row repetitions.}
  }
\details{
    
  }
\value{
  A matrix.
}
\references{
  Alexander Pilhoefer \emph{ New approaches in visualization of categorical data:
       R-package extracat}\cr
  Journal of Statistical Software, submitted March 2010
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{
This function will be advanced in the future to support more than two dimensions.
}
\seealso{
  
}
\examples{
 M = matrix(1:12, ncol=4)
 spread(M, ncol = 5, nrow = 7)
 }

