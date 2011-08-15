\name{bloma}
\alias{bloma}
\title{
random block matrices for testing}
\description{
Creates a random 2-way table with k blocks of average size ncr.
}
\usage{
bloma(k, ncr = 12, ar = 1)
}
\arguments{
  \item{k}{
Number of blocks
}
  \item{ncr}{
expected size of the blocks.
}
  \item{ar}{
ratio for the blocks
}
}

\value{
A matrix.
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



\examples{
bloma(3,3,2)
}
