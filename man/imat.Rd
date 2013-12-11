\name{imat}
\alias{imat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
indicator variables
}
\description{
converts a single categorical variable into indicator variables
}
\usage{
imat(s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
factor variable
}
}

\value{
matrix
}

\examples{
require(MASS)
imat(housing$Type)
}

