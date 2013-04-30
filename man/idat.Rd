\name{idat}
\alias{idat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
indicator dataframe
}
\description{
Converts all factor variables in a dataframe to a set of binary variables.}
\usage{
idat(x, allcat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
dataframe
}
  \item{allcat}{
leave last category out
}
}

\value{
dataframe
}


\seealso{
\link[extracat]{Burt}, \link[extracat]{imat}
}
\examples{
idat(housing)
}

