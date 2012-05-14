\name{arsim}
\alias{arsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
block-structured arrays}
\description{
Generates an array or matrix that includes k fully separated block-clusters.
}
\usage{
arsim(n, dim, k, noise = 0, shuffle = TRUE, v = 0.1, minc = 1, exp.prop = NULL, min.prop = 1/dim/4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of observations in the array.
}
  \item{dim}{
The dimension of the array.
}
  \item{k}{
The number of clusters. 1 for no clusters.
}
  \item{noise}{
The proportion of the \code{n} observations which is white noise.
}
  \item{shuffle}{
Whether or not to disrupt the block structure by randomly mixing up the category orders.
}
  \item{v}{
The variability with which the observations are assigned to the clusters.
}
  \item{minc}{
The minimum number of categories each cluster must have in each variable. E.g. \code{minc = 2} means, that each cluster is at least 2 categories thick in all dimensions.
}
  \item{exp.prop}{
expected proportions of observations in the clusters.
}
  \item{min.prop}{
minimum proportion of observations in each cluster.
}
}
\details{
Not a very sophisticated way of generating random arrays but it serves for tests and illustrations of the other functions.
}
\value{
A simulated data array.
}
\examples{
A <- arsim(1000, c(12,12), 3, shuffle = FALSE)
fluctile(A)

\dontrun{
A2<- arsim(1000, c(12,12,12), 3, shuffle = FALSE)
fluctile3d(A2, shape ="oct")
}

}

