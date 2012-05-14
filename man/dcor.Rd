\name{dcor}
\alias{dcor}
\title{
Distance correlation for tables
}
\description{
This computes the distance correlation for two categorical variables.
}
\usage{
dcor(x)
}
\arguments{
  \item{x}{
A data crosstable.
}
}

\value{
The correlation value which is between 0 and 1.
}
\references{
Szekely, G. J. Rizzo, M. L. and Bakirov, N. K. (2007). "Measuring and testing independence by correlation of distances", Annals of Statistics, 35/6, 2769-2794
}
\author{
Alexander Pilhoefer
}




\examples{
require(MASS)
p1 <- (p1<-runif(8, 0.5,1))/sum(p1)
p2 <- (p2<-runif(12,0.5,1))/sum(p2)

rmat <- matrix(rmultinom(1,1000,prob = as.vector(outer(p1,p2))), ncol=12)

dcor(rmat)
rmat2 <- optile(rmat)
dcor(rmat2 <- optile(rmat, fun = "distcor"))

fluctile(rmat)

bmat <- arsim(1000, c(12,12), 3)

dcor(bmat)
dcor(bmat2 <- optile(bmat, fun = "distcor"))

fluctile(bmat)
cfluctile(bmat2)

}

