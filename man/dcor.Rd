\name{wdcor}
\alias{wdcor}
\alias{wdcor.table}
\alias{wdcor.default}
\alias{wdcor.data.frame}

\title{
Weighted Distance Correlation
}
\description{
An efficient implementation of the distance correlation for two variables with the additional option to weight the observations.
The main application for the weights is to use frequencies according to ordinal variables which can be represented by a contingency table (\code{wdcor.table)}).
 Another idea is to make the distance correlation more robust by assigning small weights to observations which are far from the rest of the data.
 
 For large datasets the distance correlation is often said to be too inefficient to be of any great use. The function \link[extracat]{approx.dcor} offers a pretty good approximation of the distance correlation via binning and \code{wdcor.table)}.
 
\code{wdcor.data.frame} computes a distance correlation matrix. Factor variables are transformed to integer via \code{data.matrix}.
}
\usage{
wdcor(x,\dots)
\method{wdcor}{default}(x,y,w = NULL,ep = 1,\dots)
\method{wdcor}{table}(x,ep = 1,\dots)
\method{wdcor}{data.frame}(x,approx = TRUE,\dots)
}

\arguments{
  \item{x}{
A data crosstable or a numeric vector.
}
  \item{y}{
A numeric vector.
}
\item{w}{
Weights, typically frequencies. The default weights all cases the same, which leads to the standard distance correlation.
}
  \item{ep}{
The euclidean (absolute) distances can be taken to the power of \code{ep}.
}
\item{approx}{
Whether or not to use \link[extracat]{approx.dcor} instead of \link[extracat]{wdcor} for the \code{data.frame} method.
}
\item{\dots}{
dots.
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

\seealso{\link[extracat]{approx.dcor}}


\examples{

# repeat and change N for different results and computation times.
N <- 2000
x1 <- rnorm(N,mean=10,sd=3)
x2 <- runif(N,0,40)
x3 <- rnorm(N,mean=30,sd=4)
x <- sample(c(x1,x2,x3),N)		
			
y <- rnorm(1,sd=0.0001)*(x-mean(x))^4+ rnorm(1,sd=0.01)*(x-mean(x))^3 + rnorm(1,sd=0.1)*(x-mean(x))^2 + rnorm(1)*(x-mean(x))+rnorm(N,sd=runif(N,3,10))+ runif(N,0,20)*sin(abs(scale(x))*2*pi)

plot(x,y,pch=19,col=alpha("black",0.2))
system.time(dd<-wdcor(x,y))

y2 <- runif(2000)
system.time(dde<-wdcor(x,y2))
dd
dde

}

