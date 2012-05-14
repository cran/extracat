\name{classcrit}
\alias{classcrit}
\alias{hammcrit}
\alias{ihcrit}
\alias{iccrit}
\alias{shcrit}
\alias{sccrit}

\title{
The classification criterion
}
\description{
These functions compute the classification criterion, the weighted classification criterion and their scaled versions.
}
\usage{
classcrit(x, concordant = TRUE)
hammcrit(x)
ihcrit(x)
iccrit(x)
shcrit(x)
sccrit(x, concordant = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix, table or array. The hamming criterion is currently only implemented for two-dimensional tables and arrays.
}
  \item{concordant}{
Whether to return the number of concordant or discordant observation pairs.
}
}
\details{
\code{classcrit} and \code{hammcrit} compute the criteria itself. \code{iccrit} and  \code{ihcrit} compute the worst possible value after optimization for the given marginal distributions.  \code{sccrit} and  \code{shcrit} are equivalent to  \code{classcrit(x, concordant = FALSE))/iccrit(x)} and  \code{hammcrit(x)/ihcrit(x)} respectively.
}
\value{
The criterion value.
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{kendalls}
}
\examples{
M <-arsim(1000, c(12,12), 3)
classcrit(M)
hammcrit(M)
shcrit(M)
sccrit(M)

M2 <- optile(M, iter = 100)

classcrit(M2)
hammcrit(M2)
shcrit(M2)
sccrit(M2)
}

