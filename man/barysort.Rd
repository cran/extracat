\name{barysort}
\alias{barysort}
\title{
row and column moment reordering
}
\description{
An iterative  row and column reordering procedure based on the barycenter heuristic.
}
\usage{
barysort(x, vs = 1)
}

\arguments{
  \item{x}{
The data matrix.
}
  \item{vs}{
version. no effect.
}
}
\value{
The reordered matrix.
}
\examples{
# a good and quick solution:
a <- arsim(2000,c(14,14),3, noise=0.3)
fluctile(a2<-barysort(a))
BCI(a2)

# which is near
a3 <- optile(a, iter=100)
BCI(a3)

#but IBCC is even better:
a4 <- optile(a, fun="IBCC",iter=100)
BCI(a4)
}
