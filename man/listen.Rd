\name{listen}
\alias{listen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
auto-update interactivity for cpcp plots
}
\description{
Makes the \code{cpcp} plot listen to selection changes and calls \link{resort} to update the display.
}
\usage{listen(init = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{ logical. If \code{FALSE} no message will be printed. It is usually not necessary for the user to change this parameter.}
}
\details{
Choose BREAK in the File menu in order to return back to the console.
}
\value{
No return value.
}
\references{
Alexander Pilhoefer, Antony Unwin (2013). New Approaches in Visualization of Categorical Data: R Package
  extracat. Journal of Statistical Software, 53(7), 1-25. URL \url{http://www.jstatsoft.org/v53/i07/}
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}

\seealso{
\link[iplots]{ievent.wait}
\link[extracat]{cpcp}
\link[extracat]{resort}
}
\examples{
\dontrun{
cpcp(housing,ord = c(4,3,2,1),gap.space=0.15)
listen()
}
# now make your selections and see what happens.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ event }
\keyword{ listen }% __ONLY ONE__ keyword per line
