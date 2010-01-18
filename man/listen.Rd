\name{listen}
\Rdversion{1.1}
\alias{listen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
auto-update interactivity for cpcp plots
}
\description{
Makes the \code{cpcp} plot listen to selection changes and calls \link{resort} to update the display.
}
\usage{listen()}
%- maybe also 'usage' for other objects documented here.
\details{
Choose BREAK in the plot menu in order to get back to the console.
}
\value{
No return value.
}
\references{
  Alexander Pilhoefer \emph{ New approaches in visualization of categorical data:
       R-package extracat}\cr
  Journal of Statistical Software, submitted Jan 2010
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[iplots]{ievent.wait}
\link[extracat]{cpcp}
\link[extracat]{resort}
}
\examples{
cpcp(housing,ord = c(4,3,2,1),gap.space=0.15)
## Not run:
listen()
## End(Not run)
# now make your selections and see what happens.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ event }
\keyword{ listen }% __ONLY ONE__ keyword per line
