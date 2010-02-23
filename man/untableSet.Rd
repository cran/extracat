\name{untableSet}
\Rdversion{1.1}
\alias{untableSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
data.frame conversion
}
\description{
Converts a frequency table into a raw data.frame.
}
\usage{untableSet(dset, freqvar = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dset}{
The \code{data.frame} including a frequency variable "Freq".
}

\item{freqvar}{
    Optional name of a frequency variable in \code{V}. 
    If \code{dset} contains a variable called \code{"Freq"} (see \link[stats]{ftable}) 
    it will be defined as frequency variable if \code{freqvar} is unspecified.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{data.frame}.
}
\references{
  Alexander Pilhoefer \emph{ New approaches in visualization of categorical data:
       R-package extracat}\cr
  Journal of Statistical Software, submitted March 2010
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
hs2 = untableSet(housing)
summary(hs2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ frequency table }

