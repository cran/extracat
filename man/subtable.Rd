\name{subtable}
\Rdversion{1.1}
\alias{subtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
data.frame reduction
}
\description{
Reduces a dataframe into a frequency table with prespecified entries.
}
\usage{subtable(dset, cols, keep.zero = TRUE, allfactor = TRUE, freqvar = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dset}{
The \code{data.frame} to reduce.
}
  \item{cols}{
An ordered integer vector containing the indices of the colums to keep.
}
  \item{keep.zero}{
A logical indicating whether to include zero-cases in the output.
}
  \item{allfactor}{
A logical indicating whether to convert all variables into factor variables. Integer variables will be applied a fitting (non-lexicographic) level order.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
hs2 = subtable(housing,c(3,1))
summary(hs2)
}
